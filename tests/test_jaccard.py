#!/usr/bin/env python
import importlib.util
import math
import os
import shutil
import sys
import tempfile
import types
import unittest
from unittest import mock

import numpy as np
import pandas as pd


def load_jaccard_score_module():
    """
    Loads bin/jaccard_score.py. pybedtools (and the bedtools binary it wraps)
    is not needed by these tests, so it is stubbed out when it cannot be
    imported.
    """
    try:
        import pybedtools  # noqa: F401
    except Exception:
        stub = types.ModuleType("pybedtools")
        stub.BedTool = mock.MagicMock(name="BedTool")
        sys.modules["pybedtools"] = stub
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    module_path = os.path.join(repo_root, "bin", "jaccard_score.py")
    spec = importlib.util.spec_from_file_location("jaccard_score", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


##########################################
# Reference implementation, independent of bedtools
def read_bed(path):
    """Reads a bed file into {chrom: [merged (start, end), ...]}"""
    intervals = {}
    with open(path, "r") as bed:
        for line in bed:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            intervals.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    merged = {}
    for chrom, chrom_intervals in intervals.items():
        for start, end in sorted(chrom_intervals):
            if merged.get(chrom) and start <= merged[chrom][-1][1]:
                last_start, last_end = merged[chrom][-1]
                merged[chrom][-1] = (last_start, max(last_end, end))
            else:
                merged.setdefault(chrom, []).append((start, end))
    return merged


def total_bp(intervals):
    return sum(end - start for chrom in intervals for start, end in intervals[chrom])


def reference_jaccard(fileA, fileB):
    """
    Computes what bedtools jaccard reports for two bed files: intersection bp,
    union bp, jaccard (intersection/union) and the number of intersecting
    interval pairs. An empty union gives an undefined (NaN) jaccard, exactly
    as bedtools does.
    """
    a, b = read_bed(fileA), read_bed(fileB)
    intersection, n_intersections = 0, 0
    for chrom in set(a) & set(b):
        for a_start, a_end in a[chrom]:
            for b_start, b_end in b[chrom]:
                overlap = min(a_end, b_end) - max(a_start, b_start)
                if overlap > 0:
                    intersection += overlap
                    n_intersections += 1
    union = total_bp(a) + total_bp(b) - intersection
    jaccard = float("nan") if union == 0 else intersection / union
    return intersection, union, jaccard, n_intersections


def fake_run_jaccard(fileA, fileB, genomefile):
    """
    Drop-in replacement for jaccard_score.run_jaccard that mirrors the key
    ordering and stringification of the bedtools-backed version.
    """
    intersection, union, jaccard, n_intersections = reference_jaccard(fileA, fileB)
    record = {
        "intersection": intersection,
        "union-intersection": union,
        "jaccard": jaccard,
        "n_intersections": n_intersections,
        "fileA": os.path.basename(fileA),
        "fileB": os.path.basename(fileB),
    }
    keylist = sorted(record.keys())
    return ([str(record[key]) for key in keylist], keylist)


class TestParseScore(unittest.TestCase):
    def setUp(self):
        self.mod = load_jaccard_score_module()

    def test_undefined_scores_stay_nan(self):
        for score in ("nan", "NaN", float("nan"), "", ".", None):
            self.assertTrue(math.isnan(self.mod.parse_score(score)))

    def test_numpy_and_pandas_nans_stay_nan(self):
        for score in (
            np.nan,
            np.float64("nan"),
            np.float32("nan"),
            np.float64(np.nan),
            pd.NA,
            pd.NaT,
            pd.Series([np.nan])[0],
            pd.array([None], dtype="Float64")[0],
        ):
            self.assertTrue(math.isnan(self.mod.parse_score(score)))

    def test_non_scalar_scores_stay_nan(self):
        """pd.isna on a non-scalar is ambiguous, it must not raise"""
        self.assertTrue(math.isnan(self.mod.parse_score(np.array([np.nan, 1.0]))))
        self.assertTrue(math.isnan(self.mod.parse_score([np.nan])))

    def test_real_scores_are_preserved(self):
        self.assertEqual(self.mod.parse_score("0.3333"), 0.3333)
        self.assertEqual(self.mod.parse_score(0.5), 0.5)
        self.assertEqual(self.mod.parse_score("0"), 0.0)
        self.assertEqual(self.mod.parse_score("1"), 1.0)
        # numpy/pandas scalars that are not missing values
        self.assertEqual(self.mod.parse_score(np.float64(0.25)), 0.25)
        self.assertEqual(self.mod.parse_score(np.float32(0.5)), 0.5)
        self.assertEqual(self.mod.parse_score(pd.Series([0.75])[0]), 0.75)


class TestRunJaccard(unittest.TestCase):
    def setUp(self):
        self.mod = load_jaccard_score_module()

    def test_bedtools_keys_match_the_table_columns(self):
        """
        The tabular output is written with the module's own column order, so
        the sorted bedtools keys have to line up with it.
        """
        bedtool = mock.MagicMock()
        bedtool.sort.return_value = bedtool
        bedtool.jaccard.return_value = {
            "intersection": 50,
            "union-intersection": 150,
            "jaccard": 0.333333,
            "n_intersections": 1,
        }
        with mock.patch.object(self.mod, "BedTool", return_value=bedtool):
            data, keylist = self.mod.run_jaccard("/tmp/a.bed", "/tmp/b.bed", "genome.txt")
        self.assertEqual(keylist, list(self.mod.TABLE_COLUMNS))
        self.assertEqual(float(data[keylist.index("jaccard")]), 0.333333)
        self.assertEqual(data[keylist.index("fileA")], "a.bed")
        self.assertEqual(data[keylist.index("fileB")], "b.bed")


class PeakFileTestCase(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.mod = load_jaccard_score_module()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def write_peaks(self, name, intervals):
        path = os.path.join(self.test_dir, name)
        with open(path, "w") as peaks:
            for chrom, start, end in intervals:
                peaks.write("\t".join([chrom, str(start), str(end)]) + "\n")
        return path

    def table_rows(self, outTable):
        """Parses the tabular output into a list of {column: value} rows"""
        header = outTable[0].split("\t")
        return [dict(zip(header, row.split("\t"))) for row in outTable[1:]]


class TestLoopJaccard(PeakFileTestCase):
    def setUp(self):
        super().setUp()
        self.patcher = mock.patch.object(self.mod, "run_jaccard", fake_run_jaccard)
        self.patcher.start()

    def tearDown(self):
        self.patcher.stop()
        super().tearDown()

    def test_accurate_jaccard_scores(self):
        # A: 100bp, B: 100bp overlapping A by 50bp, C: 100bp disjoint from both
        fileA = self.write_peaks("A_peaks.bed", [("chr1", 0, 100)])
        fileB = self.write_peaks("B_peaks.bed", [("chr1", 50, 150)])
        fileC = self.write_peaks("C_peaks.bed", [("chr1", 200, 300)])
        outTable, out, snames = self.mod.loop_jaccard([fileA, fileB, fileC], "genome.txt")
        self.assertEqual(snames, ["A", "B", "C"])
        # 50bp intersection over a 150bp union
        self.assertAlmostEqual(out.loc["A", "B"], 50 / 150)
        # no intersection at all
        self.assertAlmostEqual(out.loc["A", "C"], 0.0)
        self.assertAlmostEqual(out.loc["B", "C"], 0.0)
        # self comparisons and symmetry
        for sname in snames:
            self.assertEqual(out.loc[sname, sname], 1.0)
        for left in snames:
            for right in snames:
                self.assertEqual(out.loc[left, right], out.loc[right, left])
        # header line plus one line per pairwise comparison
        self.assertEqual(outTable[0].split("\t"), list(self.mod.TABLE_COLUMNS))
        rows = self.table_rows(outTable)
        self.assertEqual(len(rows), 3)
        self.assertAlmostEqual(float(rows[0]["jaccard"]), 50 / 150)
        self.assertEqual(rows[0]["intersection"], "50")
        self.assertEqual(rows[0]["union-intersection"], "150")

    def test_multi_interval_jaccard_scores(self):
        fileA = self.write_peaks(
            "A_peaks.bed", [("chr1", 0, 100), ("chr1", 200, 300), ("chr2", 0, 100)]
        )
        fileB = self.write_peaks("B_peaks.bed", [("chr1", 75, 250), ("chr3", 0, 100)])
        outTable, out, snames = self.mod.loop_jaccard([fileA, fileB], "genome.txt")
        # intersection: 25bp (chr1 75-100) + 50bp (chr1 200-250) = 75bp
        # union: 300bp + 275bp - 75bp = 500bp
        self.assertAlmostEqual(out.loc["A", "B"], 75 / 500)

    def test_samples_with_no_peaks_are_na(self):
        fileA = self.write_peaks("A_peaks.bed", [("chr1", 0, 100)])
        fileB = self.write_peaks("B_peaks.bed", [("chr1", 50, 150)])
        empty = self.write_peaks("C_peaks.bed", [])
        outTable, out, snames = self.mod.loop_jaccard(
            [fileA, fileB, empty], "genome.txt"
        )
        # the sample with no peaks is NA everywhere, itself included
        self.assertTrue(out.loc["A", "C"] != out.loc["A", "C"])  # NaN
        self.assertTrue(math.isnan(out.loc["C", "A"]))
        self.assertTrue(math.isnan(out.loc["C", "C"]))
        # the samples with peaks are untouched
        self.assertAlmostEqual(out.loc["A", "B"], 50 / 150)
        self.assertEqual(out.loc["A", "A"], 1.0)
        # every comparison is still reported, with NA scores
        rows = self.table_rows(outTable)
        self.assertEqual(len(rows), 3)
        na_rows = [row for row in rows if "C_peaks.bed" in (row["fileA"], row["fileB"])]
        self.assertEqual(len(na_rows), 2)
        for row in na_rows:
            for column in ("intersection", "jaccard", "n_intersections", "union-intersection"):
                self.assertEqual(row[column], "NA")

    def test_bedtools_is_not_run_on_files_without_peaks(self):
        """bedtools cannot score an empty peak file, so it is never asked to"""
        fileA = self.write_peaks("A_peaks.bed", [("chr1", 0, 100)])
        empty = self.write_peaks("B_peaks.bed", [])
        calls = []

        def recording_run_jaccard(fileA, fileB, genomefile):
            calls.append((fileA, fileB))
            return fake_run_jaccard(fileA, fileB, genomefile)

        with mock.patch.object(self.mod, "run_jaccard", recording_run_jaccard):
            self.mod.loop_jaccard([fileA, empty], "genome.txt")
        self.assertEqual(calls, [])

    def test_missing_peak_file_is_na(self):
        fileA = self.write_peaks("A_peaks.bed", [("chr1", 0, 100)])
        missing = os.path.join(self.test_dir, "B_peaks.bed")
        outTable, out, snames = self.mod.loop_jaccard([fileA, missing], "genome.txt")
        self.assertTrue(math.isnan(out.loc["A", "B"]))
        self.assertEqual(self.table_rows(outTable)[0]["jaccard"], "NA")

    def test_undefined_bedtools_score_is_na(self):
        """
        A NaN handed back by bedtools for two non-empty files is reported as
        NA too, whether it arrives as a string or as a numpy float
        """
        fileA = self.write_peaks("A_peaks.bed", [("chr1", 0, 100)])
        fileB = self.write_peaks("B_peaks.bed", [("chr1", 50, 150)])

        for nan_score in ("nan", np.str_(np.nan)):
            def nan_run_jaccard(fileA, fileB, genomefile, nan_score=nan_score):
                data, keylist = fake_run_jaccard(fileA, fileB, genomefile)
                data[keylist.index("jaccard")] = nan_score
                return (data, keylist)

            with mock.patch.object(self.mod, "run_jaccard", nan_run_jaccard):
                outTable, out, snames = self.mod.loop_jaccard(
                    [fileA, fileB], "genome.txt"
                )
            self.assertTrue(math.isnan(out.loc["A", "B"]))
            self.assertEqual(self.table_rows(outTable)[0]["jaccard"], "NA")

    def test_single_file_matrix(self):
        fileA = self.write_peaks("A_peaks.bed", [("chr1", 0, 100)])
        outTable, out, snames = self.mod.loop_jaccard([fileA], "genome.txt")
        self.assertEqual(out.shape, (1, 1))
        self.assertEqual(out.loc["A", "A"], 1.0)
        # header only, there is nothing to compare
        self.assertEqual(outTable, ["\t".join(self.mod.TABLE_COLUMNS)])


class TestDropNaSamples(unittest.TestCase):
    def setUp(self):
        self.mod = load_jaccard_score_module()

    def matrix(self):
        nan = float("nan")
        return pd.DataFrame(
            [
                [1.0, 0.5, nan],
                [0.5, 1.0, nan],
                [nan, nan, nan],
            ],
            columns=["A", "B", "C"],
            index=["A", "B", "C"],
        )

    def test_na_samples_are_dropped(self):
        out, snames = self.mod.drop_na_samples(self.matrix(), ["A", "B", "C"])
        self.assertEqual(list(out.columns), ["A", "B"])
        self.assertEqual(list(out.index), ["A", "B"])
        self.assertEqual(snames, ["A", "B"])
        self.assertFalse(out.isna().any().any())

    def test_complete_matrix_is_untouched(self):
        full = self.matrix().loc[["A", "B"], ["A", "B"]]
        out, snames = self.mod.drop_na_samples(full, ["A", "B"])
        self.assertEqual(list(out.columns), ["A", "B"])
        self.assertEqual(snames, ["A", "B"])

    def test_all_na_matrix(self):
        nan = float("nan")
        all_na = pd.DataFrame(
            [[nan, nan], [nan, nan]], columns=["A", "B"], index=["A", "B"]
        )
        out, snames = self.mod.drop_na_samples(all_na, ["A", "B"])
        self.assertEqual(out.shape, (0, 0))
        self.assertEqual(snames, [])


class TestPlotsExcludeNaSamples(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.mod = load_jaccard_score_module()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def path(self, name):
        return os.path.join(self.test_dir, name)

    def matrix(self, nvalid=2):
        """Square matrix of nvalid scored samples plus one NA sample"""
        snames = ["S%d" % i for i in range(1, nvalid + 1)] + ["NAsample"]
        out = pd.DataFrame(
            float("nan"), columns=snames, index=snames, dtype="float"
        )
        for left in snames[:nvalid]:
            for right in snames[:nvalid]:
                out.loc[left, right] = 1.0 if left == right else 0.4
        return (out, snames)

    def test_pca_tab_holds_only_scored_samples(self):
        out, snames = self.matrix()
        pcatab, pcaplot = self.path("pca.tsv"), self.path("pca.pdf")
        self.mod.pca_plot(out, snames, "macsNarrow", pcatab, pcaplot)
        PCAdata = pd.read_csv(pcatab, sep="\t")
        self.assertEqual(sorted(PCAdata["sample_name"]), ["S1", "S2"])
        self.assertNotIn("NAsample", list(PCAdata["sample_name"]))
        self.assertFalse(PCAdata[["PC1", "PC2"]].isna().any().any())
        self.assertTrue(os.path.getsize(pcaplot) > 0)

    def test_heatmap_tab_reports_na_samples(self):
        out, snames = self.matrix()
        heatmap_tab, heatmap = self.path("hm.tsv"), self.path("hm.pdf")
        self.mod.plot_heatmap(out, heatmap, "macsNarrow", heatmap_tab, snames)
        # the table keeps every sample, so its columns stay aligned across
        # peak callers, with NA for the samples that have no peaks
        with open(heatmap_tab, "r") as tab:
            lines = tab.read().strip().split("\n")
        self.assertEqual(lines[0].split("\t"), ["S1", "S2", "NAsample", "peakcaller"])
        self.assertEqual(lines[-1].split("\t"), ["NA", "NA", "NA", "macsNarrow"])
        hm = pd.read_csv(heatmap_tab, sep="\t")
        self.assertTrue(hm["NAsample"].isna().all())
        # but the plot itself is clustered on the scored samples only
        self.assertTrue(os.path.getsize(heatmap) > 0)

    def test_placeholder_plots_when_too_few_scored_samples(self):
        out, snames = self.matrix(nvalid=1)
        pcatab, pcaplot = self.path("pca.tsv"), self.path("pca.pdf")
        heatmap_tab, heatmap = self.path("hm.tsv"), self.path("hm.pdf")
        self.mod.pca_plot(out, snames, "macsNarrow", pcatab, pcaplot)
        self.mod.plot_heatmap(out, heatmap, "macsNarrow", heatmap_tab, snames)
        PCAdata = pd.read_csv(pcatab, sep="\t")
        self.assertEqual(list(PCAdata["sample_name"]), ["S1"])
        self.assertTrue(os.path.getsize(pcaplot) > 0)
        self.assertTrue(os.path.getsize(heatmap) > 0)
        hm = pd.read_csv(heatmap_tab, sep="\t")
        self.assertEqual(list(hm.columns), ["S1", "NAsample", "peakcaller"])


class TestPeakFileHasIntervals(PeakFileTestCase):
    def write_file(self, name, contents):
        path = os.path.join(self.test_dir, name)
        with open(path, "w") as fh:
            fh.write(contents)
        return path

    def test_missing_or_blank_paths(self):
        self.assertFalse(self.mod.peak_file_has_intervals(""))
        self.assertFalse(self.mod.peak_file_has_intervals(None))
        self.assertFalse(
            self.mod.peak_file_has_intervals(os.path.join(self.test_dir, "nope.bed"))
        )

    def test_empty_and_header_only_files(self):
        self.assertFalse(self.mod.peak_file_has_intervals(self.write_file("a.bed", "")))
        self.assertFalse(
            self.mod.peak_file_has_intervals(self.write_file("b.bed", "\n\n"))
        )
        header_only = self.write_file(
            "c.bed", '# comment\ntrack name="peaks"\nbrowser position chr1\n'
        )
        self.assertFalse(self.mod.peak_file_has_intervals(header_only))

    def test_file_with_intervals(self):
        with_peaks = self.write_file("d.bed", "# comment\nchr1\t10\t20\n")
        self.assertTrue(self.mod.peak_file_has_intervals(with_peaks))


if __name__ == "__main__":
    unittest.main()
