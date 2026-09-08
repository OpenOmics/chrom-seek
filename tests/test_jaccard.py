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
        "union": union,
        "jaccard": jaccard,
        "n_intersections": n_intersections,
        "fileA": os.path.basename(fileA),
        "fileB": os.path.basename(fileB),
    }
    keylist = sorted(record.keys())
    return ([str(record[key]) for key in keylist], keylist)


class TestNanToZero(unittest.TestCase):
    def setUp(self):
        self.mod = load_jaccard_score_module()

    def test_undefined_scores_are_zero(self):
        for score in ("nan", "NaN", float("nan"), "", ".", None):
            self.assertEqual(self.mod.nan_to_zero(score), 0.0)

    def test_numpy_and_pandas_nans_are_zero(self):
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
            self.assertEqual(self.mod.nan_to_zero(score), 0.0)

    def test_real_scores_are_preserved(self):
        self.assertEqual(self.mod.nan_to_zero("0.3333"), 0.3333)
        self.assertEqual(self.mod.nan_to_zero(0.5), 0.5)
        self.assertEqual(self.mod.nan_to_zero("0"), 0.0)
        self.assertEqual(self.mod.nan_to_zero("1"), 1.0)
        # numpy/pandas scalars that are not missing values
        self.assertEqual(self.mod.nan_to_zero(np.float64(0.25)), 0.25)
        self.assertEqual(self.mod.nan_to_zero(np.float32(0.5)), 0.5)
        self.assertEqual(self.mod.nan_to_zero(pd.Series([0.75])[0]), 0.75)

    def test_non_scalar_scores_are_zero(self):
        """pd.isna on a non-scalar is ambiguous, it must not raise"""
        self.assertEqual(self.mod.nan_to_zero(np.array([np.nan, 1.0])), 0.0)
        self.assertEqual(self.mod.nan_to_zero([np.nan]), 0.0)


class TestRunJaccard(unittest.TestCase):
    def setUp(self):
        self.mod = load_jaccard_score_module()

    def test_jaccard_is_the_fourth_column(self):
        """
        loop_jaccard pulls the score out of data[3], which only holds the
        jaccard value because the bedtools keys are sorted alphabetically.
        """
        bedtool = mock.MagicMock()
        bedtool.sort.return_value = bedtool
        bedtool.jaccard.return_value = {
            "intersection": 50,
            "union": 150,
            "jaccard": 0.333333,
            "n_intersections": 1,
        }
        with mock.patch.object(self.mod, "BedTool", return_value=bedtool):
            data, keylist = self.mod.run_jaccard("/tmp/a.bed", "/tmp/b.bed", "genome.txt")
        self.assertEqual(keylist[3], "jaccard")
        self.assertEqual(float(data[3]), 0.333333)
        self.assertEqual(data[keylist.index("fileA")], "a.bed")
        self.assertEqual(data[keylist.index("fileB")], "b.bed")


class TestLoopJaccard(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.mod = load_jaccard_score_module()
        self.patcher = mock.patch.object(self.mod, "run_jaccard", fake_run_jaccard)
        self.patcher.start()

    def tearDown(self):
        self.patcher.stop()
        shutil.rmtree(self.test_dir)

    def write_peaks(self, name, intervals):
        path = os.path.join(self.test_dir, name)
        with open(path, "w") as peaks:
            for chrom, start, end in intervals:
                peaks.write("\t".join([chrom, str(start), str(end)]) + "\n")
        return path

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
        # one header line plus one line per pairwise comparison
        self.assertEqual(len(outTable), 1 + 3)
        self.assertEqual(outTable[0].split("\t")[3], "jaccard")

    def test_multi_interval_jaccard_scores(self):
        fileA = self.write_peaks(
            "A_peaks.bed", [("chr1", 0, 100), ("chr1", 200, 300), ("chr2", 0, 100)]
        )
        fileB = self.write_peaks("B_peaks.bed", [("chr1", 75, 250), ("chr3", 0, 100)])
        outTable, out, snames = self.mod.loop_jaccard([fileA, fileB], "genome.txt")
        # intersection: 25bp (chr1 75-100) + 50bp (chr1 200-250) = 75bp
        # union: 300bp + 275bp - 75bp = 500bp
        self.assertAlmostEqual(out.loc["A", "B"], 75 / 500)

    def test_empty_peak_file_scores_as_zero(self):
        """An empty peak file gives an undefined (NaN) jaccard, recorded as 0"""
        fileA = self.write_peaks("A_peaks.bed", [("chr1", 0, 100)])
        empty = self.write_peaks("B_peaks.bed", [])
        both_empty = self.write_peaks("C_peaks.bed", [])
        # bedtools itself reports NaN for these comparisons
        self.assertTrue(math.isnan(reference_jaccard(empty, both_empty)[2]))
        outTable, out, snames = self.mod.loop_jaccard(
            [fileA, empty, both_empty], "genome.txt"
        )
        self.assertFalse(out.isna().any().any())
        self.assertEqual(out.loc["A", "B"], 0.0)
        self.assertEqual(out.loc["B", "C"], 0.0)
        self.assertEqual(out.loc["B", "B"], 1.0)

    def test_numpy_nan_score_never_reaches_the_matrix(self):
        """A numpy NaN handed back by run_jaccard is still scored as 0"""
        fileA = self.write_peaks("A_peaks.bed", [("chr1", 0, 100)])
        fileB = self.write_peaks("B_peaks.bed", [("chr1", 50, 150)])

        def numpy_nan_run_jaccard(fileA, fileB, genomefile):
            # str(numpy.nan), which is what a numpy NaN looks like once
            # run_jaccard has stringified the bedtools record
            data, keylist = fake_run_jaccard(fileA, fileB, genomefile)
            data[keylist.index("jaccard")] = np.str_(np.nan)
            return (data, keylist)

        with mock.patch.object(self.mod, "run_jaccard", numpy_nan_run_jaccard):
            outTable, out, snames = self.mod.loop_jaccard([fileA, fileB], "genome.txt")
        self.assertFalse(out.isna().any().any())
        self.assertEqual(out.loc["A", "B"], 0.0)
        # the raw bedtools record is still reported verbatim
        self.assertEqual(outTable[1].split("\t")[3], "nan")

    def test_single_file_matrix(self):
        fileA = self.write_peaks("A_peaks.bed", [("chr1", 0, 100)])
        outTable, out, snames = self.mod.loop_jaccard([fileA], "genome.txt")
        self.assertEqual(out.shape, (1, 1))
        self.assertEqual(out.loc["A", "A"], 1.0)
        self.assertEqual(outTable, [])


class TestPeakFileHasIntervals(unittest.TestCase):
    """
    main() drops empty/missing peak files before any comparison is run, so
    these files never reach the NaN handling in loop_jaccard.
    """

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.mod = load_jaccard_score_module()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

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
