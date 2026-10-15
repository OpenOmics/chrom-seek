import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "merge_uropa_diffbind.py"
DIFFBIND_HEADER = "seqnames\tstart\tend\tstrand\tConc\tFold\tFDR\n"
UROPA_TABLE = (
    "peak_chr\tpeak_start\tpeak_end\tpeak_id\tpeak_score\tpeak_strand\tgene_id\n"
    "chr1\t100\t150\tPeak1\t0\t.\tgene1\n"
    "chr1\t200\t250\tPeak2\t0\t.\tgene2\n"
)


class TestMergeUropaDiffbind(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.diffbind = Path(self.temp_dir.name) / "diffbind.tab"
        self.uropa = Path(self.temp_dir.name) / "uropa.txt"
        self.output = Path(self.temp_dir.name) / "merged.txt"
        self.uropa.write_text(UROPA_TABLE)

    def run_script(self, rows, *options):
        self.diffbind.write_text(DIFFBIND_HEADER + rows)
        return subprocess.run(
            [
                sys.executable, str(SCRIPT),
                "--diffbind", str(self.diffbind),
                "--uropa", str(self.uropa),
                "--output", str(self.output),
                *options,
            ],
            capture_output=True,
            text=True,
        )

    def test_no_differential_peaks_reports_analysis_result(self):
        cases = [
            ("header_only", "", [], "0.05", "0"),
            ("fdr_filtered", "chr1\t101\t150\t*\t5\t2\t0.1\n", [], "0.05", "0"),
            (
                "fold_filtered", "chr1\t101\t150\t*\t5\t-0.5\t0.01\n",
                ["--fdr", "0.02", "--fold", "1"], "0.02", "1.0",
            ),
        ]
        for name, rows, options, fdr, fold in cases:
            with self.subTest(name=name):
                result = self.run_script(rows, *options)
                self.assertEqual(result.returncode, 1)
                self.assertIn("No differential peaks identified", result.stderr)
                self.assertIn(str(self.diffbind), result.stderr)
                self.assertIn(f"FDR <= {fdr}", result.stderr)
                self.assertIn(f"absolute fold-change >= {fold}", result.stderr)
                self.assertIn("not a software error", result.stderr)
                self.assertNotIn("Traceback", result.stderr)
                self.assertIn("0 total peaks remaining", result.stdout)
                self.assertFalse(self.output.exists())

    def test_passing_peaks_are_annotated_and_sorted(self):
        result = self.run_script(
            "chr1\t101\t150\t*\t5\t-1\t0.05\n"
            "chr1\t201\t250\t*\t6\t2\t0.01\n"
            "chr1\t301\t350\t*\t5\t0.5\t0.01\n"
            "chr1\t401\t450\t*\t5\t2\t0.1\n",
            "--fold", "1",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("No differential peaks identified", result.stderr)
        self.assertIn("2 total peaks remaining", result.stdout)
        with self.output.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["peak_id"] for row in rows], ["Peak2", "Peak1"])
        self.assertEqual([row["gene_id"] for row in rows], ["gene2", "gene1"])
        self.assertEqual([float(row["FDR"]) for row in rows], [0.01, 0.05])
        self.assertEqual([float(row["Fold"]) for row in rows], [2, -1])
        self.assertNotIn("strand", rows[0])
        self.assertNotIn("Conc", rows[0])


if __name__ == "__main__":
    unittest.main()
