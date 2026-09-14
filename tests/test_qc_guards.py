#!/usr/bin/env python
import csv
import os
import shutil
import tempfile
import unittest
import importlib.util
from argparse import Namespace


def load_prep_diffbindqc_module():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    module_path = os.path.join(repo_root, "bin", "prep_diffbindQC.py")
    spec = importlib.util.spec_from_file_location("prep_diffbindQC", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestPrepDiffbindQCGuards(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.mod = load_prep_diffbindqc_module()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_nan_like(self):
        self.assertTrue(self.mod.nan_like(None))
        self.assertTrue(self.mod.nan_like(""))
        self.assertTrue(self.mod.nan_like("NaN"))
        self.assertFalse(self.mod.nan_like("S1"))

    def test_main_drops_empty_and_nan_peak_rows(self):
        # one valid sample, one empty peak file, one missing/NaN-like peak path
        cfg = {
            "project": {
                "peaks": {
                    "inputs": {"S1": "", "S2": "", "S3": ""},
                    "chips": ["S1", "S2", "S3"],
                },
                "groups": {"G1": ["S1", "S2", "S3"]},
            }
        }

        bam_s1 = os.path.join(self.test_dir, "S1.Q5DD.bam")
        bam_s2 = os.path.join(self.test_dir, "S2.Q5DD.bam")
        bam_s3 = os.path.join(self.test_dir, "S3.Q5DD.bam")
        open(bam_s1, "a").close()
        open(bam_s2, "a").close()
        open(bam_s3, "a").close()

        peak_s1 = os.path.join(self.test_dir, "S1.narrowPeak")
        peak_s2 = os.path.join(self.test_dir, "S2.narrowPeak")
        with open(peak_s1, "w") as fh:
            fh.write("chr1\t10\t20\n")
        with open(peak_s2, "w") as fh:
            fh.write("")

        out_csv = os.path.join(self.test_dir, "AllSamples-DiffBind_prep.csv")
        args = Namespace(
            pktool="narrow",
            output=out_csv,
            peaks=[peak_s1, peak_s2],
            bams=[bam_s1, bam_s2, bam_s3],
            cfg=cfg,
        )

        self.mod.main(args)

        with open(out_csv, "r") as fh:
            rows = list(csv.DictReader(fh))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["SampleID"], "S1")
        self.assertEqual(rows[0]["Peaks"], peak_s1)


if __name__ == "__main__":
    unittest.main()
