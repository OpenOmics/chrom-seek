#!/usr/bin/env python
import unittest
import shutil
import tempfile
import os
from src.files import peakcalls


class TestParsePeakCall(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_full(self):
        str2file = \
            "Sample\tInputControl\tGroup\tBlock\n" + \
            "one\t\tg1\t\n" + \
            "two\tone\tg2\t\n" + \
            "three\tone\tg3,g4\ttwo"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        chip2input, groups, blocks = peakcalls(path)
        # chips and input mapping
        self.assertEqual([v for v in chip2input.values()], ['', 'one', 'one'])    
        # grouping
        for g in ('g1', 'g2', 'g3', 'g4'):
            self.assertIn(g, groups.keys())
        self.assertEqual(groups['g4'], ['three'])
        # blocking
        for b in ('one', 'two', 'three'):
            self.assertIn(b, blocks.keys())
        self.assertEqual(blocks['two'], '')
        self.assertEqual(blocks['three'], 'two')
    

    def test_optional_one(self):
        str2file = \
            "Sample\tGroup\tBlock\n" + \
            "one\tg1\t\n" + \
            "two\tg2\t\n" + \
            "three\tg3,g4\ttwo"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        chip2input, groups, blocks = peakcalls(path)
        input_vals = list(set([v for v in chip2input.values()]))
        # only Nones for input values
        self.assertEqual(input_vals, [''])
        self.assertEqual(blocks['three'], 'two')
        self.assertEqual(blocks['two'], '')


    def test_optional_two(self):
        str2file = \
            "Sample\tGroup\n" + \
            "one\tg1\n" + \
            "two\tg2\n" + \
            "three\tg3,g4"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        chip2input, groups, blocks = peakcalls(path)
        input_vals = list(set([v for v in chip2input.values()]))
        # only Nones for input values
        self.assertEqual(input_vals, [''])
        # only Nones for blocks values
        self.assertEqual(set(blocks.values()), {''})


    def test_bad_group_labels(self):
        str2file = \
            "Sample\tGroup\n" + \
            "one\tg_1\n" + \
            "two\tg-2\n" + \
            "three\tg*3,g4" + \
            "four\tg4"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        with self.assertRaisesRegex(ValueError, r"Group\(s\)\: g_1, g-2, g\*3"):
            chip2input, groups, blocks = peakcalls(path)


    def test_bad_headers(self):
        str2file = \
            "Sample\tGroup\tETC\n" + \
            "one\tg_1\tone\n" + \
            "two\tg-2\tone\n" + \
            "three\tg*3,g4\tone\n" + \
            "four\tg4\tone\n"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        with self.assertRaisesRegex(ValueError, r"Peakcall file has unsupported headers!"):
            chip2input, groups, blocks = peakcalls(path)


    def test_with_space(self):
        str2file = \
            "Sample\tGroup\tETC\n" + \
            "one\tg_1  one\n" + \
            "two\tg-2\tone\n" + \
            "three\tg*3,g4\tone\n" + \
            "four\tg4\tone\n"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        with self.assertRaisesRegex(ValueError, r"Spaces exist in peakcall file"):
            chip2input, groups, blocks = peakcalls(path)


    def test_reps_no_input_col(self):
        str2file = "Sample\tGroup\n" + \
                    "WT_S1\tG1\n" + \
                    "WT_S2\tG1\n" + \
                    "WT_S3\tG2\n" + \
                    "WT_S4\tG2\n"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        chip2input, groups, blocks = peakcalls(path)
        self.assertEqual(set(chip2input.values()), {''})


    def test_atac_with_input_control(self):
        """Test that ATAC assay raises error when InputControl column is not blank"""
        str2file = \
            "Sample\tInputControl\tGroup\tBlock\n" + \
            "ATAC_S1\tInput_S1\tG1\tB1\n" + \
            "ATAC_S2\tInput_S2\tG2\tB2\n"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        with self.assertRaisesRegex(ValueError, r"ATAC-seq assay does not support InputControls, please leave InputColumn empty"):
            chip2input, groups, blocks = peakcalls(path, assay="atac")


    def test_atac_without_input_control(self):
        """Test that ATAC assay works correctly when InputControl column is blank"""
        str2file = \
            "Sample\tInputControl\tGroup\tBlock\n" + \
            "ATAC_S1\t\tG1\tB1\n" + \
            "ATAC_S2\t\tG2\tB2\n"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        # This should not raise an error
        chip2input, groups, blocks = peakcalls(path, assay="atac")
        self.assertEqual(chip2input['ATAC_S1'], '')
        self.assertEqual(chip2input['ATAC_S2'], '')


    def test_input_control_equals_sample_positive(self):
        """Test that error is raised when InputControl equals Sample"""
        str2file = \
            "Sample\tInputControl\tGroup\tBlock\n" + \
            "Sample1\tSample1\tG1\tB1\n" + \
            "Sample2\tInput2\tG2\tB2\n"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        with self.assertRaisesRegex(ValueError, r"InputControl cannot be the same as Sample. All these samples have an error: Sample1"):
            chip2input, groups, blocks = peakcalls(path)


    def test_input_control_equals_sample_negative(self):
        """Test that no error is raised when InputControl does not equal Sample"""
        str2file = \
            "Sample\tInputControl\tGroup\tBlock\n" + \
            "Sample1\tInput1\tG1\tB1\n" + \
            "Sample2\tInput2\tG2\tB2\n" + \
            "Sample3\t\tG3\tB3\n"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        # This should not raise an error
        chip2input, groups, blocks = peakcalls(path)
        self.assertEqual(chip2input['Sample1'], 'Input1')
        self.assertEqual(chip2input['Sample2'], 'Input2')
        self.assertEqual(chip2input['Sample3'], '')


    def test_block_equals_sample_positive(self):
        """Test that error is raised when Block equals Sample"""
        str2file = \
            "Sample\tInputControl\tGroup\tBlock\n" + \
            "Sample1\tInput1\tG1\tSample1\n" + \
            "Sample2\tInput2\tG2\tB1\n"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        with self.assertRaisesRegex(ValueError, r"Block cannot be the same as Sample. All these samples have an error: Sample1"):
            chip2input, groups, blocks = peakcalls(path)


    def test_block_equals_sample_negative(self):
        """Test that no error is raised when Block does not equal Sample"""
        str2file = \
            "Sample\tInputControl\tGroup\tBlock\n" + \
            "Sample1\tInput1\tG1\tB1\n" + \
            "Sample2\tInput2\tG2\tB2\n" + \
            "Sample3\tInput3\tG3\t\n"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        # This should not raise an error
        chip2input, groups, blocks = peakcalls(path)
        self.assertEqual(blocks['Sample1'], 'B1')
        self.assertEqual(blocks['Sample2'], 'B2')
        self.assertEqual(blocks['Sample3'], '')


    def test_group_equals_sample_positive(self):
        """Test that error is raised when Group equals Sample"""
        str2file = \
            "Sample\tInputControl\tGroup\tBlock\n" + \
            "Sample1\tInput1\tSample1\tB1\n" + \
            "Sample2\tInput2\tG2\tB2\n"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        with self.assertRaisesRegex(ValueError, r"Group cannot be the same as Sample. All these samples have an error: Sample1"):
            chip2input, groups, blocks = peakcalls(path)


    def test_group_equals_sample_negative(self):
        """Test that no error is raised when Group does not equal Sample"""
        str2file = \
            "Sample\tInputControl\tGroup\tBlock\n" + \
            "Sample1\tInput1\tG1\tB1\n" + \
            "Sample2\tInput2\tG2\tB2\n" + \
            "Sample3\tInput3\tG3,G4\tB3\n"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        # This should not raise an error
        chip2input, groups, blocks = peakcalls(path)
        self.assertIn('G1', groups)
        self.assertIn('G2', groups)
        self.assertIn('G3', groups)
        self.assertIn('G4', groups)
        self.assertEqual(groups['G1'], ['Sample1'])
        self.assertEqual(groups['G2'], ['Sample2'])
        self.assertEqual(groups['G3'], ['Sample3'])


    def test_multiple_validation_errors(self):
        """Test that multiple samples can violate the same rule"""
        str2file = \
            "Sample\tInputControl\tGroup\tBlock\n" + \
            "Sample1\tSample1\tG1\tB1\n" + \
            "Sample2\tSample2\tG2\tB2\n"
        path = os.path.join(self.test_dir, 'test.txt')
        with open(path, 'w') as f:
            f.write(str2file)
            f.close()
        # Should raise error mentioning both samples
        with self.assertRaisesRegex(ValueError, r"InputControl cannot be the same as Sample"):
            chip2input, groups, blocks = peakcalls(path)


if __name__ == '__main__':
    unittest.main()