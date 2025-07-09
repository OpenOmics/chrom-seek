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
        self.assertIsNone(blocks['two'])
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
        self.assertEqual(input_vals, [None])
        self.assertEqual(blocks['three'], 'two')
        self.assertIsNone(blocks['two'])


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
        self.assertEqual(input_vals, [None])
        # only Nones for blocks values
        self.assertEqual(set(blocks.values()), {None})


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



if __name__ == '__main__':
    unittest.main()