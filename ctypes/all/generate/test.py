#!/usr/bin/env python3
"""
Makes ctypes interface to the C functions in libmath.so.
Generate: 0 inputs, 1 output
For each function, the shape and dtype of the output is specified.
"""

import unittest
import ctypes
from ctypes import c_size_t
import numpy as np
from time import time
from constants import *


class TestFunctions(unittest.TestCase):

    def test_zeros(self):
        x = zeros((2, 3), dtype=np.float32)
        y = np.zeros((2, 3), dtype=np.float32)
        self.assertTrue((x == y).all())
    
    def test_ones(self):
        x = ones((2, 3), dtype=np.complex64)
        y = np.ones((2, 3), dtype=np.complex64)
        self.assertTrue((x == y).all())


if __name__ == '__main__':
    unittest.main()
