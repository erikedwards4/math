#!/usr/bin/env python3
"""
Makes ctypes interface to the C functions in libmath.so.
Vec2scalar: 1 input, 1 output with vectors reduced to scalars along dim.
Each function works vector-wise (1 vector at a time).
"""

from time import time
import ctypes
from ctypes import c_size_t
import numpy as np

__all__ = []

FLT_DTYPES = (np.float32, np.complex64)
REAL_DTYPES = (np.float32, np.float64)
CPLX_DTYPES = (np.complex64, np.complex128)
DTYPES = REAL_DTYPES + CPLX_DTYPES
CLIB = ctypes.cdll.LoadLibrary("libmath.so")
C_FLT_PTR = ctypes.POINTER(ctypes.c_float)
C_DBL_PTR = ctypes.POINTER(ctypes.c_double)



# Main
def main() -> None:
    """
    Only used for quick command-line test
    """
    # x1 = np.random.randn(2, 3)
    # x2 = np.random.randn(2, 3)
    # x1 = x1 + 1j*x1
    # x2 = x2 + 1j*x2


    tic = time()
    # y = np.
    print(time()-tic)
    # print(y)

    tic = time()
    # y = new(x1, x2)
    print(time()-tic)
    # print(y)


if __name__ == "__main__":
    main()
