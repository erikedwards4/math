#!/usr/bin/env python3
"""
Makes ctypes interface to the C functions in libmath.so.
Elementwise1: 1 input, 1 output with same shape as input
Each function works element-wise (1 element at a time).
Operators: meant to emulate x++, x-- and -x.
"""

from time import time
import ctypes
import numpy as np

__all__ = ['plusplus', 'minusminus', 'neg']

# Global constants
REAL_DTYPES = (np.float32, np.float64)
CPLX_DTYPES = (np.complex64, np.complex128)
DTYPES = REAL_DTYPES + CPLX_DTYPES
CLIB = ctypes.cdll.LoadLibrary("libmath.so")
C_FLT_PTR = ctypes.POINTER(ctypes.c_float)
C_DBL_PTR = ctypes.POINTER(ctypes.c_double)


def plusplus(x) -> np.ndarray:
    """
    Elementwise1: Operators: plusplus
    This always operates inplace (like x++).
    """
    assert x.dtype in REAL_DTYPES, f"input data type must be in {REAL_DTYPES}"
    N = x.size

    if x.dtype == np.float32:
        X = x.ctypes.data_as(C_FLT_PTR)
        ret = CLIB.plusplus_s(X, N)
    else:
        X = x.ctypes.data_as(C_DBL_PTR)
        ret = CLIB.plusplus_d(X, N)

    assert ret == 0, "error during call to C function"

    return x


def minusminus(x) -> np.ndarray:
    """
    Elementwise1: Operators: minusminus
    This always operates inplace (like x--).
    """
    assert x.dtype in REAL_DTYPES, f"input data type must be in {REAL_DTYPES}"
    N = x.size

    if x.dtype == np.float32:
        X = x.ctypes.data_as(C_FLT_PTR)
        ret = CLIB.minusminus_s(X, N)
    else:
        X = x.ctypes.data_as(C_DBL_PTR)
        ret = CLIB.minusminus_d(X, N)

    assert ret == 0, "error during call to C function"

    return x


def neg(x, inplace=False) -> np.ndarray:
    """
    Elementwise1: Operators: neg
    This is like -x.
    """
    assert x.dtype in DTYPES, f"input data type must be in {DTYPES}"
    N = x.size

    if inplace:
        if x.dtype == np.float32:
            X = x.ctypes.data_as(C_FLT_PTR)
            ret = CLIB.neg_inplace_s(X, N)
        elif x.dtype == np.float64:
            X = x.ctypes.data_as(C_DBL_PTR)
            ret = CLIB.neg_inplace_d(X, N)
        elif x.dtype == np.complex64:
            X = x.ctypes.data_as(C_FLT_PTR)
            ret = CLIB.neg_inplace_c(X, N)
        else:
            X = x.ctypes.data_as(C_DBL_PTR)
            ret = CLIB.neg_inplace_z(X, N)
        assert ret == 0, "error during call to C function"
        return x

    else:
        y = np.empty_like(x)
        if x.dtype == np.float32:
            X = x.ctypes.data_as(C_FLT_PTR)
            Y = y.ctypes.data_as(C_FLT_PTR)
            ret = CLIB.neg_s(Y, X, N)
        elif x.dtype == np.float64:
            X = x.ctypes.data_as(C_DBL_PTR)
            Y = y.ctypes.data_as(C_DBL_PTR)
            ret = CLIB.neg_d(Y, X, N)
        elif x.dtype == np.complex64:
            X = x.ctypes.data_as(C_FLT_PTR)
            Y = y.ctypes.data_as(C_FLT_PTR)
            ret = CLIB.neg_c(Y, X, N)
        else:
            X = x.ctypes.data_as(C_DBL_PTR)
            Y = y.ctypes.data_as(C_DBL_PTR)
            ret = CLIB.neg_z(Y, X, N)
        assert ret == 0, "error during call to C function"
        return y


# Main
def main() -> None:
    """
    Only used for quick command-line test
    """
    x = np.random.randn(2, 3)
    x = x + 1j*x
    print(x)

    tic = time()
    x += 1
    print(time()-tic)
    print(x)

    tic = time()
    plusplus(x)
    print(time()-tic)
    print(x)


if __name__ == "__main__":
    main()
