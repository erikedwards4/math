#!/usr/bin/env python3
"""
Makes ctypes interface to the C functions in libmath.so.
Elementwise1: 1 input, 1 output with same shape as input
Each function works element-wise (1 element at a time).
Nonlin: various static nonlinearities.
"""

from time import time
import ctypes
import numpy as np

__all__ = ['abs', 'square']

REAL_DTYPES = (np.float32, np.float64)
CPLX_DTYPES = (np.complex64, np.complex128)
DTYPES = REAL_DTYPES + CPLX_DTYPES
CLIB = ctypes.cdll.LoadLibrary("libmath.so")
C_FLT_PTR = ctypes.POINTER(ctypes.c_float)
C_DBL_PTR = ctypes.POINTER(ctypes.c_double)


def abs(x, inplace=False) -> np.ndarray:
    """
    Elementwise1: Nonlin: abs
    """
    N = x.size

    if inplace:
        assert x.dtype in REAL_DTYPES, f"input data type must be in {REAL_DTYPES}"
        if x.dtype == np.float32:
            X = x.ctypes.data_as(C_FLT_PTR)
            ret = CLIB.abs_inplace_s(X, N)
        elif x.dtype == np.float64:
            X = x.ctypes.data_as(C_DBL_PTR)
            ret = CLIB.abs_inplace_d(X, N)
        assert ret == 0, "error during call to C function"
        return x

    else:
        assert x.dtype in DTYPES, f"input data type must be in {DTYPES}"
        if x.dtype == np.float32:
            y = np.empty_like(x)
            X = x.ctypes.data_as(C_FLT_PTR)
            Y = y.ctypes.data_as(C_FLT_PTR)
            ret = CLIB.abs_s(Y, X, N)
        elif x.dtype == np.float64:
            y = np.empty_like(x)
            X = x.ctypes.data_as(C_DBL_PTR)
            Y = y.ctypes.data_as(C_DBL_PTR)
            ret = CLIB.abs_d(Y, X, N)
        elif x.dtype == np.complex64:
            y = np.empty_like(x, dtype=np.float32)
            X = x.ctypes.data_as(C_FLT_PTR)
            Y = y.ctypes.data_as(C_FLT_PTR)
            ret = CLIB.abs_c(Y, X, N)
        else:
            y = np.empty_like(x, dtype=np.float64)
            X = x.ctypes.data_as(C_DBL_PTR)
            Y = y.ctypes.data_as(C_DBL_PTR)
            ret = CLIB.abs_z(Y, X, N)
        assert ret == 0, "error during call to C function"
        return y


def square(x, inplace=False) -> np.ndarray:
    """
    Elementwise1: Nonlin: square
    """
    N = x.size

    if inplace:
        assert x.dtype in REAL_DTYPES, f"input data type must be in {REAL_DTYPES}"
        if x.dtype == np.float32:
            X = x.ctypes.data_as(C_FLT_PTR)
            ret = CLIB.square_inplace_s(X, N)
        elif x.dtype == np.float64:
            X = x.ctypes.data_as(C_DBL_PTR)
            ret = CLIB.square_inplace_d(X, N)
        assert ret == 0, "error during call to C function"
        return x

    else:
        assert x.dtype in DTYPES, f"input data type must be in {DTYPES}"
        if x.dtype == np.float32:
            y = np.empty_like(x)
            X = x.ctypes.data_as(C_FLT_PTR)
            Y = y.ctypes.data_as(C_FLT_PTR)
            ret = CLIB.square_s(Y, X, N)
        elif x.dtype == np.float64:
            y = np.empty_like(x)
            X = x.ctypes.data_as(C_DBL_PTR)
            Y = y.ctypes.data_as(C_DBL_PTR)
            ret = CLIB.square_d(Y, X, N)
        elif x.dtype == np.complex64:
            y = np.empty_like(x, dtype=np.float32)
            X = x.ctypes.data_as(C_FLT_PTR)
            Y = y.ctypes.data_as(C_FLT_PTR)
            ret = CLIB.square_c(Y, X, N)
        else:
            y = np.empty_like(x, dtype=np.float64)
            X = x.ctypes.data_as(C_DBL_PTR)
            Y = y.ctypes.data_as(C_DBL_PTR)
            ret = CLIB.square_z(Y, X, N)
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
    y = np.square(x)
    print(time()-tic)
    print(y)

    tic = time()
    y = square(x, inplace=False)
    print(time()-tic)
    print(y)


if __name__ == "__main__":
    main()
