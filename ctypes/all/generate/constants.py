#!/usr/bin/env python3
"""
Makes ctypes interface to the C functions in libmath.so.
Generate: 0 inputs, 1 output
For each function, the shape and dtype of the output is specified.
"""

import ctypes
from ctypes import c_size_t
import numpy as np
from time import time

__all__ = ['zeros', 'ones', 'twos',
           'e', 'ln2', 'ln10', 'log2e', 'log10e',
           'sqrt2', 'isqrt2',
           'pi', 'ipi', 'pi_2', 'pi_4',
           'eps', 'realmin', 'realmax',
           'inf', 'nan', 'fill']

FLT_DTYPES = (np.float32, np.complex64)
DTYPES = (np.float32, np.float64, np.complex64, np.complex128)
CLIB = ctypes.cdll.LoadLibrary("libmath.so")
C_FLT_PTR = ctypes.POINTER(ctypes.c_float)
C_DBL_PTR = ctypes.POINTER(ctypes.c_double)


def zeros(shape=(1), dtype=np.float32) -> np.ndarray:
    """
    Generate: Constants: zeros
    """
    assert dtype in DTYPES, f"data type must be in {DTYPES}"

    y = np.empty(shape, dtype=dtype)
    N = y.size

    if y.dtype == np.float32:
        Y = y.ctypes.data_as(C_FLT_PTR)
        ret = CLIB.zeros_s(Y, N)
    elif y.dtype == np.float64:
        Y = y.ctypes.data_as(C_DBL_PTR)
        ret = CLIB.zeros_d(Y, N)
    elif y.dtype == np.complex64:
        Y = y.ctypes.data_as(C_FLT_PTR)
        ret = CLIB.zeros_c(Y, N)
    else:
        Y = y.ctypes.data_as(C_DBL_PTR)
        ret = CLIB.zeros_z(Y, N)

    assert ret == 0, "error during call to C function"

    return y


def ones(shape=(1), dtype=np.float32) -> np.ndarray:
    """
    Generate: Constants: ones
    """
    assert dtype in DTYPES, f"data type must be in {DTYPES}"

    y = np.empty(shape, dtype=dtype)
    N = y.size

    if y.dtype == np.float32:
        Y = y.ctypes.data_as(C_FLT_PTR)
        ret = CLIB.ones_s(Y, N)
    elif y.dtype == np.float64:
        Y = y.ctypes.data_as(C_DBL_PTR)
        ret = CLIB.ones_d(Y, N)
    elif y.dtype == np.complex64:
        Y = y.ctypes.data_as(C_FLT_PTR)
        ret = CLIB.ones_c(Y, N)
    else:
        Y = y.ctypes.data_as(C_DBL_PTR)
        ret = CLIB.ones_z(Y, N)

    assert ret == 0, "error during call to C function"

    return y


def twos(shape=(1), dtype=np.float32) -> np.ndarray:
    """
    Generate: Constants: twos
    """
    assert dtype in DTYPES, f"data type must be in {DTYPES}"

    y = np.empty(shape, dtype=dtype)
    N = y.size

    if y.dtype == np.float32:
        Y = y.ctypes.data_as(C_FLT_PTR)
        ret = CLIB.twos_s(Y, N)
    elif y.dtype == np.float64:
        Y = y.ctypes.data_as(C_DBL_PTR)
        ret = CLIB.twos_d(Y, N)
    elif y.dtype == np.complex64:
        Y = y.ctypes.data_as(C_FLT_PTR)
        ret = CLIB.twos_c(Y, N)
    else:
        Y = y.ctypes.data_as(C_DBL_PTR)
        ret = CLIB.twos_z(Y, N)

    assert ret == 0, "error during call to C function"

    return y


def e(shape=(1), dtype=np.float32) -> np.ndarray:
    """
    Generate: Constants: e
    """
    assert dtype in DTYPES, f"data type must be in {DTYPES}"

    y = np.empty(shape, dtype=dtype)
    N = y.size

    if y.dtype == np.float32:
        Y = y.ctypes.data_as(C_FLT_PTR)
        ret = CLIB.e_s(Y, N)
    elif y.dtype == np.float64:
        Y = y.ctypes.data_as(C_DBL_PTR)
        ret = CLIB.e_d(Y, N)
    elif y.dtype == np.complex64:
        Y = y.ctypes.data_as(C_FLT_PTR)
        ret = CLIB.e_c(Y, N)
    else:
        Y = y.ctypes.data_as(C_DBL_PTR)
        ret = CLIB.e_z(Y, N)

    assert ret == 0, "error during call to C function"

    return y


def ln2(shape=(1), dtype=np.float32) -> np.ndarray:
    """
    Generate: Constants: ln2
    """
    assert dtype in DTYPES, f"data type must be in {DTYPES}"

    y = np.empty(shape, dtype=dtype)
    N = y.size

    C_PTR = C_FLT_PTR if dtype in FLT_DTYPES else C_DBL_PTR
    Y = y.ctypes.data_as(C_PTR)

    if y.dtype == np.float32:
        # Y = y.ctypes.data_as(C_FLT_PTR)
        ret = CLIB.ln2_s(Y, N)
    elif y.dtype == np.float64:
        # Y = y.ctypes.data_as(C_DBL_PTR)
        ret = CLIB.ln2_d(Y, N)
    elif y.dtype == np.complex64:
        # Y = y.ctypes.data_as(C_FLT_PTR)
        ret = CLIB.ln2_c(Y, N)
    else:
        # Y = y.ctypes.data_as(C_DBL_PTR)
        ret = CLIB.ln2_z(Y, N)

    assert ret == 0, "error during call to C function"

    return y


def ln10(shape=(1), dtype=np.float32) -> np.ndarray:
    pass


def log2e(shape=(1), dtype=np.float32) -> np.ndarray:
    pass


def log10e(shape=(1), dtype=np.float32) -> np.ndarray:
    pass


def sqrt2(shape=(1), dtype=np.float32) -> np.ndarray:
    pass


def isqrt2(shape=(1), dtype=np.float32) -> np.ndarray:
    pass


def pi(shape=(1), dtype=np.float32) -> np.ndarray:
    pass


def ipi(shape=(1), dtype=np.float32) -> np.ndarray:
    pass


def pi_2(shape=(1), dtype=np.float32) -> np.ndarray:
    pass


def pi_4(shape=(1), dtype=np.float32) -> np.ndarray:
    pass


def eps(shape=(1), dtype=np.float32) -> np.ndarray:
    pass


def realmin(shape=(1), dtype=np.float32) -> np.ndarray:
    pass


def realmax(shape=(1), dtype=np.float32) -> np.ndarray:
    pass


def inf(shape=(1), dtype=np.float32) -> np.ndarray:
    pass


def nan(shape=(1), dtype=np.float32) -> np.ndarray:
    pass


def fill(x=0, shape=(1), dtype=np.float32) -> np.ndarray:
    pass


# Main
def main() -> np.ndarray:
    """
    Only used for quick command-line test
    """
    shape = (2, 3)
    # shape = (8000, 300)

    tic = time()
    y = np.ones(shape)
    print(time()-tic)
    print(y)

    tic = time()
    y = ln2(shape)
    print(time()-tic)
    print(y)

    return y


if __name__ == "__main__":
    main()
