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

__all__ = ['range', 'iqr', 'idr']

FLT_DTYPES = (np.float32, np.complex64)
REAL_DTYPES = (np.float32, np.float64)
CPLX_DTYPES = (np.complex64, np.complex128)
DTYPES = REAL_DTYPES + CPLX_DTYPES
CLIB = ctypes.cdll.LoadLibrary("libmath.so")
C_FLT_PTR = ctypes.POINTER(ctypes.c_float)
C_DBL_PTR = ctypes.POINTER(ctypes.c_double)


def range(x, axis=0) -> np.ndarray:
    """
    Vec2scalar: Ranges: range
    """
    # Check
    DTYPE = x.dtype
    COLMAJOR = not x.flags['C_CONTIGUOUS']
    assert DTYPE in REAL_DTYPES, f"input data type must be in {REAL_DTYPES}"
    assert x.ndim < 5, "input x must have ndim < 5"
    dim = x.ndim + axis if axis < 0 else axis
    assert 0 <= dim < x.ndim, "axis out of range [0 x.ndim)"

    # Get input/output shapes
    R = x.shape[0] if x.ndim > 0 else 1
    C = x.shape[1] if x.ndim > 1 else 1
    S = x.shape[2] if x.ndim > 2 else 1
    H = x.shape[3] if x.ndim > 3 else 1
    Ry = 1 if dim == 0 else R
    Cy = 1 if dim == 1 else C
    Sy = 1 if dim == 2 else S
    Hy = 1 if dim == 3 else H
    R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)
    DIM = c_size_t(dim)

    # Make pointer
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X = x.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((Ry, Cy, Sy, Hy), dtype=DTYPE)
    while y.ndim > x.ndim:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.range_s(Y, X, R, C, S, H, COLMAJOR, DIM)
    else:
        ret = CLIB.range_d(Y, X, R, C, S, H, COLMAJOR, DIM)
    assert ret == 0, "error during call to C function"

    return y


def iqr(x, axis=0) -> np.ndarray:
    """
    Vec2scalar: Ranges: iqr
    """
    # Check
    DTYPE = x.dtype
    COLMAJOR = not x.flags['C_CONTIGUOUS']
    assert DTYPE in REAL_DTYPES, f"input data type must be in {REAL_DTYPES}"
    assert x.ndim < 5, "input x must have ndim < 5"
    dim = x.ndim + axis if axis < 0 else axis
    assert 0 <= dim < x.ndim, "axis out of range [0 x.ndim)"

    # Get input/output shapes
    R = x.shape[0] if x.ndim > 0 else 1
    C = x.shape[1] if x.ndim > 1 else 1
    S = x.shape[2] if x.ndim > 2 else 1
    H = x.shape[3] if x.ndim > 3 else 1
    Ry = 1 if dim == 0 else R
    Cy = 1 if dim == 1 else C
    Sy = 1 if dim == 2 else S
    Hy = 1 if dim == 3 else H
    R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)
    DIM = c_size_t(dim)

    # Make pointer
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X = x.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((Ry, Cy, Sy, Hy), dtype=DTYPE)
    while y.ndim > x.ndim:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.iqr_s(Y, X, R, C, S, H, COLMAJOR, DIM)
    else:
        ret = CLIB.iqr_d(Y, X, R, C, S, H, COLMAJOR, DIM)
    assert ret == 0, "error during call to C function"

    return y


def idr(x, axis=0) -> np.ndarray:
    """
    Vec2scalar: Ranges: idr
    """
    # Check
    DTYPE = x.dtype
    COLMAJOR = not x.flags['C_CONTIGUOUS']
    assert DTYPE in REAL_DTYPES, f"input data type must be in {REAL_DTYPES}"
    assert x.ndim < 5, "input x must have ndim < 5"
    dim = x.ndim + axis if axis < 0 else axis
    assert 0 <= dim < x.ndim, "axis out of range [0 x.ndim)"

    # Get input/output shapes
    R = x.shape[0] if x.ndim > 0 else 1
    C = x.shape[1] if x.ndim > 1 else 1
    S = x.shape[2] if x.ndim > 2 else 1
    H = x.shape[3] if x.ndim > 3 else 1
    Ry = 1 if dim == 0 else R
    Cy = 1 if dim == 1 else C
    Sy = 1 if dim == 2 else S
    Hy = 1 if dim == 3 else H
    R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)
    DIM = c_size_t(dim)

    # Make pointer
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X = x.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((Ry, Cy, Sy, Hy), dtype=DTYPE)
    while y.ndim > x.ndim:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.idr_s(Y, X, R, C, S, H, COLMAJOR, DIM)
    else:
        ret = CLIB.idr_d(Y, X, R, C, S, H, COLMAJOR, DIM)
    assert ret == 0, "error during call to C function"

    return y


# Main
def main() -> None:
    """
    Only used for quick command-line test
    """
    x = np.random.randn(2, 3, 1, 4)
    # x = x + 1j*x
    print(x)

    tic = time()
    # y = np.max(x, axis=0) - np.min(x, axis=0)
    y = np.percentile(x, q=75, axis=0) - np.percentile(x, q=25, axis=0)
    # y = np.percentile(x, q=90, axis=0) - np.percentile(x, q=10, axis=0)
    print(time()-tic)
    print(y)

    tic = time()
    # y = range(x, axis=0)
    y = iqr(x, axis=0)
    # y = idr(x, axis=0)
    print(time()-tic)
    print(y)


if __name__ == "__main__":
    main()
