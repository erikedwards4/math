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

__all__ = ['iprctile', 'imed', 'imax', 'imin', 'iamax']

FLT_DTYPES = (np.float32, np.complex64)
REAL_DTYPES = (np.float32, np.float64)
CPLX_DTYPES = (np.complex64, np.complex128)
DTYPES = REAL_DTYPES + CPLX_DTYPES
CLIB = ctypes.cdll.LoadLibrary("libmath.so")
C_FLT_PTR = ctypes.POINTER(ctypes.c_float)
C_DBL_PTR = ctypes.POINTER(ctypes.c_double)


def iprctile(x, axis=0, p=50.0) -> np.ndarray:
    """
    Vec2scalar: Iprctiles: prctile
    """
    # Check
    DTYPE = x.dtype
    COLMAJOR = not x.flags['C_CONTIGUOUS']
    assert DTYPE in REAL_DTYPES, f"input data type must be in {REAL_DTYPES}"
    assert x.ndim < 5, "input x must have ndim < 5"
    dim = x.ndim + axis if axis < 0 else axis
    assert 0 <= dim < x.ndim, "axis out of range [0 x.ndim)"
    assert 0.0 <= p <= 100.0, "p (percentile) must be in [0 100]"

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
        P = ctypes.c_float(p)
        ret = CLIB.iprctile_s(Y, X, R, C, S, H, COLMAJOR, DIM, P)
    else:
        P = ctypes.c_double(p)
        ret = CLIB.iprctile_d(Y, X, R, C, S, H, COLMAJOR, DIM, P)
    assert ret == 0, "error during call to C function"

    return y


def imed(x, axis=0) -> np.ndarray:
    """
    Vec2scalar: Iprctiles: imed
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
        ret = CLIB.imed_s(Y, X, R, C, S, H, COLMAJOR, DIM)
    else:
        ret = CLIB.imed_d(Y, X, R, C, S, H, COLMAJOR, DIM)
    assert ret == 0, "error during call to C function"

    return y


def imax(x, axis=0) -> np.ndarray:
    """
    Vec2scalar: Iprctiles: imax
    """
    # Check
    DTYPE = x.dtype
    COLMAJOR = not x.flags['C_CONTIGUOUS']
    assert DTYPE in DTYPES, f"input data type must be in {DTYPES}"
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
    y = np.empty((Ry, Cy, Sy, Hy), dtype=DTYPE) if DTYPE in REAL_DTYPES \
        else np.empty((Ry, Cy, Sy, Hy), dtype=np.float32) if DTYPE == np.complex64 \
        else np.empty((Ry, Cy, Sy, Hy), dtype=np.float64)
    while y.ndim > x.ndim:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.imax_s(Y, X, R, C, S, H, COLMAJOR, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.imax_d(Y, X, R, C, S, H, COLMAJOR, DIM)
    elif DTYPE == np.complex64:
        ret = CLIB.imax_c(Y, X, R, C, S, H, COLMAJOR, DIM)
    else:
        ret = CLIB.imax_z(Y, X, R, C, S, H, COLMAJOR, DIM)
    assert ret == 0, "error during call to C function"

    return y


def imin(x, axis=0) -> np.ndarray:
    """
    Vec2scalar: Iprctiles: imin
    """
    # Check
    DTYPE = x.dtype
    COLMAJOR = not x.flags['C_CONTIGUOUS']
    assert DTYPE in DTYPES, f"input data type must be in {DTYPES}"
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
    y = np.empty((Ry, Cy, Sy, Hy), dtype=DTYPE) if DTYPE in REAL_DTYPES \
        else np.empty((Ry, Cy, Sy, Hy), dtype=np.float32) if DTYPE == np.complex64 \
        else np.empty((Ry, Cy, Sy, Hy), dtype=np.float64)
    while y.ndim > x.ndim:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.imin_s(Y, X, R, C, S, H, COLMAJOR, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.imin_d(Y, X, R, C, S, H, COLMAJOR, DIM)
    elif DTYPE == np.complex64:
        ret = CLIB.imin_c(Y, X, R, C, S, H, COLMAJOR, DIM)
    else:
        ret = CLIB.imin_z(Y, X, R, C, S, H, COLMAJOR, DIM)
    assert ret == 0, "error during call to C function"

    return y


def iamax(x, axis=0) -> np.ndarray:
    """
    Vec2scalar: Iprctiles: iamax
    """
    # Check
    DTYPE = x.dtype
    COLMAJOR = not x.flags['C_CONTIGUOUS']
    assert DTYPE in DTYPES, f"input data type must be in {DTYPES}"
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
    y = np.empty((Ry, Cy, Sy, Hy), dtype=DTYPE) if DTYPE in REAL_DTYPES \
        else np.empty((Ry, Cy, Sy, Hy), dtype=np.float32) if DTYPE == np.complex64 \
        else np.empty((Ry, Cy, Sy, Hy), dtype=np.float64)
    while y.ndim > x.ndim:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.iamax_s(Y, X, R, C, S, H, COLMAJOR, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.iamax_d(Y, X, R, C, S, H, COLMAJOR, DIM)
    elif DTYPE == np.complex64:
        ret = CLIB.iamax_c(Y, X, R, C, S, H, COLMAJOR, DIM)
    else:
        ret = CLIB.iamax_z(Y, X, R, C, S, H, COLMAJOR, DIM)
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
    # y = np.max(np.abs(x), axis=0)
    y = np.percentile(x, axis=0, q=50.0)
    print(time()-tic)
    print(y)

    tic = time()
    # y = iamax(x, axis=0)
    y = iprctile(x, axis=3, p=50.0)
    print(time()-tic)
    print(y)


if __name__ == "__main__":
    main()
