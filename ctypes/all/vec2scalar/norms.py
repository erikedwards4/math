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

__all__ = ['norm0', 'norm1', 'norm2', 'normp']

FLT_DTYPES = (np.float32, np.complex64)
REAL_DTYPES = (np.float32, np.float64)
CPLX_DTYPES = (np.complex64, np.complex128)
DTYPES = REAL_DTYPES + CPLX_DTYPES
CLIB = ctypes.cdll.LoadLibrary("libmath.so")
C_FLT_PTR = ctypes.POINTER(ctypes.c_float)
C_DBL_PTR = ctypes.POINTER(ctypes.c_double)


def norm0(x, axis=0) -> np.ndarray:
    """
    Vec2scalar: Norms: norm0
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
        ret = CLIB.norm0_s(Y, X, R, C, S, H, COLMAJOR, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.norm0_d(Y, X, R, C, S, H, COLMAJOR, DIM)
    elif DTYPE == np.complex64:
        ret = CLIB.norm0_c(Y, X, R, C, S, H, COLMAJOR, DIM)
    else:
        ret = CLIB.norm0_z(Y, X, R, C, S, H, COLMAJOR, DIM)
    assert ret == 0, "error during call to C function"

    return y


def norm1(x, axis=0) -> np.ndarray:
    """
    Vec2scalar: Norms: norm1
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
        ret = CLIB.norm1_s(Y, X, R, C, S, H, COLMAJOR, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.norm1_d(Y, X, R, C, S, H, COLMAJOR, DIM)
    elif DTYPE == np.complex64:
        ret = CLIB.norm1_c(Y, X, R, C, S, H, COLMAJOR, DIM)
    else:
        ret = CLIB.norm1_z(Y, X, R, C, S, H, COLMAJOR, DIM)
    assert ret == 0, "error during call to C function"

    return y


def norm2(x, axis=0) -> np.ndarray:
    """
    Vec2scalar: Norms: norm2
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
        ret = CLIB.norm2_s(Y, X, R, C, S, H, COLMAJOR, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.norm2_d(Y, X, R, C, S, H, COLMAJOR, DIM)
    elif DTYPE == np.complex64:
        ret = CLIB.norm2_c(Y, X, R, C, S, H, COLMAJOR, DIM)
    else:
        ret = CLIB.norm2_z(Y, X, R, C, S, H, COLMAJOR, DIM)
    assert ret == 0, "error during call to C function"

    return y


def normp(x, axis=0, p=2.0) -> np.ndarray:
    """
    Vec2scalar: Norms: normp
    """
    # Check
    DTYPE = x.dtype
    COLMAJOR = not x.flags['C_CONTIGUOUS']
    assert DTYPE in DTYPES, f"input data type must be in {DTYPES}"
    assert x.ndim < 5, "input x must have ndim < 5"
    dim = x.ndim + axis if axis < 0 else axis
    assert 0 <= dim < x.ndim, "axis out of range [0 x.ndim)"
    assert p > 0.0, "p must be positive"

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
        P = ctypes.c_float(p)
        ret = CLIB.normp_s(Y, X, R, C, S, H, COLMAJOR, DIM, P)
    elif DTYPE == np.float64:
        P = ctypes.c_double(p)
        ret = CLIB.normp_d(Y, X, R, C, S, H, COLMAJOR, DIM, P)
    elif DTYPE == np.complex64:
        P = ctypes.c_float(p)
        ret = CLIB.normp_c(Y, X, R, C, S, H, COLMAJOR, DIM, P)
    else:
        P = ctypes.c_double(p)
        ret = CLIB.normp_z(Y, X, R, C, S, H, COLMAJOR, DIM, P)
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
    y = np.linalg.norm(x, axis=0, ord=1.5)
    print(time()-tic)
    print(y)

    tic = time()
    y = normp(x, axis=0, p=1.5)
    print(time()-tic)
    print(y)


if __name__ == "__main__":
    main()
