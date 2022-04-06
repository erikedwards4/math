#!/usr/bin/env python3
"""
Makes ctypes interface to the C functions in libmath.so.
Vecs2scalar: 2 inputs, 1 output with shape reduced along dim (1 scalar per vector pair).
Each function works vector-wise (1 pair of vectors at a time).

The broadcasting only works for ndim < 3.
"""

from time import time
import ctypes
from ctypes import c_size_t
import numpy as np

__all__ = ['dist0', 'dist1', 'dist2', 'distp', 'dist_cos2']

FLT_DTYPES = (np.float32, np.complex64)
REAL_DTYPES = (np.float32, np.float64)
CPLX_DTYPES = (np.complex64, np.complex128)
DTYPES = REAL_DTYPES + CPLX_DTYPES
CLIB = ctypes.cdll.LoadLibrary("libmath.so")
C_FLT_PTR = ctypes.POINTER(ctypes.c_float)
C_DBL_PTR = ctypes.POINTER(ctypes.c_double)


def _isvec(x) -> bool:
    """
    Checks if x is a vector (only 1 axis with length > 1).
    """
    return (max(x.shape) == x.size)

def dist0(x1, x2, axis=0) -> np.ndarray:
    """
    Vecs2scalar: Distance: dist0
    """
    # Check
    DTYPE = x1.dtype
    COLMAJOR1 = not x1.flags['C_CONTIGUOUS']
    COLMAJOR2 = not x2.flags['C_CONTIGUOUS']
    assert DTYPE in DTYPES, f"input data type must be in {DTYPES}"
    assert x2.dtype == DTYPE, "inputs must have the same data type"
    assert COLMAJOR1 == COLMAJOR2, "inputs must have the same row/col order"
    assert x1.ndim < 5, "input x1 must have ndim < 5"
    assert x2.ndim < 5, "input x2 must have ndim < 5"
    NDIM = max(x1.ndim, x2.ndim)
    DIM = NDIM + axis if axis < 0 else axis
    assert 0 <= DIM < NDIM, "axis out of range [0 NDIM)"
    if x1.shape != x2.shape:
        ISVEC1 = (x1.size == max(x1.shape))
        ISVEC2 = (x2.size == max(x2.shape))
        assert ISVEC1 or ISVEC2, "one input must be a vector for broadcasting"

    # Get input/output shapes
    R1 = x1.shape[0] if x1.ndim > 0 else 1
    C1 = x1.shape[1] if x1.ndim > 1 else 1
    S1 = x1.shape[2] if x1.ndim > 2 else 1
    H1 = x1.shape[3] if x1.ndim > 3 else 1
    R2 = x2.shape[0] if x2.ndim > 0 else 1
    C2 = x2.shape[1] if x2.ndim > 1 else 1
    S2 = x2.shape[2] if x2.ndim > 2 else 1
    H2 = x2.shape[3] if x2.ndim > 3 else 1
    R = 1 if DIM == 0 else max(R1, R2)
    C = 1 if DIM == 1 else max(C1, C2)
    S = 1 if DIM == 2 else max(S1, S2)
    H = 1 if DIM == 3 else max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((R, C, S, H), dtype=DTYPE) if DTYPE in REAL_DTYPES \
        else np.empty((R, C, S, H), dtype=np.float32) if DTYPE == np.complex64 \
        else np.empty((R, C, S, H), dtype=np.float64)
    while y.ndim > NDIM:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)
    R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.dist0_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.dist0_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.complex64:
        ret = CLIB.dist0_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    else:
        ret = CLIB.dist0_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    assert ret == 0, "error during call to C function"

    return y


def dist1(x1, x2, axis=0) -> np.ndarray:
    """
    Vecs2scalar: Distance: dist1
    """
    # Check
    DTYPE = x1.dtype
    COLMAJOR1 = not x1.flags['C_CONTIGUOUS']
    COLMAJOR2 = not x2.flags['C_CONTIGUOUS']
    assert DTYPE in DTYPES, f"input data type must be in {DTYPES}"
    assert x2.dtype == DTYPE, "inputs must have the same data type"
    assert COLMAJOR1 == COLMAJOR2, "inputs must have the same row/col order"
    assert x1.ndim < 5, "input x1 must have ndim < 5"
    assert x2.ndim < 5, "input x2 must have ndim < 5"
    NDIM = max(x1.ndim, x2.ndim)
    DIM = NDIM + axis if axis < 0 else axis
    assert 0 <= DIM < NDIM, "axis out of range [0 NDIM)"
    if x1.shape != x2.shape:
        ISVEC1 = (x1.size == max(x1.shape))
        ISVEC2 = (x2.size == max(x2.shape))
        assert ISVEC1 or ISVEC2, "one input must be a vector for broadcasting"

    # Get input/output shapes
    R1 = x1.shape[0] if x1.ndim > 0 else 1
    C1 = x1.shape[1] if x1.ndim > 1 else 1
    S1 = x1.shape[2] if x1.ndim > 2 else 1
    H1 = x1.shape[3] if x1.ndim > 3 else 1
    R2 = x2.shape[0] if x2.ndim > 0 else 1
    C2 = x2.shape[1] if x2.ndim > 1 else 1
    S2 = x2.shape[2] if x2.ndim > 2 else 1
    H2 = x2.shape[3] if x2.ndim > 3 else 1
    R = 1 if DIM == 0 else max(R1, R2)
    C = 1 if DIM == 1 else max(C1, C2)
    S = 1 if DIM == 2 else max(S1, S2)
    H = 1 if DIM == 3 else max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((R, C, S, H), dtype=DTYPE) if DTYPE in REAL_DTYPES \
        else np.empty((R, C, S, H), dtype=np.float32) if DTYPE == np.complex64 \
        else np.empty((R, C, S, H), dtype=np.float64)
    while y.ndim > NDIM:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)
    R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.dist1_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.dist1_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.complex64:
        ret = CLIB.dist1_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    else:
        ret = CLIB.dist1_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    assert ret == 0, "error during call to C function"

    return y


def dist2(x1, x2, axis=0) -> np.ndarray:
    """
    Vecs2scalar: Distance: dist2
    """
    # Check
    DTYPE = x1.dtype
    COLMAJOR1 = not x1.flags['C_CONTIGUOUS']
    COLMAJOR2 = not x2.flags['C_CONTIGUOUS']
    assert DTYPE in DTYPES, f"input data type must be in {DTYPES}"
    assert x2.dtype == DTYPE, "inputs must have the same data type"
    assert COLMAJOR1 == COLMAJOR2, "inputs must have the same row/col order"
    assert x1.ndim < 5, "input x1 must have ndim < 5"
    assert x2.ndim < 5, "input x2 must have ndim < 5"
    NDIM = max(x1.ndim, x2.ndim)
    DIM = NDIM + axis if axis < 0 else axis
    assert 0 <= DIM < NDIM, "axis out of range [0 NDIM)"
    if x1.shape != x2.shape:
        ISVEC1 = (x1.size == max(x1.shape))
        ISVEC2 = (x2.size == max(x2.shape))
        assert ISVEC1 or ISVEC2, "one input must be a vector for broadcasting"

    # Get input/output shapes
    R1 = x1.shape[0] if x1.ndim > 0 else 1
    C1 = x1.shape[1] if x1.ndim > 1 else 1
    S1 = x1.shape[2] if x1.ndim > 2 else 1
    H1 = x1.shape[3] if x1.ndim > 3 else 1
    R2 = x2.shape[0] if x2.ndim > 0 else 1
    C2 = x2.shape[1] if x2.ndim > 1 else 1
    S2 = x2.shape[2] if x2.ndim > 2 else 1
    H2 = x2.shape[3] if x2.ndim > 3 else 1
    R = 1 if DIM == 0 else max(R1, R2)
    C = 1 if DIM == 1 else max(C1, C2)
    S = 1 if DIM == 2 else max(S1, S2)
    H = 1 if DIM == 3 else max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((R, C, S, H), dtype=DTYPE) if DTYPE in REAL_DTYPES \
        else np.empty((R, C, S, H), dtype=np.float32) if DTYPE == np.complex64 \
        else np.empty((R, C, S, H), dtype=np.float64)
    while y.ndim > NDIM:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)
    R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.dist2_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.dist2_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.complex64:
        ret = CLIB.dist2_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    else:
        ret = CLIB.dist2_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    assert ret == 0, "error during call to C function"

    return y


def distp(x1, x2, axis=0, p=2.0) -> np.ndarray:
    """
    Vecs2scalar: Distance: distp
    """
    # Check
    DTYPE = x1.dtype
    COLMAJOR1 = not x1.flags['C_CONTIGUOUS']
    COLMAJOR2 = not x2.flags['C_CONTIGUOUS']
    assert DTYPE in DTYPES, f"input data type must be in {DTYPES}"
    assert x2.dtype == DTYPE, "inputs must have the same data type"
    assert COLMAJOR1 == COLMAJOR2, "inputs must have the same row/col order"
    assert x1.ndim < 5, "input x1 must have ndim < 5"
    assert x2.ndim < 5, "input x2 must have ndim < 5"
    NDIM = max(x1.ndim, x2.ndim)
    DIM = NDIM + axis if axis < 0 else axis
    assert 0 <= DIM < NDIM, "axis out of range [0 NDIM)"
    assert p > 0.0, "p (power) must be positive"
    if x1.shape != x2.shape:
        ISVEC1 = (x1.size == max(x1.shape))
        ISVEC2 = (x2.size == max(x2.shape))
        assert ISVEC1 or ISVEC2, "one input must be a vector for broadcasting"

    # Get input/output shapes
    R1 = x1.shape[0] if x1.ndim > 0 else 1
    C1 = x1.shape[1] if x1.ndim > 1 else 1
    S1 = x1.shape[2] if x1.ndim > 2 else 1
    H1 = x1.shape[3] if x1.ndim > 3 else 1
    R2 = x2.shape[0] if x2.ndim > 0 else 1
    C2 = x2.shape[1] if x2.ndim > 1 else 1
    S2 = x2.shape[2] if x2.ndim > 2 else 1
    H2 = x2.shape[3] if x2.ndim > 3 else 1
    R = 1 if DIM == 0 else max(R1, R2)
    C = 1 if DIM == 1 else max(C1, C2)
    S = 1 if DIM == 2 else max(S1, S2)
    H = 1 if DIM == 3 else max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((R, C, S, H), dtype=DTYPE) if DTYPE in REAL_DTYPES \
        else np.empty((R, C, S, H), dtype=np.float32) if DTYPE == np.complex64 \
        else np.empty((R, C, S, H), dtype=np.float64)
    while y.ndim > NDIM:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)
    R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)

    # Run
    if DTYPE == np.float32:
        P = ctypes.c_float(p)
        ret = CLIB.distp_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM, P)
    elif DTYPE == np.float64:
        P = ctypes.c_double(p)
        ret = CLIB.distp_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM, P)
    elif DTYPE == np.complex64:
        P = ctypes.c_float(p)
        ret = CLIB.distp_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM, P)
    else:
        P = ctypes.c_double(p)
        ret = CLIB.distp_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM, P)
    assert ret == 0, "error during call to C function"

    return y


def dist_cos2(x1, x2, axis=0) -> np.ndarray:
    """
    Vecs2scalar: Distance: dist_cos2
    """
    # Check
    DTYPE = x1.dtype
    COLMAJOR1 = not x1.flags['C_CONTIGUOUS']
    COLMAJOR2 = not x2.flags['C_CONTIGUOUS']
    assert DTYPE in REAL_DTYPES, f"input data type must be in {REAL_DTYPES}"
    assert x2.dtype == DTYPE, "inputs must have the same data type"
    assert COLMAJOR1 == COLMAJOR2, "inputs must have the same row/col order"
    assert x1.ndim < 5, "input x1 must have ndim < 5"
    assert x2.ndim < 5, "input x2 must have ndim < 5"
    NDIM = max(x1.ndim, x2.ndim)
    DIM = NDIM + axis if axis < 0 else axis
    assert 0 <= DIM < NDIM, "axis out of range [0 NDIM)"
    if x1.shape != x2.shape:
        ISVEC1 = (x1.size == max(x1.shape))
        ISVEC2 = (x2.size == max(x2.shape))
        assert ISVEC1 or ISVEC2, "one input must be a vector for broadcasting"

    # Get input/output shapes
    R1 = x1.shape[0] if x1.ndim > 0 else 1
    C1 = x1.shape[1] if x1.ndim > 1 else 1
    S1 = x1.shape[2] if x1.ndim > 2 else 1
    H1 = x1.shape[3] if x1.ndim > 3 else 1
    R2 = x2.shape[0] if x2.ndim > 0 else 1
    C2 = x2.shape[1] if x2.ndim > 1 else 1
    S2 = x2.shape[2] if x2.ndim > 2 else 1
    H2 = x2.shape[3] if x2.ndim > 3 else 1
    R = 1 if DIM == 0 else max(R1, R2)
    C = 1 if DIM == 1 else max(C1, C2)
    S = 1 if DIM == 2 else max(S1, S2)
    H = 1 if DIM == 3 else max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((R, C, S, H), dtype=DTYPE)
    while y.ndim > NDIM:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)
    R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.dist_cos2_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    else:
        ret = CLIB.dist_cos2_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    assert ret == 0, "error during call to C function"

    return y


# Main
def main() -> None:
    """
    Only used for quick command-line test
    """
    x1 = np.random.randn(2, 3, 1)
    x2 = np.random.randn(2, 1, 1)
    x1 = x1 + 1j*x1
    x2 = x2 + 1j*x2
    # print(f"x1 = {x1}")
    # print(f"x2 = {x2}")
    p = 0.5

    tic = time()
    y = np.power(np.sum(np.abs(x1-x2)**p, axis=0), 1.0/p)
    print(time()-tic)
    print(f"y = {y}")

    tic = time()
    y = distp(x1, x2, axis=0, p=p)
    print(time()-tic)
    print(f"y = {y}")


if __name__ == "__main__":
    main()
