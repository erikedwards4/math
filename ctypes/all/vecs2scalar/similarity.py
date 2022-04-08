#!/usr/bin/env python3
"""
Makes ctypes interface to the C functions in libmath.so.
Vecs2scalar: 2 inputs, 1 output with shape reduced along dim (1 scalar per vector pair).
Each function works vector-wise (1 pair of vectors at a time).

The broadcasting is only intended for a single vector into a larger tensor.
"""

from time import time
import ctypes
from ctypes import c_size_t
import numpy as np

__all__ = ['dot', 'cov', 'corr', 'corr_opus', 'cos2', 'cokurtosis', 'spearman', 'kendall']

FLT_DTYPES = (np.float32, np.complex64)
REAL_DTYPES = (np.float32, np.float64)
CPLX_DTYPES = (np.complex64, np.complex128)
DTYPES = REAL_DTYPES + CPLX_DTYPES
CLIB = ctypes.cdll.LoadLibrary("libmath.so")
C_FLT_PTR = ctypes.POINTER(ctypes.c_float)
C_DBL_PTR = ctypes.POINTER(ctypes.c_double)


def dot(x1, x2, axis=0) -> np.ndarray:
    """
    Vecs2scalar: Similarity: dot
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
    dim = NDIM + axis if axis < 0 else axis
    assert 0 <= dim < NDIM, "axis out of range [0 NDIM)"
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
    R = 1 if dim == 0 else max(R1, R2)
    C = 1 if dim == 1 else max(C1, C2)
    S = 1 if dim == 2 else max(S1, S2)
    H = 1 if dim == 3 else max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)
    DIM = c_size_t(dim)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((R, C, S, H), dtype=DTYPE)
    while y.ndim > NDIM:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.dot_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.dot_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.complex64:
        ret = CLIB.dot_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    else:
        ret = CLIB.dot_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    assert ret == 0, "error during call to C function"

    return y


def cov(x1, x2, axis=0) -> np.ndarray:
    """
    Vecs2scalar: Similarity: cov
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
    dim = NDIM + axis if axis < 0 else axis
    assert 0 <= dim < NDIM, "axis out of range [0 NDIM)"
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
    R = 1 if dim == 0 else max(R1, R2)
    C = 1 if dim == 1 else max(C1, C2)
    S = 1 if dim == 2 else max(S1, S2)
    H = 1 if dim == 3 else max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)
    DIM = c_size_t(dim)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((R, C, S, H), dtype=DTYPE)
    while y.ndim > NDIM:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.cov_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.cov_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.complex64:
        ret = CLIB.cov_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    else:
        ret = CLIB.cov_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    assert ret == 0, "error during call to C function"

    return y


def corr(x1, x2, axis=0) -> np.ndarray:
    """
    Vecs2scalar: Similarity: corr
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
    dim = NDIM + axis if axis < 0 else axis
    assert 0 <= dim < NDIM, "axis out of range [0 NDIM)"
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
    R = 1 if dim == 0 else max(R1, R2)
    C = 1 if dim == 1 else max(C1, C2)
    S = 1 if dim == 2 else max(S1, S2)
    H = 1 if dim == 3 else max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)
    DIM = c_size_t(dim)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((R, C, S, H), dtype=DTYPE)
    while y.ndim > NDIM:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.corr_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.corr_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.complex64:
        ret = CLIB.corr_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    else:
        ret = CLIB.corr_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    assert ret == 0, "error during call to C function"

    return y


def corr_opus(x1, x2, axis=0) -> np.ndarray:
    """
    Vecs2scalar: Similarity: corr_opus
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
    dim = NDIM + axis if axis < 0 else axis
    assert 0 <= dim < NDIM, "axis out of range [0 NDIM)"
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
    R = 1 if dim == 0 else max(R1, R2)
    C = 1 if dim == 1 else max(C1, C2)
    S = 1 if dim == 2 else max(S1, S2)
    H = 1 if dim == 3 else max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)
    DIM = c_size_t(dim)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((R, C, S, H), dtype=DTYPE)
    while y.ndim > NDIM:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.corr_opus_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.corr_opus_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.complex64:
        ret = CLIB.corr_opus_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    else:
        ret = CLIB.corr_opus_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    assert ret == 0, "error during call to C function"

    return y


def cos2(x1, x2, axis=0) -> np.ndarray:
    """
    Vecs2scalar: Similarity: cos2
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
    dim = NDIM + axis if axis < 0 else axis
    assert 0 <= dim < NDIM, "axis out of range [0 NDIM)"
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
    R = 1 if dim == 0 else max(R1, R2)
    C = 1 if dim == 1 else max(C1, C2)
    S = 1 if dim == 2 else max(S1, S2)
    H = 1 if dim == 3 else max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)
    DIM = c_size_t(dim)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((R, C, S, H), dtype=DTYPE)
    while y.ndim > NDIM:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.cos2_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.cos2_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.complex64:
        ret = CLIB.cos2_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    else:
        ret = CLIB.cos2_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    assert ret == 0, "error during call to C function"

    return y


def cokurtosis(x1, x2, axis=0) -> np.ndarray:
    """
    Vecs2scalar: Similarity: cokurtosis
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
    dim = NDIM + axis if axis < 0 else axis
    assert 0 <= dim < NDIM, "axis out of range [0 NDIM)"
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
    R = 1 if dim == 0 else max(R1, R2)
    C = 1 if dim == 1 else max(C1, C2)
    S = 1 if dim == 2 else max(S1, S2)
    H = 1 if dim == 3 else max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)
    DIM = c_size_t(dim)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((R, C, S, H), dtype=DTYPE)
    while y.ndim > NDIM:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.cokurtosis_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.cokurtosis_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    assert ret == 0, "error during call to C function"

    return y


def spearman(x1, x2, axis=0) -> np.ndarray:
    """
    Vecs2scalar: Similarity: spearman
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
    dim = NDIM + axis if axis < 0 else axis
    assert 0 <= dim < NDIM, "axis out of range [0 NDIM)"
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
    R = 1 if dim == 0 else max(R1, R2)
    C = 1 if dim == 1 else max(C1, C2)
    S = 1 if dim == 2 else max(S1, S2)
    H = 1 if dim == 3 else max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)
    DIM = c_size_t(dim)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((R, C, S, H), dtype=DTYPE)
    while y.ndim > NDIM:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.spearman_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.spearman_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    assert ret == 0, "error during call to C function"

    return y


def kendall(x1, x2, axis=0) -> np.ndarray:
    """
    Vecs2scalar: Similarity: kendall
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
    dim = NDIM + axis if axis < 0 else axis
    assert 0 <= dim < NDIM, "axis out of range [0 NDIM)"
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
    R = 1 if dim == 0 else max(R1, R2)
    C = 1 if dim == 1 else max(C1, C2)
    S = 1 if dim == 2 else max(S1, S2)
    H = 1 if dim == 3 else max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)
    DIM = c_size_t(dim)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Init output
    y = np.empty((R, C, S, H), dtype=DTYPE)
    while y.ndim > NDIM:
        y = y.squeeze(axis=-1)
    Y = y.ctypes.data_as(C_PTR)

    # Run
    if DTYPE == np.float32:
        ret = CLIB.kendall_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    elif DTYPE == np.float64:
        ret = CLIB.kendall_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1, DIM)
    assert ret == 0, "error during call to C function"

    return y


# Main
def main() -> None:
    """
    Only used for quick command-line test
    """
    x1 = np.random.randn(1, 1, 1, 2)
    x2 = np.random.randn(4, 3, 2, 2)
    # x1 = x1 + 1j*x1
    # x2 = x2 + 1j*x2

    tic = time()
    # y = np.tensordot(x1, x2, axes=(0, 0))
    print(time()-tic)
    # print(y)

    tic = time()
    y = dot(x1, x2, axis=3)
    print(time()-tic)
    print(y)
    y = dot(x1, x2, axis=-1)
    print(y)


if __name__ == "__main__":
    main()
