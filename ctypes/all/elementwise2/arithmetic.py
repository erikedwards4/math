#!/usr/bin/env python3
"""
Makes ctypes interface to the C functions in libmath.so.
Elementwise2: Arithmetic: 2 inputs, 1 output with same shape as input.
Each function works element-wise (1 pair of elements at a time).
"""

from time import time
import ctypes
from ctypes import c_size_t
import numpy as np

__all__ = ['plus', 'minus', 'times', 'rdivide', 'adiff', 'pow']

FLT_DTYPES = (np.float32, np.complex64)
REAL_DTYPES = (np.float32, np.float64)
CPLX_DTYPES = (np.complex64, np.complex128)
DTYPES = REAL_DTYPES + CPLX_DTYPES
CLIB = ctypes.cdll.LoadLibrary("libmath.so")
C_FLT_PTR = ctypes.POINTER(ctypes.c_float)
C_DBL_PTR = ctypes.POINTER(ctypes.c_double)


def plus(x1, x2, inplace=False) -> np.ndarray:
    """
    Elementwise2: Arithmetic: plus
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

    # Get input/output shapes
    R1 = x1.shape[0] if x1.ndim > 0 else 1
    C1 = x1.shape[1] if x1.ndim > 1 else 1
    S1 = x1.shape[2] if x1.ndim > 2 else 1
    H1 = x1.shape[3] if x1.ndim > 3 else 1
    R2 = x2.shape[0] if x2.ndim > 0 else 1
    C2 = x2.shape[1] if x2.ndim > 1 else 1
    S2 = x2.shape[2] if x2.ndim > 2 else 1
    H2 = x2.shape[3] if x2.ndim > 3 else 1
    R, C, S, H = max(R1, R2), max(C1, C2), max(S1, S2), max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)
    NDIM = max(x1.ndim, x2.ndim)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Run
    if inplace:
        if DTYPE == np.float32:
            ret = CLIB.plus_inplace_s(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.float64:
            ret = CLIB.plus_inplace_d(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.complex64:
            ret = CLIB.plus_inplace_c(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        else:
            ret = CLIB.plus_inplace_z(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        assert ret == 0, "error during call to C function"
        return x1
    else:
        y = np.empty((R, C, S, H), dtype=DTYPE)
        while y.ndim > NDIM:
            y = y.squeeze(axis=-1)
        Y = y.ctypes.data_as(C_PTR)
        R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)
        if DTYPE == np.float32:
            ret = CLIB.plus_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.float64:
            ret = CLIB.plus_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.complex64:
            ret = CLIB.plus_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        else:
            ret = CLIB.plus_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        assert ret == 0, "error during call to C function"
        return y


def minus(x1, x2, inplace=False) -> np.ndarray:
    """
    Elementwise2: Arithmetic: minus
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

    # Get input/output shapes
    R1 = x1.shape[0] if x1.ndim > 0 else 1
    C1 = x1.shape[1] if x1.ndim > 1 else 1
    S1 = x1.shape[2] if x1.ndim > 2 else 1
    H1 = x1.shape[3] if x1.ndim > 3 else 1
    R2 = x2.shape[0] if x2.ndim > 0 else 1
    C2 = x2.shape[1] if x2.ndim > 1 else 1
    S2 = x2.shape[2] if x2.ndim > 2 else 1
    H2 = x2.shape[3] if x2.ndim > 3 else 1
    R, C, S, H = max(R1, R2), max(C1, C2), max(S1, S2), max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Run
    if inplace:
        if DTYPE == np.float32:
            ret = CLIB.minus_inplace_s(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.float64:
            ret = CLIB.minus_inplace_d(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.complex64:
            ret = CLIB.minus_inplace_c(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        else:
            ret = CLIB.minus_inplace_z(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        assert ret == 0, "error during call to C function"
        return x1
    else:
        y = np.empty((R, C, S, H), dtype=DTYPE)
        while y.ndim > NDIM:
            y = y.squeeze(axis=-1)
        Y = y.ctypes.data_as(C_PTR)
        R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)
        if DTYPE == np.float32:
            ret = CLIB.minus_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.float64:
            ret = CLIB.minus_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.complex64:
            ret = CLIB.minus_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        else:
            ret = CLIB.minus_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        assert ret == 0, "error during call to C function"
        return y


def times(x1, x2, inplace=False) -> np.ndarray:
    """
    Elementwise2: Arithmetic: times
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

    # Get input/output shapes
    R1 = x1.shape[0] if x1.ndim > 0 else 1
    C1 = x1.shape[1] if x1.ndim > 1 else 1
    S1 = x1.shape[2] if x1.ndim > 2 else 1
    H1 = x1.shape[3] if x1.ndim > 3 else 1
    R2 = x2.shape[0] if x2.ndim > 0 else 1
    C2 = x2.shape[1] if x2.ndim > 1 else 1
    S2 = x2.shape[2] if x2.ndim > 2 else 1
    H2 = x2.shape[3] if x2.ndim > 3 else 1
    R, C, S, H = max(R1, R2), max(C1, C2), max(S1, S2), max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Run
    if inplace:
        if DTYPE == np.float32:
            ret = CLIB.times_inplace_s(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.float64:
            ret = CLIB.times_inplace_d(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.complex64:
            ret = CLIB.times_inplace_c(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        else:
            ret = CLIB.times_inplace_z(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        assert ret == 0, "error during call to C function"
        return x1
    else:
        y = np.empty((R, C, S, H), dtype=DTYPE)
        while y.ndim > NDIM:
            y = y.squeeze(axis=-1)
        Y = y.ctypes.data_as(C_PTR)
        R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)
        if DTYPE == np.float32:
            ret = CLIB.times_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.float64:
            ret = CLIB.times_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.complex64:
            ret = CLIB.times_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        else:
            ret = CLIB.times_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        assert ret == 0, "error during call to C function"
        return y


def rdivide(x1, x2, inplace=False) -> np.ndarray:
    """
    Elementwise2: Arithmetic: rdivide
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

    # Get input/output shapes
    R1 = x1.shape[0] if x1.ndim > 0 else 1
    C1 = x1.shape[1] if x1.ndim > 1 else 1
    S1 = x1.shape[2] if x1.ndim > 2 else 1
    H1 = x1.shape[3] if x1.ndim > 3 else 1
    R2 = x2.shape[0] if x2.ndim > 0 else 1
    C2 = x2.shape[1] if x2.ndim > 1 else 1
    S2 = x2.shape[2] if x2.ndim > 2 else 1
    H2 = x2.shape[3] if x2.ndim > 3 else 1
    R, C, S, H = max(R1, R2), max(C1, C2), max(S1, S2), max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Run
    if inplace:
        if DTYPE == np.float32:
            ret = CLIB.rdivide_inplace_s(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.float64:
            ret = CLIB.rdivide_inplace_d(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.complex64:
            ret = CLIB.rdivide_inplace_c(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        else:
            ret = CLIB.rdivide_inplace_z(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        assert ret == 0, "error during call to C function"
        return x1
    else:
        y = np.empty((R, C, S, H), dtype=DTYPE)
        while y.ndim > NDIM:
            y = y.squeeze(axis=-1)
        Y = y.ctypes.data_as(C_PTR)
        R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)
        if DTYPE == np.float32:
            ret = CLIB.rdivide_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.float64:
            ret = CLIB.rdivide_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.complex64:
            ret = CLIB.rdivide_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        else:
            ret = CLIB.rdivide_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        assert ret == 0, "error during call to C function"
        return y


def adiff(x1, x2, inplace=False) -> np.ndarray:
    """
    Elementwise2: Arithmetic: adiff
    """
    # Check
    DTYPE = x1.dtype
    COLMAJOR1 = not x1.flags['C_CONTIGUOUS']
    COLMAJOR2 = not x2.flags['C_CONTIGUOUS']
    if inplace:
        assert DTYPE in REAL_DTYPES, f"input data type must be in {REAL_DTYPES}"
    else:
        assert DTYPE in DTYPES, f"input data type must be in {DTYPES}"
    assert x2.dtype == DTYPE, "inputs must have the same data type"
    assert COLMAJOR1 == COLMAJOR2, "inputs must have the same row/col order"
    assert x1.ndim < 5, "input x1 must have ndim < 5"
    assert x2.ndim < 5, "input x2 must have ndim < 5"
    NDIM = max(x1.ndim, x2.ndim)

    # Get input/output shapes
    R1 = x1.shape[0] if x1.ndim > 0 else 1
    C1 = x1.shape[1] if x1.ndim > 1 else 1
    S1 = x1.shape[2] if x1.ndim > 2 else 1
    H1 = x1.shape[3] if x1.ndim > 3 else 1
    R2 = x2.shape[0] if x2.ndim > 0 else 1
    C2 = x2.shape[1] if x2.ndim > 1 else 1
    S2 = x2.shape[2] if x2.ndim > 2 else 1
    H2 = x2.shape[3] if x2.ndim > 3 else 1
    R, C, S, H = max(R1, R2), max(C1, C2), max(S1, S2), max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Run
    if inplace:
        if DTYPE == np.float32:
            ret = CLIB.adiff_inplace_s(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        else:
            ret = CLIB.adiff_inplace_d(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        assert ret == 0, "error during call to C function"
        return x1
    else:
        y = np.empty((R, C, S, H), dtype=DTYPE) if DTYPE in REAL_DTYPES \
            else np.empty((R, C, S, H), dtype=np.float32) if DTYPE == np.complex64 \
            else np.empty((R, C, S, H), dtype=np.float64)
        while y.ndim > NDIM:
            y = y.squeeze(axis=-1)
        Y = y.ctypes.data_as(C_PTR)
        R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)
        if DTYPE == np.float32:
            ret = CLIB.adiff_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.float64:
            ret = CLIB.adiff_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.complex64:
            ret = CLIB.adiff_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        else:
            ret = CLIB.adiff_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        assert ret == 0, "error during call to C function"
        return y


def pow(x1, x2, inplace=False) -> np.ndarray:
    """
    Elementwise2: Arithmetic: pow
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

    # Get input/output shapes
    R1 = x1.shape[0] if x1.ndim > 0 else 1
    C1 = x1.shape[1] if x1.ndim > 1 else 1
    S1 = x1.shape[2] if x1.ndim > 2 else 1
    H1 = x1.shape[3] if x1.ndim > 3 else 1
    R2 = x2.shape[0] if x2.ndim > 0 else 1
    C2 = x2.shape[1] if x2.ndim > 1 else 1
    S2 = x2.shape[2] if x2.ndim > 2 else 1
    H2 = x2.shape[3] if x2.ndim > 3 else 1
    R, C, S, H = max(R1, R2), max(C1, C2), max(S1, S2), max(H1, H2)
    R1, C1, S1, H1 = c_size_t(R1), c_size_t(C1), c_size_t(S1), c_size_t(H1)
    R2, C2, S2, H2 = c_size_t(R2), c_size_t(C2), c_size_t(S2), c_size_t(H2)

    # Make pointers
    C_PTR = C_FLT_PTR if DTYPE in FLT_DTYPES else C_DBL_PTR
    X1 = x1.ctypes.data_as(C_PTR)
    X2 = x2.ctypes.data_as(C_PTR)

    # Run
    if inplace:
        if DTYPE == np.float32:
            ret = CLIB.pow_inplace_s(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.float64:
            ret = CLIB.pow_inplace_d(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.complex64:
            ret = CLIB.pow_inplace_c(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        else:
            ret = CLIB.pow_inplace_z(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        assert ret == 0, "error during call to C function"
        return x1
    else:
        y = np.empty((R, C, S, H), dtype=DTYPE)
        while y.ndim > NDIM:
            y = y.squeeze(axis=-1)
        Y = y.ctypes.data_as(C_PTR)
        R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)
        if DTYPE == np.float32:
            ret = CLIB.pow_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.float64:
            ret = CLIB.pow_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.complex64:
            ret = CLIB.pow_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        else:
            ret = CLIB.pow_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        assert ret == 0, "error during call to C function"
        return y


# Main
def main() -> None:
    """
    Only used for quick command-line test
    """
    x1 = np.random.randn(2, 3)
    # x1 = x1 + 1j*x1
    x2 = np.random.randn(2, 3)
    # x2 = x2 + 1j*x2

    tic = time()
    y = np.abs(np.subtract(x1, x2))
    print(time()-tic)
    print(y)

    tic = time()
    # rdivide(x1, x2, inplace=True)
    y = adiff(x1, x2, inplace=False)
    adiff(x1, x2, inplace=True)
    print(time()-tic)
    print(y)


if __name__ == "__main__":
    main()
