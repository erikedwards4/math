#!/usr/bin/env python3
"""
Makes ctypes interface to the C functions in libmath.so.
Elementwise2: 2 inputs, 1 output with same shape as input
Each function works element-wise (1 pair of elements at a time).
"""

from time import time
import ctypes
from ctypes import c_size_t
import numpy as np

__all__ = ['hypot', 'atan2']

FLT_DTYPES = (np.float32, np.complex64)
REAL_DTYPES = (np.float32, np.float64)
CPLX_DTYPES = (np.complex64, np.complex128)
DTYPES = REAL_DTYPES + CPLX_DTYPES
CLIB = ctypes.cdll.LoadLibrary("libmath.so")
C_FLT_PTR = ctypes.POINTER(ctypes.c_float)
C_DBL_PTR = ctypes.POINTER(ctypes.c_double)


def atan2(x1, x2, inplace=False) -> np.ndarray:
    """
    Elementwise2: Trig2: atan2
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
            ret = CLIB.atan2_inplace_s(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.float64:
            ret = CLIB.atan2_inplace_d(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        assert ret == 0, "error during call to C function"
        return x1
    else:
        y = np.empty((R, C, S, H), dtype=DTYPE)
        while y.ndim > max(x1.ndim, x2.ndim):
            y = y.squeeze(axis=-1)
        Y = y.ctypes.data_as(C_PTR)
        R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)
        if DTYPE == np.float32:
            ret = CLIB.atan2_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.float64:
            ret = CLIB.atan2_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        assert ret == 0, "error during call to C function"
        return y


def hypot(x1, x2, inplace=False) -> np.ndarray:
    """
    Elementwise2: Trig2: hypot
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
            ret = CLIB.hypot_inplace_s(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        else:
            ret = CLIB.hypot_inplace_d(X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        assert ret == 0, "error during call to C function"
        return x1
    else:
        y = np.empty((R, C, S, H), dtype=DTYPE) if DTYPE in REAL_DTYPES \
            else np.empty((R, C, S, H), dtype=np.float32) if DTYPE == np.complex64 \
            else np.empty((R, C, S, H), dtype=np.float64)
        while y.ndim > max(x1.ndim, x2.ndim):
            y = y.squeeze(axis=-1)
        Y = y.ctypes.data_as(C_PTR)
        R, C, S, H = c_size_t(R), c_size_t(C), c_size_t(S), c_size_t(H)
        if DTYPE == np.float32:
            ret = CLIB.hypot_s(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.float64:
            ret = CLIB.hypot_d(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        elif DTYPE == np.complex64:
            ret = CLIB.hypot_c(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
        else:
            ret = CLIB.hypot_z(Y, X1, X2, R1, C1, S1, H1, R2, C2, S2, H2, COLMAJOR1)
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
    y = np.arctan2(x1, x2)
    # y = np.hypot(x1, x2)
    # y = np.hypot(np.abs(x1), np.abs(x2))
    print(time()-tic)
    print(y)

    tic = time()
    # y = hypot(x1, x2, inplace=False)
    # hypot(x1, x2, inplace=True)
    y = atan2(x1, x2, inplace=False)
    atan2(x1, x2, inplace=True)
    print(time()-tic)
    print(y)
    print(x1)


if __name__ == "__main__":
    main()
