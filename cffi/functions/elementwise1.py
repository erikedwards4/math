#!/usr/bin/env python3
"""
Makes cffi interface to the C functions.
Elementwise1: 1 input, 1 output with same shape as input
Each function works element-wise (1 element at a time).
"""

from time import time
import numpy as np
from _math_cffi import ffi, lib

__all__ = ['abs']

# Global constants
DTYPES = (np.float32, np.float64, float, np.complex64, np.complex128, complex)


# Elementwise1: Nonlin:

def abs(x) -> np.ndarray:
    """
    Elementwise1: Nonlin: abs
    """
    assert isinstance(x, np.ndarray), "input must be an nd.array"
    assert x.dtype in DTYPES, f"input data type must be in {DTYPES}"

    if x.dtype == np.float32:
        X = ffi.cast('float*', x.ctypes.data)
        ret = lib.abs_inplace_s(X, x.size)
        assert ret == 0, "error during call to C function"
        return x
    elif x.dtype == np.float64:
        X = ffi.cast('double*', x.ctypes.data)
        ret = lib.abs_inplace_d(X, x.size)
        assert ret == 0, "error during call to C function"
        return x
    elif x.dtype == np.complex64:
        y = np.empty_like(x, dtype=np.float32)
        X = ffi.cast('float*', x.ctypes.data)
        Y = ffi.cast('float*', y.ctypes.data)
        ret = lib.abs_c(Y, X, x.size)
        assert ret == 0, "error during call to C function"
        return y
    else:
        y = np.empty_like(x, dtype=np.float64)
        X = ffi.cast('double*', x.ctypes.data)
        Y = ffi.cast('double*', y.ctypes.data)
        ret = lib.abs_z(Y, X, x.size)
        assert ret == 0, "error during call to C function"
        return y


# Main
def main() -> np.ndarray:
    """
    Only used for quick command-line test
    """
    x = np.random.randn(20000, 300)
    x = x + 1j*x
    tic = time()
    y = abs(x)
    print(time()-tic)
    tic = time()
    y = np.abs(x)
    print(time()-tic)
    return y


if __name__ == "__main__":
    main()
