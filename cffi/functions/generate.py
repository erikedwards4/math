#!/usr/bin/env python3
"""
Makes cffi interface to the C functions.
Generate: 0 inputs, 1 output
"""

from time import time
import numpy as np
from _math_cffi import ffi, lib

__all__ = ['zeros', 'ones', 'twos']

# Global constants
DTYPES = (np.float32, np.float64, float, np.complex64, np.complex128, complex)


# Generate: Constants

def zeros(shape=(1), dtype=np.float32) -> np.ndarray:
    """
    Generate: Constants: zeros
    """
    assert dtype in DTYPES, f"data type must be in {DTYPES}"
    
    y = np.empty(shape, dtype)

    if dtype == np.float32:
        Y = ffi.cast('float*', y.ctypes.data)
        ret = lib.zeros_s(Y, y.size)
    elif dtype == np.float64 or dtype == float:
        Y = ffi.cast('double*', y.ctypes.data)
        ret = lib.zeros_d(Y, y.size)
    elif dtype == np.complex64:
        Y = ffi.cast('float*', y.ctypes.data)
        ret = lib.zeros_c(Y, y.size)
    else:
        Y = ffi.cast('double*', y.ctypes.data)
        ret = lib.zeros_z(Y, y.size)
    
    assert ret == 0, "error during call to C function"

    return y


def ones(shape=(1), dtype=np.float32) -> np.ndarray:
    """
    Generate: Constants: ones
    """
    assert dtype in DTYPES, f"data type must be in {DTYPES}"
    
    y = np.empty(shape, dtype)

    if dtype == np.float32:
        Y = ffi.cast('float*', y.ctypes.data)
        ret = lib.ones_s(Y, y.size)
    elif dtype == np.float64 or dtype == float:
        Y = ffi.cast('double*', y.ctypes.data)
        ret = lib.ones_d(Y, y.size)
    elif dtype == np.complex64:
        Y = ffi.cast('float*', y.ctypes.data)
        ret = lib.ones_c(Y, y.size)
    else:
        Y = ffi.cast('double*', y.ctypes.data)
        ret = lib.ones_z(Y, y.size)
    
    assert ret == 0, "error during call to C function"

    return y


def twos(shape=(1), dtype=np.float32) -> np.ndarray:
    """
    Generate: Constants: twos
    """
    assert dtype in DTYPES, f"data type must be in {DTYPES}"
    
    y = np.empty(shape, dtype)

    if dtype == np.float32:
        Y = ffi.cast('float*', y.ctypes.data)
        ret = lib.twos_s(Y, y.size)
    elif dtype == np.float64 or dtype == float:
        Y = ffi.cast('double*', y.ctypes.data)
        ret = lib.twos_d(Y, y.size)
    elif dtype == np.complex64:
        Y = ffi.cast('float*', y.ctypes.data)
        ret = lib.twos_c(Y, y.size)
    else:
        Y = ffi.cast('double*', y.ctypes.data)
        ret = lib.twos_z(Y, y.size)
    
    assert ret == 0, "error during call to C function"

    return y


# Main
def main() -> np.ndarray:
    """
    Only used for quick command-line test
    """
    tic = time()
    y = ones((30000, 800), dtype=np.float32)
    print(time()-tic)
    tic = time()
    y = np.ones((30000, 800), dtype=np.float32)
    print(time()-tic)
    return y


if __name__ == "__main__":
    main()
