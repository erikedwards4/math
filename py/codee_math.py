#!/usr/bin/env python3
"""
Makes ctypes interface to the C functions in math.so.
"""

# import ctypes
from cffi import FFI
ffibuilder = FFI()

# cdef() expects a single string declaring the C types, functions and
# globals needed to use the shared object. It must be in valid C syntax.
ffibuilder.cdef(
    """
    int abs_s (float *Y, const float *X, const size_t N);
    """)

# set_source() gives the name of the python extension module to
# produce, and some C source code as a string.  This C code needs
# to make the declarated functions, types and globals available,
# so it is often just the "#include".
ffibuilder.set_source("_codee_math_cffi",
"""
     #include "codee_math.h"   // the C header of the library
""",
    library_dirs = [],
    libraries = ['m'],   # library name, for the linker
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
