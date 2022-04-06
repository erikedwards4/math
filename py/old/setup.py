from setuptools import setup, Extension

# Compile *abs.cpp* into a shared library 
setup(
    #...
    ext_modules=[Extension('abs', ['abs.cpp'],),],
)
