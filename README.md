# math

math: basic math functions

================================================

C functions, with associated command-line tools in C++,
that implement many of the usual math function found in NumPy, Octave, etc.

The command-line programs are written in C++ with a consistent style and interface.
The low-level functions themselves are written in C for fastest performance (e.g., openBLAS).

The C functions are meant for the developer; the C++ command-line tools are meant for the end-user.
The interface to each C function is BLAS-like, meaning that one specifies the input and/or output dimensions,
the matrix order as row-major or column-major, and so on.

The C++ command-line programs are written in a consistent style that was developed for command-line tools in general.
All of these command-line tools use argtable2 (http://argtable.sourceforge.net/) for parsing inputs and option flags.
For any of these, use -h (--help) as a flag to get help (description and usage examples).

Input/output is supported for NumPy tensors (https://numpy.org/)
and several C++ tensor formats: Armadillo (http://arma.sourceforge.net/),
ArrayFire (https://arrayfire.com/), and my own minimal format for Eigen (http://eigen.tuxfamily.org/).


## Dependencies
Requires argtable2, openBLAS, LAPACKE.
For Ubuntu, these are available by apt-get:
```console
sudo apt-get install libargtable2-0 libblas3 libopenblas-base liblapack3 liblapacke
```


## Installation
```console
cd /opt/codee
git clone https://github.com/erikedwards4/math
cd /opt/codee/math
make
```

Each C function can also be compiled separately; see c subdirectory Makefile for details.
To make an archive library, do:
```console
cd /opt/codee/math/c
make libmath.a CC=clang
```
This creates /opt/codee/math/lib/libmath.a with all of the C object files.
This could be useful if trying to use the C functions in other applications.
Change clang to clang++ to compile for use with C++ applications.


## Usage
See each resulting command-line tool for help (use -h or --help option).
For example:
```console
/opt/codee/math/bin/log2 --help
```


## List of functions
all: Generate Construct Matsel Rearrange Split_Join Elementwise1 Elementwise2 Complex Stats  

Generate: Constants Other_Gen Random  
    Constants: zeros, ones, twos, e, ln2, ln10, log2e, log10e, sqrt2, isqrt2, pi, ipi, pi_2, pi_4, sqrt2, isqrt2, eps, realmin, realmax, inf, nan, fill  
    Other_Gen: eye, linspace, logspace, primes, randperm  
    Random: uniform_int, uniform_real, bernoulli, binomial, negative_binomial, geometric, poisson, exponential, gamma, weibull, extreme_value, normal, lognormal, chi_squared, cauchy, fisher_f, student_t, discrete, piecewise_constant, piecewise_linear  

Construct: diagmat, toeplitz, tril, triu, repmat  

Matsel: diag, row, col  

Rearrange: transpose, ctranspose, flip, sort, shift, cshift  

Split_Join: split2, split3, join2, join3  

Elementwise1: Operators Trig Exp_Log Round Special  
    Operators: plusplus, minusminus
    Trig: sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh, rad2deg, deg2rad  
    Exp_Log: exp, exp2, exp10, log, log2, log10  
    Round: floor, ceil, trunc, round  
    Special: erf, erfc, tgamma, lgamma  
    Nonlin: abs, square, cube, sqrt, cbrt, pow, deadzone  

Elementwise2: plus, minus, times, rdivide, pow, hypot, atan2  

Complex: complex, polar, real, imag, conj, arg, norm, proj

Stats: Sums Prctiles Ranges Moments Norms Other_Stats
    Sums: sum, asum, cnt  
    Prctiles: min, max, amin, amax, median, prctile, prctiles
    Ranges: range, iqr, interdecile_range  
    Moments: mean, var, skewness, kurtosis  
    Norms: norm1, norm2, normp  
    Other_Stats: std, coeff_var, mad  

## Contributing
This is currently for viewing and usage only.  
Feel free to contact the author (erik.edwards4@gmail.com) if any suggestions!


## License
[BSD 3-Clause](https://choosealicense.com/licenses/bsd-3-clause/)

