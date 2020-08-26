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

Input/output is supported for NumPy tensors (https://numpy.org/),  
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
To make an archive library:  
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
All: Generate Construct Matsel Rearrange Split_Join Elementwise1 Elementwise2 Vec2scalar Vecs2scalar Vec2vec Complex Linalg  

Generate: Constants Random Other_Generate  
Constants: zeros ones twos e ln2 ln10 log2e log10e sqrt2 isqrt2 pi ipi pi_2 pi_4 sqrt2 isqrt2 eps realmin realmax inf nan fill  
Random: uniform_int uniform_real bernoulli binomial negative_binomial geometric poisson exponential gamma weibull extreme_value normal lognormal chi_squared cauchy fisher_f student_t   discrete piecewise_constant piecewise_linear  
Other_Generate: eye linspace logspace primes randperm  
 
Construct: diagmat toeplitz tril triu repmat  

Matsel: diag row col  

Rearrange: transpose ctranspose flip sort shift cshift  

Split_Join: split2 split3 join2 join3  

Elementwise1: Operators Trig Exp_Log Round Special  
Operators: plusplus minusminus neg  
Trig: sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh rad2deg deg2rad  
Exp_Log: exp exp2 exp10 log log2 log10  
Round: floor ceil trunc round  
Special: erf erfc tgamma lgamma  
Nonlin: abs square cube sqrt cbrt reciprocal sign deadzone  

Elementwise2: plus minus times rdivide pow hypot atan2 adiff  

Complex: complex polar real imag conj arg norm proj  

Vec2scalar: Sums Prctiles Iprctiles Ranges Norms Moments Other_Means Other_Spreads Other_Stats  
Sums: sum asum cnt  
Prctiles: prctile median max min amax amin  
Iprctiles: iprctile imed imax imin iamax iamin  
Ranges: range iqr idr  
Moments: mean var skewness kurtosis  
Norms: norm0 norm1 norm2 normp  
Other_Means: trimmean winsormean geomean harmean genmean  
Other_Spreads: std geostd trimstd trimvar coeff_var mad  
Other_Stats: prod  

Vecs2scalar: Similarity Dist  
Similarity: dot cov corr cos2 cokurtosis spearman kendall  
Dist: dist0 dist1 dist2 distp dist_cos  

Vec2vec: Center Scale Normalize Reorder Other_Vec2vec  
Center: mean0 med0 geomean1  
Scale: zscore mscore gscore range1 iqr1 idr1  
Normalize: normalize1 normalize2 normalizep  
Reorder: flip shift cshift sort  
Other_Vec2vec: sorti ranks prctiles moments winsorize trim cumsum cumprod softmax betamax  

Linalg: Matmul Transform Sim_Mat Dist_Mat Other_Linalg  
Matmul: kronecker matmul1 matmul1t matmul2 matmul2t matmul3 matmul3t  
Transform: linear affine #projective  
Sim_Mat: dotmat cosmat covmat corrmat spearmat  
Dist_Mat: distmat0 distmat1 distmat2 distmatp distmat_cos  
Other_Linalg: solve chol eig svd  


## Contributing
This is currently for viewing and usage only.  
Feel free to contact the author (erik.edwards4@gmail.com) if any suggestions!


## License
[BSD 3-Clause](https://choosealicense.com/licenses/bsd-3-clause/)

