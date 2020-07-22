General notes on programming:

The functions are organized primarily according to their input/output structure.
This follows much experience with functions and programs to work with inputs and outputs --
the repeated themes in coding that arise follow primarily the numbers and types of inputs and outputs.
Therefore, the major groups of functions are:

Generate: 0 Inputs, 1 Output
Examples: zeros, ones, eye, linspace, logspace, randperm.
"Generator" or "Factory" functions.
These only take option flags, and then generate a matrix or tensor accordingly.

Elementwise1: 1 Input, 1 Output
Examples: sin, cos, rad2deg, deg2rad, exp, log, ceil, floor, erf, abs, square, sqrt.
These operate on one element at a time, and therefore do not require
any information about rows, columns, row-major vs. col-major, etc.
For real-valued data types, it was found faster to use array indexing.
For complex-valued types, it was found faster to use pointer arithmetic.

Elementwise2: 2 Inputs, 1 Output
Examples: plus, minus, times, rdivide, atan2, pow.
These operate on two elements at a time, from two separate inputs.
These all handle array broadcasting, using standard rules.
Therefore, all info about rows, columns, row-major vs. col-major, etc. must be given.
I have a consistent and efficient way to achieve the broadcasting.
Special cases are when input 1 is a scalar (N1==1), input 2 is a scalar (N2==1),
or when inputs 1 and 2 have the same size (N1==N2), in which case no broadcasting is needed.
Otherwise, the broadcasting loops are entered (different for row- vs. col-major),
and use simple pointer arithmetic for clarity and efficiency (has been well-tested).

Complex: 1-2 Inputs, 1 Output
Examples: real, imag, complex, polar, arg.
These operate just like Elementwise1 and Elementwise2 functions,
with the purpose of basic complex-valued operations not covered elsewhere.

Vec2scalar: 1 Input, 1 Output
Examples: norm1, norm2, cnt, sum, mean, std, var, median, mad, range, iqr, prctile.
A.k.a. "Stats" or "Reductions".
These operate on each vector within the input X to produce one scalar per vector.
Option flag -d (--dim) indicates which direction the vectors are oriented within X.
That is, -d0 means col vecs, -d1 means row vecs, etc.
Y has the same size as X, except with only 1 element along dim.
The length of each vector within X is denoted L.
That is, L=R for d=0, L=C for d=1, etc.
To work for tensors up to 4-D, with efficiency and with one unified code loop for many cases,
I designate G global groups (outer loop) and B blocks (inner loop) of vectors.
That is, there are V total vectors within X, each of length L (often V = G*B).
Thus, there are V*L elements in X and V elements in Y.

Vec2vec: 1 Input, 1 Output
Examples: flip, shift, cshift, sort, cumsum, cumprod, prctiles.
These operate on each vector within the input to produce one vector in Y for each vector in X.
Option flag -d (--dim) indicates which direction the vectors are oriented within X.
That is, -d0 means col vecs, -d1 means row vecs, etc.
There are V total vectors within X and within Y.
Y has the same size as X, except that the vector length along dim may be different.
The vector length within X is denoted Lx and the vector length within Y is denoted Ly.
If Lx==Ly, then these are both denoted as L.
For speed, the special case where X is a single 1-D vector is treated separately.
All other cases (tensors filled with 1-D vectors) are treated by the loops over groups and blocks.
Again, a typical program designates G groups and B blocks of vectors (often V = G*B).
For some functions, I found a faster method than iterating through each vector individually (so G and B may mean something different there).

Vec2vec2: 2 Inputs, 1 Output
Examples: dot, cross. [perhaps weighted stats]
These operate on 2 sets of input vectors (X1,X2) with equal vector lengths.
Option flag -d (--dim) indicates which direction the vectors are oriented within X.
That is, -d0 means col vecs, -d1 means row vecs, etc.
The output is a reduction (1 scalar per pair of vectors) or another set of vectors (1 vec per pair of vectors).
These iterate over G groups and B blocks of input vectors, and achieve array broadcasting.

Other: 1-3 Inputs, 1-3 Outputs.
Examples: transpose, ctranspose, join2, split2, row, col, diag, diagmat, toeplitz, repmat.
These are heterogeneous operations that combine, split, rearrange, select, etc.

Linalg: 1-3 Inputs, 1-3 Outputs.
Examples: matnorm, qr, cholesky, eig, svd.
These are heterogeneous operations from linear algebra meant for matrices.
