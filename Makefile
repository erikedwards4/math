#@author Erik Edwards
#@date 2018-present
#@license BSD 3-clause

#math is my own library of functions for basic math operations in C and C++.
#This is the makefile for the command-line tools in C++.

#Functions that use C-style _Complex must be compiled with g++-6 or g++-7,
#so there is once choice of compiler (CCC=complex CC) for those cases,
#and the usual choice (CC) applies to all the other functions.
#CC can also be used for the 2nd-stage of compilation for the CCC functions.

SHELL=/bin/bash
ss=../util/bin/srci2src
CC=clang++
CCC=g++-7

ifeq ($(CC),clang++)
	STD=-std=c++11
	WFLAG=-Weverything -Wno-c++98-compat -Wno-old-style-cast -Wno-gnu-imaginary-constant -Wno-date-time
else
	STD=-std=gnu++14
	WFLAG=-Wall -Wextra
endif

INCLS=-Ic -I../util
CFLAGS=$(WFLAG) $(STD) -O2 -ffast-math -march=native $(INCLS)
CCFLAGS=-Wall -Wextra -std=gnu++14 -O2 -ffast-math -march=native $(INCLS)
#LIBS=-largtable2 -lopenblas -llapacke -lfftw3f -lfftw3 -lm


All: all
all: Dirs Generate Matmanip Elementwise1 Elementwise2 Vec2scalar Vecs2scalar Vec2vec Complex Linalg Clean

Dirs:
	mkdir -pm 777 bin obj


#Generate: aka "Factory" functions
#These take 0 inputs (other than parameters) and generate 1 output.
#The Random functions are all done in C++
Generate: Constants Other_Gen Random

#Constants: 0 inputs, 1 output with a single constant repeated
Constants: zeros ones twos e ln2 ln10 log2e log10e sqrt2 isqrt2 pi ipi pi_2 pi_4 sqrt2 isqrt2 eps realmin realmax inf nan fill
zeros: srci/zeros.cpp c/zeros.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
ones: srci/ones.cpp c/ones.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
twos: srci/twos.cpp c/twos.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
e: srci/e.cpp c/e.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
ln2: srci/ln2.cpp c/ln2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
ln10: srci/ln10.cpp c/ln10.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
log2e: srci/log2e.cpp c/log2e.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
log10e: srci/log10e.cpp c/log10e.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
sqrt2: srci/sqrt2.cpp c/sqrt2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
isqrt2: srci/isqrt2.cpp c/isqrt2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
pi: srci/pi.cpp c/pi.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
ipi: srci/ipi.cpp c/ipi.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
pi_2: srci/pi_2.cpp c/pi_2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
pi_4: srci/pi_4.cpp c/pi_4.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
eps: srci/eps.cpp c/eps.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 
realmin: srci/realmin.cpp c/realmin.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
realmax: srci/realmax.cpp c/realmax.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
inf: srci/inf.cpp c/inf.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
nan: srci/nan.cpp c/nan.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
fill: srci/fill.cpp c/fill.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

#Other_Gen: 0 inputs, 1 output
Other_Gen: eye linspace logspace primes randperm
eye: srci/eye.cpp c/eye.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
linspace: srci/linspace.cpp c/linspace.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
logspace: srci/logspace.cpp c/logspace.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
primes: srci/primes.cpp c/primes.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
randperm: srci/randperm.cpp c/randperm.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

#Random: 0 inputs, 1 output with random numbers
#The Rand functions do not require PCG random to be installed
#These otherse can use C++ only (the random library), and were originally programmed to do so,
#but are currently coded to use PCG random here (must sudo apt-get install libpcg-cpp-dev)
#If PCG random is not available, then comment out after Rand
Random: Rand Uniform Bernoulli Poisson Normal Sampling

Rand:randi randu randn
randi: srci/randi.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
randu: srci/randu.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
randn: srci/randn.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm

Uniform: uniform_int uniform_real
uniform_int: srci/uniform_int.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
uniform_real: srci/uniform_real.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-conversion; $(CC) obj/$@.o -obin/$@ -largtable2

Bernoulli: bernoulli binomial negative_binomial geometric
bernoulli: srci/bernoulli.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
binomial: srci/binomial.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
negative_binomial: srci/negative_binomial.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
geometric: srci/geometric.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

Poisson: poisson exponential gamma weibull extreme_value
poisson: srci/poisson.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
exponential: srci/exponential.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-conversion; $(CC) obj/$@.o -obin/$@ -largtable2
gamma: srci/gamma.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-conversion; $(CC) obj/$@.o -obin/$@ -largtable2
weibull: srci/weibull.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-conversion; $(CC) obj/$@.o -obin/$@ -largtable2
extreme_value: srci/extreme_value.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-conversion; $(CC) obj/$@.o -obin/$@ -largtable2

Normal: normal lognormal chi_squared cauchy fisher_f student_t
normal: srci/normal.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-conversion; $(CC) obj/$@.o -obin/$@ -largtable2
lognormal: srci/lognormal.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-conversion; $(CC) obj/$@.o -obin/$@ -largtable2
chi_squared: srci/chi_squared.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-conversion; $(CC) obj/$@.o -obin/$@ -largtable2
cauchy: srci/cauchy.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-conversion; $(CC) obj/$@.o -obin/$@ -largtable2
fisher_f: srci/fisher_f.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-conversion; $(CC) obj/$@.o -obin/$@ -largtable2
student_t: srci/student_t.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-conversion; $(CC) obj/$@.o -obin/$@ -largtable2

Sampling: discrete piecewise_constant piecewise_linear
discrete: srci/discrete.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
piecewise_constant: srci/piecewise_constant.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
piecewise_linear: srci/piecewise_linear.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2



#Matmanip: various manipulations of matrix shape and sub-matrix elements
Matmanip: Construct Matsel Rearrange Split_Join

#Construct: 1-2 vectors or matrices input, 1 matrix output constructed from the inputs
Construct: diagmat toeplitz tril triu repmat
diagmat: srci/diagmat.cpp c/diagmat.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
toeplitz: srci/toeplitz.cpp c/toeplitz.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
tril: srci/tril.cpp c/tril.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
triu: srci/triu.cpp c/triu.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
repmat: srci/repmat.cpp c/repmat.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

#Matsel: 1 matrix input, 1 vector output selected from the matrix
Matsel: diag row col
diag: srci/diag.cpp c/diag.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
row: srci/row.cpp c/row.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
col: srci/col.cpp c/col.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

#Rearrange: 1 matrix input, 1 matrix output with elements rearranged
Rearrange: transpose ctranspose rot90
transpose: srci/transpose.cpp c/transpose.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
ctranspose: srci/ctranspose.cpp c/ctranspose.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
rot90: srci/rot90.cpp c/rot90.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

#Split_Join: utilities to split matrix into several, or join several matrices into one.
Split_Join: split2 split3 join2 join3
split2: srci/split2.cpp c/split2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
split3: srci/split3.cpp c/split3.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
join2: srci/join2.cpp c/join2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
join3: srci/join3.cpp c/join3.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2



#Elementwise1 (scalar2scalar): 1 input, 1 output functions that operate element-wise
Elementwise1: Operators Trig Exp_Log Round Special

Operators: plusplus minusminus neg
plusplus: srci/plusplus.cpp c/plusplus.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
minusminus: srci/minusminus.cpp c/minusminus.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
neg: srci/neg.cpp c/neg.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

Trig: sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh rad2deg deg2rad
sin: srci/sin.cpp c/sin.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
cos: srci/cos.cpp c/cos.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
tan: srci/tan.cpp c/tan.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
asin: srci/asin.cpp c/asin.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
acos: srci/acos.cpp c/acos.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
atan: srci/atan.cpp c/atan.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
sinh: srci/sinh.cpp c/sinh.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
cosh: srci/cosh.cpp c/cosh.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
tanh: srci/tanh.cpp c/tanh.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
asinh: srci/asinh.cpp c/asinh.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
acosh: srci/acosh.cpp c/acosh.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
atanh: srci/atanh.cpp c/atanh.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
rad2deg: srci/rad2deg.cpp c/rad2deg.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
deg2rad: srci/deg2rad.cpp c/deg2rad.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm

Exp_Log: exp exp2 exp10 log log2 log10
exp: srci/exp.cpp c/exp.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
exp2: srci/exp2.cpp c/exp2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
exp10: srci/exp10.cpp c/exp10.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
log: srci/log.cpp c/log.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
log2: srci/log2.cpp c/log2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
log10: srci/log10.cpp c/log10.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm

Round: floor ceil trunc round
floor: srci/floor.cpp c/floor.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
ceil: srci/ceil.cpp c/ceil.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
trunc: srci/trunc.cpp c/trunc.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
round: srci/round.cpp c/round.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm

Special: erf erfc tgamma lgamma #sinc bessel
erf: srci/erf.cpp c/erf.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm -lcerf
erfc: srci/erfc.cpp c/erfc.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm -lcerf
tgamma: srci/tgamma.cpp c/tgamma.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
lgamma: srci/lgamma.cpp c/lgamma.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
sinc: srci/sinc.cpp c/sinc.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm

#Nonlin: other 1 input, 1 output static nonlinearities
Nonlin: abs square cube sqrt cbrt reciprocal sign deadzone
abs: srci/abs.cpp c/abs.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
square: srci/square.cpp c/square.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
cube: srci/cube.cpp c/cube.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
sqrt: srci/sqrt.cpp c/sqrt.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
cbrt: srci/cbrt.cpp c/cbrt.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
reciprocal: srci/reciprocal.cpp c/reciprocal.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
sign: srci/sign.cpp c/sign.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
deadzone: srci/deadzone.cpp c/deadzone.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2



#Elementwise2: 2 inputs, 1 output, with array broadcasting for the 2 inputs
Elementwise2: plus minus times rdivide pow hypot atan2 adiff
plus: srci/plus.cpp c/plus.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
minus: srci/minus.cpp c/minus.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
times: srci/times.cpp c/times.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
rdivide: srci/rdivide.cpp c/rdivide.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
pow: srci/pow.cpp c/pow.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
hypot: srci/hypot.cpp c/hypot.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
atan2: srci/atan2.cpp c/atan2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
adiff: srci/adiff.cpp c/adiff.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm



#Vec2scalar: some basic statistics calculated along rows, cols, etc.
#The input consists of 1-D vectors within a tensor, and each vector is reduced to a scalar.
#Thus, these are vec2scalar or "reduction" operations.
Vec2scalar: Sums Prctiles Iprctiles Ranges Norms Moments Other_Means Other_Spreads Other_Stats

Sums: sum asum cnt
sum: srci/sum.cpp c/sum.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
asum: srci/asum.cpp c/asum.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
cnt: srci/cnt.cpp c/cnt.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

Prctiles: prctile median max min amax #amin
prctile: srci/prctile.cpp c/prctile.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke -lm
median: srci/median.cpp c/median.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke
max: srci/max.cpp c/max.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
min: srci/min.cpp c/min.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
amax: srci/amax.cpp c/amax.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
amin: srci/amin.cpp c/amin.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm

Iprctiles: iprctile imed imax imin iamax #iamin
iprctile: srci/iprctile.cpp c/iprctile.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
imed: srci/imed.cpp c/imed.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
imax: srci/imax.cpp c/imax.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
imin: srci/imin.cpp c/imin.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
iamax: srci/iamax.cpp c/iamax.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
iamin: srci/iamin.cpp c/iamin.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm

Ranges: range iqr idr
range: srci/range.cpp c/range.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
iqr: srci/iqr.cpp c/iqr.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke -lm
idr: srci/idr.cpp c/idr.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke -lm

Moments: mean var skewness kurtosis
mean: srci/mean.cpp c/mean.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
var: srci/var.cpp c/var.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
skewness: srci/skewness.cpp c/skewness.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
kurtosis: srci/kurtosis.cpp c/kurtosis.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

Norms: norm0 norm1 norm2 normp
norm0: srci/norm0.cpp c/norm0.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
norm1: srci/norm1.cpp c/norm1.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
norm2: srci/norm2.cpp c/norm2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
normp: srci/normp.cpp c/normp.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm

#Other_Means: trimmed, winsorized, geometric, harmonic, and generalized means (including special case RMS)
Other_Means: trimmean winsormean geomean harmean genmean
trimmean: srci/trimmean.cpp c/trimmean.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke -lm
winsormean: srci/winsormean.cpp c/winsormean.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke -lm
geomean: srci/geomean.cpp c/geomean.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
harmean: srci/harmean.cpp c/harmean.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
genmean: srci/genmean.cpp c/genmean.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
rms: srci/rms.cpp c/rms.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm

Other_Spreads: std geostd trimstd trimvar coeff_var mad
std: srci/std.cpp c/std.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
geostd: srci/geostd.cpp c/geostd.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
trimstd: srci/trimstd.cpp c/trimstd.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke -lm
trimvar: srci/trimvar.cpp c/trimvar.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke -lm
coeff_var: srci/coeff_var.cpp c/coeff_var.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
mad: srci/mad.cpp c/mad.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke -lm

Other_Stats: prod
prod: srci/prod.cpp c/prod.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm


#Vecs2scalar: reduction operations for 2 inputs.
#Each pair of vectors in X1 and X2 is reduced to a scalar.
Vecs2scalar: Similarity Dist #WStats

Similarity: dot cov corr cos2 cokurtosis spearman kendall
dot: srci/dot.cpp c/dot.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
cov: srci/cov.cpp c/cov.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
corr: srci/corr.cpp c/corr.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
cos2: srci/cos2.cpp c/cos2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
cokurtosis: srci/cokurtosis.cpp c/cokurtosis.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
spearman: srci/spearman.cpp c/spearman.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-padded; $(CC) obj/$@.o -obin/$@ -largtable2
kendall: srci/kendall.cpp c/kendall.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

Dist: dist0 dist1 dist2 distp dist_cos
dist0: srci/dist0.cpp c/dist0.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
dist1: srci/dist1.cpp c/dist1.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
dist2: srci/dist2.cpp c/dist2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
distp: srci/distp.cpp c/distp.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
dist_cos: srci/dist_cos.cpp c/dist_cos.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm



#Vec2vec: each vector in tensor X is transformed to a vector in tensor Y
Vec2vec: Center Scale Normalize Reorder Other_Vec2vec

Center: mean0 med0 geomean1
mean0: srci/mean0.cpp c/mean0.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
med0: srci/med0.cpp c/med0.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke
geomean1: srci/geomean1.cpp c/geomean1.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

Scale: zscore mscore gscore range1 iqr1 idr1
zscore: srci/zscore.cpp c/zscore.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
mscore: srci/mscore.cpp c/mscore.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke
gscore: srci/gscore.cpp c/gscore.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
range1: srci/range1.cpp c/range1.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
iqr1: srci/iqr1.cpp c/iqr1.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke
idr1: srci/idr1.cpp c/idr1.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke

Normalize: normalize1 normalize2 normalizep
normalize1: srci/normalize1.cpp c/normalize1.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
normalize2: srci/normalize2.cpp c/normalize2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
normalizep: srci/normalizep.cpp c/normalizep.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm

Reorder: flip shift cshift sort
flip: srci/flip.cpp c/flip.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
shift: srci/shift.cpp c/shift.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
cshift: srci/cshift.cpp c/cshift.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
sort: srci/sort.cpp c/sort.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke -lm

Other_Vec2vec: sorti ranks prctiles moments winsorize trim cumsum cumprod softmax betamax
sorti: srci/sorti.cpp c/sorti.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
ranks: srci/ranks.cpp c/ranks.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-padded; $(CC) obj/$@.o -obin/$@ -largtable2 -lm
prctiles: srci/prctiles.cpp c/prctiles.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke -lm
moments: srci/moments.cpp c/moments.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
winsorize: srci/winsorize.cpp c/winsorize.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke -lm
trim: srci/trim.cpp c/trim.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke -lm
cumsum: srci/cumsum.cpp c/cumsum.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
cumprod: srci/cumprod.cpp c/cumprod.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
prepad: srci/prepad.cpp c/prepad.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
postpad: srci/postpad.cpp c/postpad.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
softmax: srci/softmax.cpp c/softmax.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
betamax: srci/betamax.cpp c/betamax.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm



#Complex: 1-2 inputs, 1 output, for complex-valued operations
Complex: complex polar real imag conj arg proj #norm
complex: srci/complex.cpp c/complex.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
polar: srci/polar.cpp c/polar.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
real: srci/real.cpp c/real.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
imag: srci/imag.cpp c/imag.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
conj: srci/conj.cpp c/conj.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
arg: srci/arg.cpp c/arg.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
proj: srci/proj.cpp c/proj.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
#norm is same as element-wise square



#Linalg: linear algebra routines
#Also see LAPACKE ?large for U*D*U'
Linalg: Matmul Transform Sim_Mat Dist_Mat Other_Linalg

Matmul: kronecker matmul1 matmul1t matmul2 matmul2t matmul3 matmul3t
kronecker: srci/kronecker.cpp c/kronecker.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
matmul1: srci/matmul1.cpp c/matmul1.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
matmul1t: srci/matmul1t.cpp c/matmul1t.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
matmul2: srci/matmul2.cpp c/matmul2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
matmul2t: srci/matmul2t.cpp c/matmul2t.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
matmul3: srci/matmul3.cpp c/matmul3.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
matmul3t: srci/matmul3t.cpp c/matmul3t.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm

Transform: linear affine #projective
linear: srci/linear.cpp c/linear.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
affine: srci/affine.cpp c/affine.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
projective: srci/projective.cpp c/projective.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm

Sim_Mat: dotmat #cosmat covmat corrmat spearmat
dotmat: srci/dotmat.cpp c/dotmat.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
cosmat: srci/cosmat.cpp c/cosmat.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
covmat: srci/covmat.cpp c/covmat.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
corrmat: srci/corrmat.cpp c/corrmat.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
spearmat: srci/spearmat.cpp c/spearmat.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

Dist_Mat: #distmat0 distmat1 distmat2 distmatp distmat_cos
distmat0: srci/distmat0.cpp c/distmat0.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
distmat1: srci/distmat1.cpp c/distmat1.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
distmat2: srci/distmat2.cpp c/distmat2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
distmatp: srci/distmatp.cpp c/distmatp.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
distmat_cos: srci/distmat_cos.cpp c/distmat_cos.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

Other_Linalg: solve chol eig svd #qr matnorm
solve: srci/solve.cpp c/solve.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -llapacke
chol: srci/chol.cpp c/chol.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -llapacke
eig: srci/eig.cpp c/eig.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -llapacke
svd: srci/svd.cpp c/svd.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CCC) -c src/$@.cpp -oobj/$@.o $(CCFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -llapacke
qr: srci/qr.cpp c/qr.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke
matnorm: srci/matnorm.cpp c/matnorm.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas


#make clean
Clean: clean
clean:
	find ./obj -type f -name *.o | xargs rm -f
	rm -f 7 X* Y* x* y*
