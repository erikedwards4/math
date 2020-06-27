#@author Erik Edwards
#@date 2019-2020

#math is my own library of functions for basic math operations in C and C++.
#This is the makefile for the command-line tools in C++.

#Project is organized as:

SHELL=/bin/bash

ss=../util/bin/srci2src

CC=clang++


ifeq ($(CC),clang++)
	STD=-std=c++11
	WFLAG=-Weverything -Wno-c++98-compat -Wno-padded -Wno-old-style-cast -Wno-gnu-imaginary-constant
else
	STD=-std=gnu++14
	WFLAG=-Wall -Wextra
endif

CFLAGS=$(WFLAG) -O3 $(STD) -march=native -Ic
#LIBS=-largtable2 -lopenblas -llapacke -llapack -lfftw3f -lfftw3 -lm


all: Generate Construct Matsel Rearrange Split_Join Elementwise1 Elementwise2 Complex Stats Linalg
	rm -f 7 obj/*.o


#Generate: aka "Factory" functions
#These take 0 inputs (other than parameters) and generate 1 output.
#The Random functions are all done in C++
Generate: Constants Other_Generate Random

#Constants: 0 inputs, 1 output with a single constant repeated
Constants: zeros ones twos e ln2 ln10 log2e log10e sqrt2 isqrt2 pi ipi pi_2 pi_4 sqrt2 isqrt2 eps realmin realmax inf nan fill
zeros: srci/zeros.cpp c/zeros.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
ones: srci/ones.cpp c/ones.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
twos: srci/twos.cpp c/twos.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
e: srci/e.cpp c/e.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
ln2: srci/ln2.cpp c/ln2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
ln10: srci/ln10.cpp c/ln10.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
log2e: srci/log2e.cpp c/log2e.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
log10e: srci/log10e.cpp c/log10e.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
sqrt2: srci/sqrt2.cpp c/sqrt2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
isqrt2: srci/isqrt2.cpp c/isqrt2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
pi: srci/pi.cpp c/pi.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
ipi: srci/ipi.cpp c/ipi.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
pi_2: srci/pi_2.cpp c/pi_2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
pi_4: srci/pi_4.cpp c/pi_4.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
eps: srci/eps.cpp c/eps.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
realmin: srci/realmin.cpp c/realmin.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
realmax: srci/realmax.cpp c/realmax.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
inf: srci/inf.cpp c/inf.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
nan: srci/nan.cpp c/nan.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
fill: srci/fill.cpp c/fill.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas

#Other_Generate: 0 inputs, 1 output
Other_Generate: eye linspace logspace primes randperm
eye: srci/eye.cpp c/eye.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
linspace: srci/linspace.cpp c/linspace.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
logspace: srci/logspace.cpp c/logspace.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
primes: srci/primes.cpp c/primes.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
randperm: srci/randperm.cpp c/randperm.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

#Random: 0 inputs, 1 output with random numbers
#These use C++ only (the random library for now; PCG random in future)
Random: Uniform Bernoulli Poisson Normal Sampling Randperm

Uniform: uniform_int uniform_real
uniform_int: srci/uniform_int.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
uniform_real: srci/uniform_real.cpp; $(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-conversion; $(CC) obj/$@.o -obin/$@ -largtable2

Bernouli: bernoulli binomial negative_binomial geometric
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



#Construct: 1-2 vectors or matrices input, 1 matrix output constructed from the vectors
Construct: diagmat toeplitz tril triu repmat #prepad postpad
diagmat: srci/diagmat.cpp c/diagmat.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
toeplitz: srci/toeplitz.cpp c/toeplitz.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
tril: srci/tril.cpp c/tril.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
triu: srci/triu.cpp c/triu.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
prepad: srci/prepad.cpp c/prepad.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
postpad: srci/postpad.cpp c/postpad.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
repmat: srci/repmat.cpp c/repmat.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas


#Matsel: 1 matrix input, 1 vector output selected from the matrix
Matsel: diag row col
diag: srci/diag.cpp c/diag.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
row: srci/row.cpp c/row.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
col: srci/col.cpp c/col.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas


#Rearrange: 1 matrix input, 1 matrix output with elements rearranged
Rearrange: transpose ctranspose flip sort shift cshift #rot90
transpose: srci/transpose.cpp c/transpose.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
ctranspose: srci/ctranspose.cpp c/ctranspose.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
flip: srci/flip.cpp c/flip.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
sort: srci/sort.cpp c/sort.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -llapacke
shift: srci/shift.cpp c/shift.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
cshift: srci/cshift.cpp c/cshift.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
rot90: srci/rot90.cpp c/rot90.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2


#Split_Join: utilities to split matrix into several, or join several matrices into one.
Split_Join: split2 split3 join2 join3
split2: srci/split2.cpp c/split2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
split3: srci/split3.cpp c/split3.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
join2: srci/join2.cpp c/join2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
join3: srci/join3.cpp c/join3.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas



#Elementwise1: 1 input, 1 output functions that operate element-wise
Elementwise1: Operators Trig Exp_Log Round Special

Operators: plusplus minusminus
plusplus: srci/plusplus.cpp c/plusplus.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
minusminus: srci/minusminus.cpp c/minusminus.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

Trig: sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh rad2deg deg2rad
sin: srci/sin.cpp c/sin.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
cos: srci/cos.cpp c/cos.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
tan: srci/tan.cpp c/tan.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
asin: srci/asin.cpp c/asin.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
acos: srci/acos.cpp c/acos.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
atan: srci/atan.cpp c/atan.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
sinh: srci/sinh.cpp c/sinh.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
cosh: srci/cosh.cpp c/cosh.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
tanh: srci/tanh.cpp c/tanh.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
asinh: srci/asinh.cpp c/asinh.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
acosh: srci/acosh.cpp c/acosh.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
atanh: srci/atanh.cpp c/atanh.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
rad2deg: srci/rad2deg.cpp c/rad2deg.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
deg2rad: srci/deg2rad.cpp c/deg2rad.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm

Exp_Log: exp exp2 exp10 log log2 log10
exp: srci/exp.cpp c/exp.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
exp2: srci/exp2.cpp c/exp2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
exp10: srci/exp10.cpp c/exp10.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
log: srci/log.cpp c/log.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
log2: srci/log2.cpp c/log2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
log10: srci/log10.cpp c/log10.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm

Round: floor ceil trunc round
floor: srci/floor.cpp c/floor.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
ceil: srci/ceil.cpp c/ceil.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
trunc: srci/trunc.cpp c/trunc.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
round: srci/round.cpp c/round.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm

Special: erf erfc tgamma lgamma
erf: srci/erf.cpp c/erf.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
erfc: srci/erfc.cpp c/erfc.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
tgamma: srci/tgamma.cpp c/tgamma.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
lgamma: srci/lgamma.cpp c/lgamma.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm

#Nonlin: other 1 input, 1 output static nonlinearities
Nonlin: abs square cube sqrt cbrt pow deadzone
abs: srci/abs.cpp c/abs.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
square: srci/square.cpp c/square.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
cube: srci/cube.cpp c/cube.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
sqrt: srci/sqrt.cpp c/sqrt.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
cbrt: srci/cbrt.cpp c/cbrt.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
deadzone: srci/deadzone.cpp c/deadzone.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2



#Elementwise2: 2 inputs, 1 output, with array broadcasting for the 2 inputs
Elementwise2: plus minus times rdivide pow hypot atan2
plus: srci/plus.cpp c/plus.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
minus: srci/minus.cpp c/minus.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
times: srci/times.cpp c/times.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
rdivide: srci/rdivide.cpp c/rdivide.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
pow: srci/pow.cpp c/pow.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
hypot: srci/hypot.cpp c/hypot.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
atan2: srci/atan2.cpp c/atan2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm


#Complex: 1-2 inputs, 1 output, for complex operations
Complex: complex polar real imag conj arg norm proj
complex: srci/complex.cpp c/complex.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
polar: srci/polar.cpp c/polar.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
real: srci/real.cpp c/real.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
imag: srci/imag.cpp c/imag.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
conj: srci/conj.cpp c/conj.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
arg: srci/arg.cpp c/arg.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
norm: srci/norm.cpp c/norm.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
proj: srci/proj.cpp c/proj.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2



#Stats: some basic statistics calculated row-wise or column-wise
Stats: Sums Prctiles Ranges Moments Norms Other_Stats

Sums: sum asum cnt
sum: srci/sum.cpp c/sum.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
asum: srci/asum.cpp c/asum.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
cnt: srci/cnt.cpp c/cnt.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas

Prctiles: min max amin amax median prctile prctiles
min: srci/min.cpp c/min.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
max: srci/max.cpp c/max.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
amin: srci/amin.cpp c/amin.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
amax: srci/amax.cpp c/amax.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
median: srci/median.cpp c/median.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
prctile: srci/prctile.cpp c/prctile.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas

Ranges: range iqr interdecile_range
range: srci/range.cpp c/range.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
iqr: srci/iqr.cpp c/iqr.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
interdecile_range: srci/interdecile_range.cpp c/interdecile_range.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm

Moments: mean var skewness kurtosis
mean: srci/mean.cpp c/mean.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
var: srci/var.cpp c/var.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
skewness: srci/skewness.cpp c/skewness.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
kurtosis: srci/kurtosis.cpp c/kurtosis.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm

Norms: norm1 norm2 normp
norm1: srci/norm1.cpp c/norm1.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
norm2: srci/norm2.cpp c/norm2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
normp: srci/normp.cpp c/normp.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm

Other_Stats: std coeff_var mad
std: srci/std.cpp c/std.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
coeff_var: srci/coeff_var.cpp c/coeff_var.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
mad: srci/mad.cpp c/mad.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm



#Linalg: linear algebra routines
#Also see LAPACKE ?large for U*D*U'
Linalg: matmul matmul3 kronecker matnorm lq qr lu

clean:
	find ./obj -type f -name *.o | xargs rm -f
