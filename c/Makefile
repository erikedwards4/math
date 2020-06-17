#@author Erik Edwards
#@date 2017-2020

#math is my own library of C functions for basic math operations in C.

#This is the Makefile used during development of the C functions.
#Since these are not stand-alone programs (no main functions), the only outputs are object files to obj.

SHELL=/bin/bash

CC=clang++

ifeq ($(CC),$(filter $(CC),clang gcc))
	STD=-std=c11
else
	STD=-std=gnu++14
endif

ifeq ($(CC),clang++)
	STD=-std=c++11
	WFLAG=-Weverything -Wno-c++98-compat -Wno-padded -Wno-old-style-cast -Wno-gnu-imaginary-constant -Wno-deprecated
else
	WFLAG=-Wall -Wextra
endif

CFLAGS=$(WFLAG) -O3 $(STD) -march=native -fPIC


all: Generate Construct Matsel Split_Join Nonlin Stats Conv
all: Constants Generate Random Construct Matsel Split_Join Nonlin Stats Conv

#Generate: aka "Factory" functions
#These take 0 inputs (other than parameters) and generate 1 output.
#The Random functions are all done in C++
Generate: Constants Other_Generate #Random

#Constants: 0 inputs, 1 output with a single constant repeated
Constants: zeros ones twos e ln2 ln10 log2e log10e sqrt2 isqrt2 pi ipi pi_2 pi_4 eps realmin realmax inf nan fill
zeros: zeros.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ones: ones.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
twos: twos.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
e: e.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ln2: ln2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ln10: ln10.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
log2e: log2e.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
log10e: log10e.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sqrt2: sqrt2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
isqrt2: isqrt2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
pi: pi.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ipi: ipi.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
pi_2: pi_2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
pi_4: pi_4.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
eps: eps.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
realmin: realmin.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
realmax: realmax.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
inf: inf.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
nan: nan.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
fill: fill.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#Other_Generate: 0 inputs, 1 output
Other_Generate: eye linspace logspace primes
eye: eye.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
linspace: linspace.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
logspace: logspace.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
primes: primes.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Construct: 1-2 vectors input, 1 matrix output constructed from the vectors
Construct: diagmat toeplitz

#Matsel: 1 matrix input, 1 vector output selected from the matrix
Matsel: diag row col

#Split_Join: split 1 matrix into 2-3 matrices, or join 2-3 matrices into 1 matrix
#The join functions accomplish horzcat and vertcat.
Split_Join: split2 split3 join2 join3
split2: split2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
split3: split3.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
join2: join2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
join3: join3.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Elementwise1: 1 input, 1 output element-wise functions
#These include many common static nonlinearities.
Elementwise1: Trig Exp_Log Round Special Nonlin

Trig: sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh
sin: sin.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cos: cos.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
tan: tan.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
asin: asin.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
acos: acos.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
atan: atan.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sinh: sinh.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cosh: cosh.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
tanh: tanh.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
asinh: asinh.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
acosh: acosh.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
atanh: atanh.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Exp_Log: exp exp2 exp10 log log2 log10
exp: exp.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
exp2: exp2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
exp10: exp10.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
log: log.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
log2: log2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
log10: log10.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Round: floor ceil round trunc
floor: floor.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ceil: ceil.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
round: round.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
trunc: trunc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Special: erf tgamma lgamma
erf: erf.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
tgamma: tgamma.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
lgamma: lgamma.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Nonlin: abs square cube sqrt cbrt pow deadzone
abs: abs.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
square: square.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cube: cube.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sqrt: sqrt.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cbrt: cbrt.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
pow: pow.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
deadzone: deadzone.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Stats: some basic statistics calculated row-wise or column-wise
Stats: sum asum cnt min max range iqr prctile median mad mean std var skewness kurtosis coeff_var norm1 norm2 normp
sum: sum.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
asum: asum.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cnt: cnt.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
min: min.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
max: max.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
range: range.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
iqr: iqr.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
prctile: prctile.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
median: median.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
mad: mad.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
mean: mean.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
std: std.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
var: var.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
skewness: skewness.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
kurtosis: kurtosis.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
coeff_var: coeff_var.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
norm1: norm1.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
norm2: norm2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
normp: normp.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#This works with clang++
libmath.a: all
	ar crs ../lib/libmath.a obj/*.o


#make clean
clean:
	find ./obj -type f -name *.o | xargs rm -f

