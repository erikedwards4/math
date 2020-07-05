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
	WFLAG=-Weverything -Wno-c++98-compat -Wno-padded -Wno-gnu-imaginary-constant -Wno-deprecated -Wno-old-style-cast
else
	WFLAG=-Wall -Wextra
endif

CFLAGS=$(WFLAG) -O3 $(STD) -march=native -fPIC


all: Generate Construct Matsel Rearrange Split_Join Elementwise1 Elementwise2 Complex Stats
	rm -f 7


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
Other_Generate: eye linspace logspace primes randperm
eye: eye.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
linspace: linspace.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
logspace: logspace.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
primes: primes.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
randperm: randperm.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)



#Construct: 1-2 vectors or matrices input, 1 matrix output constructed from the input(s)
Construct: diagmat toeplitz tril triu repmat #prepad postpad
diagmat: diagmat.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
toeplitz: toeplitz.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
tril: tril.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
triu: triu.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
prepad: prepad.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
postpad: postpad.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
repmat: repmat.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Matsel: 1 matrix input, 1 vector output selected from the matrix
Matsel: diag row col
diag: diag.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
row: row.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
col: col.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Split_Join: split 1 matrix into 2-3 matrices, or join 2-3 matrices into 1 matrix
#The join functions accomplish horzcat and vertcat.
Split_Join: split2 split3 join2 join3
split2: split2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
split3: split3.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
join2: join2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
join3: join3.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Rearrange: 1 matrix input, 1 matrix output with element rearranged
Rearrange: transpose ctranspose flip shift cshift sort #rot90
transpose: transpose.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ctranspose: ctranspose.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
flip: flip.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
shift: shift.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cshift: cshift.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sort: sort.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
rot90: rot90.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)



#Elementwise1: 1 input, 1 output element-wise functions
#These include many common static nonlinearities.
Elementwise1: Operators Trig Exp_Log Round Special Nonlin

Operators: plusplus minusminus
plusplus: plusplus.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
minusminus: minusminus.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Trig: sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh #rad2deg deg2rad
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
rad2deg: rad2deg.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
deg2rad: deg2rad.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

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

Special: erf erfc tgamma lgamma
erf: erf.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
erfc: erfc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
tgamma: tgamma.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
lgamma: lgamma.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Nonlin: abs square cube sqrt cbrt deadzone
abs: abs.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
square: square.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cube: cube.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sqrt: sqrt.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cbrt: cbrt.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
deadzone: deadzone.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Elementwise2: 2 inputs, 1 output, with array broadcasting for the 2 inputs
Elementwise2: plus minus times rdivide pow hypot atan2
plus: plus.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
minus: minus.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
times: times.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
rdivide: rdivide.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
pow: pow.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
hypot: hypot.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
atan2: atan2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)



#Complex: 1-2 inputs, 1 output, for complex operations
Complex: complex polar real imag conj arg norm proj
complex: complex.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
polar: polar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
real: real.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
imag: imag.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
conj: conj.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
arg: arg.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
norm: norm.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
proj: proj.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)



#Stats: some basic statistics calculated row-wise or column-wise
Stats: Sums Prctiles Ranges Moments Norms Other_Stats

Sums: sum asum cnt
sum: sum.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
asum: asum.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cnt: cnt.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Prctiles: min max amin amax median prctile prctiles
min: min.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
max: max.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
amin: amin.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
amax: amax.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
median: median.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
prctile: prctile.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Ranges: range iqr interdecile_range
range: range.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
iqr: iqr.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
interdecile_range: interdecile_range.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Moments: mean var skewness kurtosis
mean: mean.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
var: var.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
skewness: skewness.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
kurtosis: kurtosis.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Norms: norm1 norm2 normp
norm1: norm1.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
norm2: norm2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
normp: normp.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Other_Stats: std coeff_var mad
std: std.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
coeff_var: coeff_var.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
mad: mad.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)



#This works with clang++
libmath.a: all
	ar crs ../lib/libmath.a obj/*.o


#make clean
clean:
	find ./obj -type f -name *.o | xargs rm -f

