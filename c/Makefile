#@author Erik Edwards
#@date 2017-present
#@license BSD 3-clause

#My own library of C functions for general math operations.

#This is the Makefile used during development of the C functions.
#Since these are not stand-alone programs (no main functions), the only outputs are object files to obj.
#It is not necessary to run this when making the C++ command-line tools.

#Needed external libraries: cblas, lapacke,
#Other external libraries: libcerf-dev (complex error function, but I commented these out for now)
#To Do: weighted statistics


SHELL=/bin/bash
CC=clang

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

CFLAGS=$(WFLAG) $(STD) -O2 -ffast-math -march=native -mfpmath=sse -fPIC


All: all
all: Dirs Generate Matmanip Elementwise1 Elementwise2 Vec2scalar Vec2vec Vecs2scalar Linalg Lib Clean

Dirs:
	mkdir -pm 777 obj


#Generate: aka "Factory" functions
#These take 0 inputs (other than parameters) and generate 1 output.
#The Random functions here are stand-alone (based on PCG random, but no dependencies).
Generate: Constants Rand Other_Gen

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

#Rand: 0 inputs, 1 output
#These are the fastest and go-to versions.
#They use modified code from PCG random, but are stand-alone (no install of PCG required)
Rand: randi randu randn
randi: randi.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
randu: randu.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
randn: randn.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#Other_Gen: 0 inputs, 1 output
Other_Gen: eye linspace logspace primes randperm
eye: eye.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
linspace: linspace.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
logspace: logspace.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
primes: primes.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
randperm: randperm.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)



#Matmanip:
Matmanip: Construct Matsel Rearrange Split_Join

#Construct: 1-2 vectors or matrices input, 1 matrix output constructed from the input(s)
Construct: diagmat toeplitz tril triu repmat
diagmat: diagmat.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
toeplitz: toeplitz.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
tril: tril.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
triu: triu.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
repmat: repmat.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#Matsel: 1 matrix input, 1 vector output selected from the matrix
Matsel: diag row col
diag: diag.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
row: row.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
col: col.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#Rearrange: 1 matrix input, 1 matrix output with elements rearranged
Rearrange: transpose ctranspose rot90
transpose: transpose.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ctranspose: ctranspose.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
rot90: rot90.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#Split_Join: split 1 matrix into 2-3 matrices, or join 2-3 matrices into 1 matrix
#The join functions accomplish horzcat and vertcat.
Split_Join: split2 split3 join2 join3
split2: split2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
split3: split3.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
join2: join2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
join3: join3.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)



#Elementwise1 (scalar2scalar): 1 input, 1 output element-wise functions
#These include many common static nonlinearities.
Elementwise1: Operators Trig Exp_Log Round Special Complex Nonlin

Operators: plusplus minusminus neg
plusplus: plusplus.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
minusminus: minusminus.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
neg: neg.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Trig: sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh rad2deg deg2rad
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

Special: erf erfc tgamma lgamma #sinc bessel
erf: erf.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
erfc: erfc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
tgamma: tgamma.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
lgamma: lgamma.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
# sinc: sinc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Complex: real imag arg conj proj #norm
real: real.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
imag: imag.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
arg: arg.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
conj: conj.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
proj: proj.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
#norm is same as element-wise square

Nonlin: abs square cube sqrt cbrt nextpow2 reciprocal sign deadzone
abs: abs.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
square: square.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cube: cube.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sqrt: sqrt.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cbrt: cbrt.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
nextpow2: nextpow2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
reciprocal: reciprocal.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sign: sign.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
deadzone: deadzone.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Elementwise2: 2 inputs, 1 output, with array broadcasting for the 2 inputs
Elementwise2: Arithmetic Trig2 Complex2

Arithmetic: plus minus times rdivide adiff pow
plus: plus.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
minus: minus.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
times: times.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
rdivide: rdivide.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
adiff: adiff.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
pow: pow.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Trig2: hypot atan2
hypot: hypot.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
atan2: atan2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Complex2: complex polar
complex: complex.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
polar: polar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)



#Vec2scalar (Stats): 1 input, 1 output
#Some basic statistics calculated along rows, cols, etc.
#The input consists of 1-D vectors within a tensor, and each vector is reduced to a scalar.
#Thus, these are vec2scalar or "reduction" operations.
Vec2scalar: Sums Prctiles Iprctiles Ranges Norms Moments Other_Means Other_Spreads Other_Stats

Sums: sum asum cnt
sum: sum.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
asum: asum.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cnt: cnt.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Prctiles: prctile median max min amax #amin
prctile: prctile.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
median: median.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
max: max.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
min: min.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
amax: amax.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
amin: amin.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Iprctiles: iprctile imed imax imin iamax #iamin
iprctile: iprctile.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
imed: imed.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
imax: imax.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
imin: imin.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
iamax: iamax.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
iamin: iamin.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Ranges: range iqr idr
range: range.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
iqr: iqr.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
idr: idr.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Norms: norm0 norm1 norm2 normp
norm0: norm0.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
norm1: norm1.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
norm2: norm2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
normp: normp.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Moments: mean var skewness kurtosis
mean: mean.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
var: var.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
skewness: skewness.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
kurtosis: kurtosis.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#Other_Means: trimmed_mean geometric_mean harmonic_mean generalized_mean root_mean_square
Other_Means: trimmean winsormean geomean harmean genmean rms
trimmean: trimmean.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
winsormean: winsormean.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
geomean: geomean.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
harmean: harmean.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
genmean: genmean.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
rms: rms.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Other_Spreads: std geostd trimstd trimvar coeff_var mad
std: std.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
geostd: geostd.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
trimstd: trimstd.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
trimvar: trimvar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
coeff_var: coeff_var.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
mad: mad.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Other_Stats: prod qselect
prod: prod.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
qselect: qselect.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)



#Vec2vec: 1 input, 1 output
#Operates on each 1-D vector in X to produce a new 1-D vector in Y.
#Y has the same size as X, except perhaps a change in the vector length.
Vec2vec: Center Scale Normalize Sort_Help Sorti_Help Sortif_Help Reorder Isort Other_Vec2vec

Center: mean0 med0 geomean1
mean0: mean0.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
med0: med0.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
geomean1: geomean1.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Scale: zscore mscore gscore range1 iqr1 idr1
zscore: zscore.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
mscore: mscore.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
gscore: gscore.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
range1: range1.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
iqr1: iqr1.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
idr1: idr1.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Normalize: normalize1 normalize2 normalizep
normalize1: normalize1.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
normalize2: normalize2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
normalizep: normalizep.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Sort_Help: cmp extremum lomuto_partition hoare_partition insertion_sort quicksort partial_sort kselect
cmp: cmp.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
extremum: extremum.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
lomuto_partition: lomuto_partition.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
hoare_partition: hoare_partition.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
insertion_sort: insertion_sort.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
quicksort: quicksort.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
partial_sort: partial_sort.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
kselect: kselect.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Sorti_Help: cmpi extremumi lomuto_partitioni hoare_partitioni insertion_sorti quicksorti partial_sorti kselecti
cmpi: cmpi.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
extremumi: extremumi.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
lomuto_partitioni: lomuto_partitioni.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
hoare_partitioni: hoare_partitioni.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
insertion_sorti: insertion_sorti.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
quicksorti: quicksorti.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
partial_sorti: partial_sorti.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
kselecti: kselecti.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Sortif_Help: cmpif extremumif lomuto_partitionif hoare_partitionif insertion_sortif quicksortif partial_sortif kselectif
cmpif: cmpif.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
extremumif: extremumif.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
lomuto_partitionif: lomuto_partitionif.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
hoare_partitionif: hoare_partitionif.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
insertion_sortif: insertion_sortif.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
quicksortif: quicksortif.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
partial_sortif: partial_sortif.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
kselectif: kselectif.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Reorder: flip shift cshift insert_sort sort lapacke_sort partsort
flip: flip.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
shift: shift.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cshift: cshift.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
insert_sort: insert_sort.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sort: sort.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
lapacke_sort: lapacke_sort.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
partsort: partsort.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Isort: insert_sorti qsorti sorti insert_ranks qranks ranks
insert_sorti: insert_sorti.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
qsorti: qsorti.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sorti: sorti.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
insert_ranks: insert_ranks.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
qranks: qranks.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ranks: ranks.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Other_Vec2vec: prctiles moments winsorize trim cumsum cumprod prepad postpad softmax betamax
prctiles: prctiles.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
moments: moments.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
winsorize: winsorize.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
trim: trim.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cumsum: cumsum.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cumprod: cumprod.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
prepad: prepad.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
postpad: postpad.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
softmax: softmax.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
betamax: betamax.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)



#Vecs2scalar: 2 inputs, 1 output
#Operates on each pair of 1-D vectors in X1, X2 to produce a scalar in Y.
#This also allows vector broadcasting (like Elementwise2, but limited to single vectors).
Vecs2scalar: Similarity Dist WStats

#Similarity or related measures for 2 vectors
Similarity: dot cov corr corr_opus cos2 cokurtosis spearman
dot: dot.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cov: cov.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
corr: corr.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
corr_opus: corr_opus.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cos2: cos2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cokurtosis: cokurtosis.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
spearman: spearman.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
kendall: kendall.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#Distance metrics or related for 2 vectors
Dist: dist0 dist1 dist2 distp dist_cos
dist0: dist0.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
dist1: dist1.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
dist2: dist2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
distp: distp.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
dist_cos: dist_cos.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#Weighted statistics
WStats: #WSums WNorms WPrctiles WRanges WMoments WOther

WSums: wsum wasum wcnt
wsum: wsum.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
wasum: wasum.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
wcnt: wcnt.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

WNorms: wnorm0 wnorm1 wnorm2 wnormp
wnorm0: wnorm0.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
wnorm1: wnorm1.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
wnorm2: wnorm2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
wnormp: wnormp.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

WPrctiles: wprctile wmedian
wprctile: wprctile.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
wmedian: wmedian.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

WRanges: wiqr widr
wiqr: wiqr.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
widr: widr.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

WMoments: wmean wvar wskewness wkurtosis
wmean: wmean.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
wvar: wvar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
wskewness: wskewness.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
wkurtosis: wkurtosis.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

WOther: wstd wcoeff_var wmad
wstd: wstd.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
wcoeff_var: wcoeff_var.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
wmad: wmad.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)



#Linalg: linear algebra
Linalg: Matmul Transform Sim_Mat Dist_Mat Other_Linalg

Matmul: matmul1 matmul1t matmul2 matmul2t matmul3 matmul3t kronecker
matmul1: matmul1.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
matmul1t: matmul1t.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
matmul2: matmul2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
matmul2t: matmul2t.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
matmul3: matmul3.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
matmul3t: matmul3t.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
kronecker: kronecker.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Transform: linear affine #projective
linear: linear.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
affine: affine.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
projective: projective.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Sim_Mat: dotmat cosmat #covmat corrmat
dotmat: dotmat.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cosmat: cosmat.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
covmat: covmat.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
corrmat: corrmat.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Dist_Mat: #distmat0 distmat1 distmat2 distmatp distmat_cos
distmat0: distmat0.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
distmat1: distmat1.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
distmat2: distmat2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
distmatp: distmatp.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
distmat_cos: distmat_cos.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

Other_Linalg: solve chol eig svd #qr matnorm
solve: solve.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
chol: chol.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
eig: eig.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
svd: svd.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
qr: qr.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
matnorm: matnorm.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#If want to make a library archive for C++ use, this works with clang++
Lib: lib
lib:
	ar crs ../lib/libmath.a obj/*.o
	gcc -fPIC -shared -o ../lib/libmath.so obj/*.o -lopenblas -llapacke -lfftw3f -lfftw3 -lm


#make clean
#find ./obj -type f -name *.o | xargs rm -f
Clean: clean
clean:
	rm -f *.s tmp.*
