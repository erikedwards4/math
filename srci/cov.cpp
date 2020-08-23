//Includes
#include "cov.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 2u, O = 1u;
size_t dim;
char b;

//Description
string descr;
descr += "Vecs2scalar function for 2 inputs with broadcasting.\n";
descr += "Gets covariance for each pair of vectors in X1, X2 along dim.\n";
descr += "Thus, output Y has length 1 along dim.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
descr += "Use -d0 to get cov along cols.\n";
descr += "Use -d1 to get cov along rows.\n";
descr += "Use -d2 to get cov along slices.\n";
descr += "Use -d3 to get cov along hyperslices.\n";
descr += "\n";
descr += "X1 and X2 must have the same size and data type, or\n";
descr += "either X1 or X2 can be a single vector of appropriate length.\n";
descr += "In that case, the single vector will be broadcast to the other input.\n";
descr += "\n";
descr += "Include -b (--biased) to use the biased denominator [default is unbiased].\n";
descr += "The biased denominator is N, and the unbiased denominator is N-1.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ cov X1 X2 -o Y \n";
descr += "$ cov -d1 X1 X2 > Y \n";
descr += "$ cat X1 | cov -d2 - X2 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
struct arg_lit    *a_b = arg_litn("b","biased",0,1,"use biased (N) denominator [default=N-1]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1u,2u,3}" << endl; return 1; }

//Get b
b = (a_b->count>0);

//Checks
if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }
if (!bcast_compat(i1,i2)) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have same size or broadcast-compatible sizes" << endl; return 1; }
if (!major_compat(i1,i2)) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same row/col major format unless vectors" << endl; return 1; }
size_t L1 = (dim==0u) ? i1.R : (dim==1u) ? i1.C : (dim==2u) ? i1.S : i1.H;
size_t L2 = (dim==0u) ? i2.R : (dim==1u) ? i2.C : (dim==2u) ? i2.S : i2.H;
if (L1!=L2) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have same length along dim" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0u) ? 1 : (i1.R>i2.R) ? i1.R : i2.R;
o1.C = (dim==1u) ? 1 : (i1.C>i2.C) ? i1.C : i2.C;
o1.S = (dim==2u) ? 1 : (i1.S>i2.S) ? i1.S : i2.S;
o1.H = (dim==3u) ? 1 : (i1.H>i2.H) ? i1.H : i2.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X1, *X2, *Y;
    try { X1 = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
    try { X2 = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
    if (codee::cov_s(Y,X1,X2,i1.R,i1.C,i1.S,i1.H,i2.R,i2.C,i2.S,i2.H,o1.iscolmajor(),dim,b))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2; delete[] Y;
}
else if (i1.T==101u)
{
    float *X1, *X2, *Y;
    try { X1 = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
    try { X2 = new float[2u*i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
    if (codee::cov_c(Y,X1,X2,i1.R,i1.C,i1.S,i1.H,i2.R,i2.C,i2.S,i2.H,o1.iscolmajor(),dim,b))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2; delete[] Y;
}

//Finish
