//Includes
#include "eig.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 1u, O = 2u;
size_t K;

//Description
string descr;
descr += "Linear algebra function.\n";
descr += "Eigendecomposition of square, Hermitian-symmetric matrix X.\n";
descr += "Only the lower-triangular part of X is accessed.\n";
descr += "This returns U and V such that X = U * diagmat(V) * U'.\n";
descr += "\n";
descr += "The outputs U and V are the eigenvectors and eigenvalues, respectively.\n";
descr += "List each output filename following a -o opt, in the order U then V.\n";
descr += "\n";
descr += "Use -k (--K) to give the number of eigenvals and eigenvecs to keep [default=all].\n";
descr += "The K largest eigenvalues and corresponding eigenvectors are returned.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ eig X -o U -o V \n";
descr += "$ eig -k7 X -o U > V \n";
descr += "$ cat X | eig -o - -o V > U \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_k = arg_intn("k","K","<uint>",0,1,"num eigencomponents to keep [default=all]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output files (U,V)");

//Get options

//Get K
if (a_k->count==0) { K = i1.R; }
else if (a_k->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "K must be positive" << endl; return 1; }
else { K = size_t(a_k->ival[0]); }
if (K>i1.R) { cerr << progstr+": " << __LINE__ << errstr << "K must be <= R (nrows in X)" << endl; return 1; }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be a matrix" << endl; return 1; }
if (!i1.issquare()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be a square matrix" << endl; return 1; }

//Set output header info
o1.F = o2.F = i1.F;
o1.T = i1.T;
o2.T = i1.isreal() ? i1.T : i1.T-100u;
o1.R = i1.R; o2.R = K;
o1.C = K; o2.C = 1u;
o1.S = o2.S = i1.S;
o1.H = o2.H = i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X, *V;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { V = new float[i1.R]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (V)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::eig_inplace_s(X,V,i1.R,i1.iscolmajor(),K))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 1 (U)" << endl; return 1; }
    }
    if (wo2)
    {
        try { ofs2.write(reinterpret_cast<char*>(V),o2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 2 (V)" << endl; return 1; }
    }
    delete[] X; delete[] V;
}
else if (i1.T==101u)
{
    float *X, *V;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { V = new float[2u*i1.R]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (V)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::eig_inplace_c(X,V,i1.R,i1.iscolmajor(),K))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 1 (U)" << endl; return 1; }
    }
    if (wo2)
    {
        try { ofs2.write(reinterpret_cast<char*>(V),o2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 2 (V)" << endl; return 1; }
    }
    delete[] X; delete[] V;
}

//Finish
