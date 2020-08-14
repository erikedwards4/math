//Includes
#include "svd.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 1, O = 3;
size_t Kmax, K;

//Description
string descr;
descr += "Linear algebra function.\n";
descr += "Singular decomposition of matrix X.\n";
descr += "\n";
descr += "The outputs U and V are the left and right singular vectors, respectively.\n";
descr += "The output S is a vector of singular values, such that: X = U*diagmat(S)*V'.\n";
descr += "\n";
descr += "Use -k (--K) to give the number of singular vals and vecs to keep [default=all].\n";
descr += "The K largest singular values and corresponding vectors are returned.\n";
descr += "The default for K is K_all = min(R,C). \n";
descr += "\n";
descr += "List each output filename following a -o opt, in the order U, S, V.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ svd X -o U -o S -o V \n";
descr += "$ svd X -o U -o - -o V > S \n";
descr += "$ cat X | svd -o U -o S -o V \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_k = arg_intn("k","K","<uint>",0,1,"num svdencomponents to keep [default=all]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output files (U,S,V)");

//Get options

//Get K
Kmax = (i1.R<i1.C) ? i1.R : i1.C;
if (a_k->count==0) { K = Kmax; }
else if (a_k->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "K must be positive" << endl; return 1; }
else { K = size_t(a_k->ival[0]); }
if (K>Kmax) { cerr << progstr+": " << __LINE__ << errstr << "K must be <= min(R,C)" << endl; return 1; }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be a matrix" << endl; return 1; }

//Set output header info
o1.F = o2.F = o3.F = i1.F;
o1.T = o3.T = i1.T;
o2.T = (i1.isreal()) ? i1.T : i1.T-100;
o1.R = i1.R;
o1.C = K;
o2.R = K;
o2.C = 1u;
o3.R = K;
o3.C = i1.C;
o1.S = o2.S = o3.S = i1.S;
o1.H = o2.H = o3.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X, *U, *S, *Vt;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { U = new float[o1.R*Kmax]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (U)" << endl; return 1; }
    try { S = new float[Kmax]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (S)" << endl; return 1; }
    try { Vt = new float[Kmax*o3.C]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 3 (Vt)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::svd_inplace_s(U,S,Vt,X,i1.R,i1.C,i1.iscolmajor(),K))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(U),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 1 (U)" << endl; return 1; }
    }
    if (wo2)
    {
        try { ofs2.write(reinterpret_cast<char*>(S),o2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 2 (S)" << endl; return 1; }
    }
    if (wo3)
    {
        try { ofs3.write(reinterpret_cast<char*>(Vt),o3.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 3 (Vt)" << endl; return 1; }
    }
    delete[] X; delete[] U; delete[] S; delete[] Vt;
}
else if (i1.T==101)
{
    float *X, *U, *S, *Vt;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { U = new float[2u*o1.R*Kmax]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (U)" << endl; return 1; }
    try { S = new float[Kmax]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (S)" << endl; return 1; }
    try { Vt = new float[2u*Kmax*o3.C]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 3 (Vt)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::svd_inplace_c(U,S,Vt,X,i1.R,i1.C,i1.iscolmajor(),K))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(U),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 1 (U)" << endl; return 1; }
    }
    if (wo2)
    {
        try { ofs2.write(reinterpret_cast<char*>(S),o2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 2 (S)" << endl; return 1; }
    }
    if (wo3)
    {
        try { ofs3.write(reinterpret_cast<char*>(Vt),o3.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 3 (Vt)" << endl; return 1; }
    }
    delete[] X; delete[] U; delete[] S; delete[] Vt;
}

//Finish
