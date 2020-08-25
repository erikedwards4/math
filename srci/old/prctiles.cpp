//Includes
#include "prctiles.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u};
const size_t I = 2u, O = 1u;
size_t dim;

//Description
string descr;
descr += "Vec2vec operation.\n";
descr += "Gets percentiles in P along dim of X.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
descr += "Use -d0 to get percentiles along cols.\n";
descr += "Use -d1 to get percentiles along rows.\n";
descr += "Use -d2 to get percentiles along slices.\n";
descr += "Use -d3 to get percentiles along hyperslices.\n";
descr += "\n";
descr += "The 2nd input P is a vector of length Ly.\n";
descr += "The output Y has length Ly along dim, but otherwise the same size as X.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ prctiles X P -o Y \n";
descr += "$ prctiles X P > Y \n";
descr += "$ prctiles -d1 X P > Y \n";
descr += "$ cat P | prctiles X > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X,P)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isrowvec() ? 1u : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (P) found to be empty" << endl; return 1; }
if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (P) must be a vector" << endl; return 1; }
if (i2.T!=i1.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0u) ? i2.N() : i1.R;
o1.C = (dim==1u) ? i2.N() : i1.C;
o1.S = (dim==2u) ? i2.N() : i1.S;
o1.H = (dim==3u) ? i2.N() : i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X, *P, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { P = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (P)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(P),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (P)" << endl; return 1; }
    if (codee::prctiles_s(Y,X,P,i2.N(),i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim))
    //if (codee::prctiles_inplace_s(Y,X,P,i2.N(),i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] P; delete[] Y;
}

//Finish
