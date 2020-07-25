//Includes
#include "dot.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 2, O = 1;
size_t dim;

//Description
string descr;
descr += "Vecs2scalar function for 2 inputs with broadcasting.\n";
descr += "Gets dot product between each pair of vectors in X1, X2 along dim.\n";
descr += "Thus, output Y has length 1 along dim.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
descr += "Use -d0 to get dot product along cols.\n";
descr += "Use -d1 to get dot product along rows.\n";
descr += "Use -d2 to get dot product along slices.\n";
descr += "Use -d3 to get dot product along hyperslices.\n";
descr += "\n";
descr += "X1 and X2 must have the same size and data type, or\n";
descr += "either X1 or X2 can be a single vector of appropriate length.\n";
descr += "In that case, the single vector will be broadcast to the other input.\n";
descr += "\n";
descr += "For complex inputs, this uses the unconjugated dot product.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ dot X1 X2 -o Y \n";
descr += "$ dot -d1 X1 X2 > Y \n";
descr += "$ cat X1 | dot -d2 - X2 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Checks
if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }
if (!bcast_compat(i1,i2)) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have same size or broadcast-compatible sizes" << endl; return 1; }
if (!major_compat(i1,i2)) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same row/col major format unless vectors" << endl; return 1; }
size_t L1 = (dim==0) ? i1.R : (dim==1) ? i1.C : (dim==2) ? i1.S : i1.H;
size_t L2 = (dim==0) ? i2.R : (dim==1) ? i2.C : (dim==2) ? i2.S : i2.H;
if (L1!=L2) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have same length along dim" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0) ? 1 : (i1.R>i2.R) ? i1.R : i2.R;
o1.C = (dim==1) ? 1 : (i1.C>i2.C) ? i1.C : i2.C;
o1.S = (dim==2) ? 1 : (i1.S>i2.S) ? i1.S : i2.S;
o1.H = (dim==3) ? 1 : (i1.H>i2.H) ? i1.H : i2.H;

//Other prep

//Process
if (i1.T==1)
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
    if (codee::dot_s(Y,X1,X2,i1.R,i1.C,i1.S,i1.H,i2.R,i2.C,i2.S,i2.H,o1.iscolmajor(),dim))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2; delete[] Y;
}
else if (i1.T==101)
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
    if (codee::dot_c(Y,X1,X2,i1.R,i1.C,i1.S,i1.H,i2.R,i2.C,i2.S,i2.H,o1.iscolmajor(),dim))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2; delete[] Y;
}

//Finish
