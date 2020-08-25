//Includes
#include "join3.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 3u, O = 1u;
size_t dim;

//Description
string descr;
descr += "Joins 3 inputs X1, X2, X3 into 1 output Y.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
descr += "Use -d0 to work along cols --> Y has R1+R2+R3 rows.\n";
descr += "Use -d1 to work along rows --> Y has C1+C2+C3 cols.\n";
descr += "Use -d2 to work along slices --> Y has S1+S2+S3 slices.\n";
descr += "Use -d3 to work along hyperslices --> Y has H1+H2+H3 hyperslices.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ join3 X1 X2 X3 -o Y \n";
descr += "$ join3 X1 X2 X3 > Y \n";
descr += "$ join3 -d1 X1 X2 X3 -o Y \n";
descr += "$ cat X1 | join3 - X2 X3 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2,X3)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Checks
if (i1.T!=i2.T || i1.T!=i3.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (!(i1.isvec() && i2.isvec() && i3.isvec()) && (i1.iscolmajor()!=i2.iscolmajor() || i1.iscolmajor()!=i2.iscolmajor()))
{ cerr << progstr+": " << __LINE__ << errstr << "all inputs must have the same row/col major format" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }
if (i3.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 3 (X3) found to be empty" << endl; return 1; }
if (dim!=0u && (i1.R!=i2.R || i1.R!=i3.R)) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have same num rows for dim!=0" << endl; return 1; }
if (dim!=1u && (i1.C!=i2.C || i1.C!=i3.C)) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have same num cols for dim!=1" << endl; return 1; }
if (dim!=2u && (i1.S!=i2.S || i1.S!=i3.S)) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have same num slices for dim!=2" << endl; return 1; }
if (dim!=3u && (i1.H!=i2.H || i1.H!=i3.H)) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have same num hyperslices for dim!=3" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0u) ? i1.R+i2.R+i3.R : i1.R;
o1.C = (dim==1u) ? i1.C+i2.C+i3.C : i1.C;
o1.S = (dim==2u) ? i1.S+i2.S+i3.S : i1.S;
o1.H = (dim==3u) ? i1.H+i2.H+i3.H : i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X1, *X2, *X3, *Y;
    try { X1 = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
    try { X2 = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
    try { X3 = new float[i3.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 3 (X3)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
    try { ifs3.read(reinterpret_cast<char*>(X3),i3.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 3 (X3)" << endl; return 1; }
    if (codee::join3_s(Y,X1,X2,X3,i1.R,i1.C,i1.S,i1.H,i2.R,i2.C,i2.S,i2.H,i3.R,i3.C,i3.S,i3.H,i1.iscolmajor(),dim))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2; delete[] X3; delete[] Y;
}
else if (i1.T==101u)
{
    float *X1, *X2, *X3, *Y;
    try { X1 = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
    try { X2 = new float[2u*i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
    try { X3 = new float[2u*i3.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 3 (X3)" << endl; return 1; }
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
    try { ifs3.read(reinterpret_cast<char*>(X3),i3.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 3 (X3)" << endl; return 1; }
    if (codee::join3_c(Y,X1,X2,X3,i1.R,i1.C,i1.S,i1.H,i2.R,i2.C,i2.S,i2.H,i3.R,i3.C,i3.S,i3.H,i1.iscolmajor(),dim))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2; delete[] X3; delete[] Y;
}

//Finish
