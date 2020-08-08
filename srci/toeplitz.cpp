//Includes
#include "toeplitz.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 2, O = 1;

//Description
string descr;
descr += "Matrix construct function.\n";
descr += "Makes Toeplitz matrix Y from vec X1 (1st col of Y),\n";
descr += "and from vec X2 (1st row of Y).\n";
descr += "\n";
descr += "If X2 is empty, then default X2 = X1' is used.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ toeplitz X1 X2 -o Y \n";
descr += "$ toeplitz X1 X2 > Y \n";
descr += "$ cat X2 | toeplitz X1 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2)");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) must be a vector" << endl; return 1; }
if (!i2.isvec() && !i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) must be empty or a vector" << endl; return 1; }
if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = i1.N();
o1.C = (i2.N()>0u) ? i2.N() : i1.N();
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X1, *X2, *Y;
    try { X1 = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
    try { X2 = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
    try { Y = new float[o1.N()](); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
    if (i2.N()>0u)
    {
        if (codee::toeplitz2_s(Y,X1,X2,i1.N(),i2.N(),o1.iscolmajor()))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    }
    else
    {
        if (codee::toeplitz1_s(Y,X1,i1.N()))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    }
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
    try { Y = new float[2u*o1.N()](); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
    if (i2.N()>0u)
    {
        if (codee::toeplitz2_c(Y,X1,X2,i1.N(),i2.N(),o1.iscolmajor()))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    }
    else
    {
        if (codee::toeplitz1_c(Y,X1,i1.N()))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2; delete[] Y;
}

//Finish
