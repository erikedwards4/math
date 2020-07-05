//Includes
#include "complex.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 2, O = 1;

//Description
string descr;
descr += "Elementwise function for 2 inputs with broadcasting.\n";
descr += "Makes complex output Y from real (X1) and imaginary (X2) parts.\n";
descr += "Both inputs (X1, X2) are real-valued.\n";
descr += "\n";
descr += "X1 and X2 must have the same size or broadcast-compatible sizes.\n";
descr += "Output (Y) has size max(R1,R2) x max(C1,C2).\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ complex X1 X2 -o Y \n";
descr += "$ complex X1 X2 > Y \n";
descr += "$ cat X1 | complex - X2 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2)");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }
if (!i2.isvec() && i2.iscolmajor()!=i1.iscolmajor()) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same row/col major format" << endl; return 1; }
if (i2.T!=i1.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (!bcast_compat(i1,i2)) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have same size or broadcast-compatible sizes" << endl; return 1; }
if (!i1.isreal() || !i2.isreal()) { cerr << progstr+": " << __LINE__ << errstr << "both inputs must be real-valued" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T + 100;
o1.R = (i1.R>i2.R) ? i1.R : i2.R;
o1.C = (i1.C>i2.C) ? i1.C : i2.C;
o1.S = (i1.S>i2.S) ? i1.S : i2.S;
o1.H = (i1.H>i2.H) ? i1.H : i2.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X1, *X2, *Y;
    try { X1 = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
    try { X2 = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
    if (codee::complex_s(Y,X1,X2,i1.R,i1.C,i1.S,i1.H,i2.R,i2.C,i2.S,i2.H,o1.iscolmajor()))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2; delete[] Y;
}

//Finish
