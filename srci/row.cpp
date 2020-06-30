//Includes
#include "row.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 1, O = 1;
int r;

//Description
string descr;
descr += "Matsel (matrix select) function.\n";
descr += "Gets one row of X as a row vector Y.\n";
descr += "If X is a 3D or 4D tensor, then \n";
descr += "Y is also a 3D or 4D tensor, but with 1 row.\n";
descr += "\n";
descr += "Use -r (--row) to specify the row number [default=0].\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ row X -o Y \n";
descr += "$ row -r2 X > Y \n";
descr += "$ cat X | row -r2 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_r = arg_intn("r","row","<uint>",0,1,"row number to get [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get r
if (a_r->count==0) { r = 0; }
else if (a_r->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "r must be nonnegative" << endl; return 1; }
else { r = a_r->ival[0]; }
if (r>=int(i1.R)) { cerr << progstr+": " << __LINE__ << errstr << "r must be int in [0 R-1]" << endl; return 1; }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = 1u; o1.C = i1.C;
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::row_inplace_s(X,int(i1.R),int(i1.C),int(i1.S),int(i1.H),i1.iscolmajor(),r))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X;
}
else if (i1.T==101)
{
    float *X;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::row_inplace_c(X,int(i1.R),int(i1.C),int(i1.S),int(i1.H),i1.iscolmajor(),r))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X;
}

//Finish

