//Includes
#include "matmul3t.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 3u, O = 1u;
int tr;

//Description
string descr;
descr += "Linear algebra function for 3 inputs.\n";
descr += "Matrix multiplication of inputs X1, X2, X3,\n";
descr += "with transpose taken of 3rd input (X3).\n";
descr += "\n";
descr += "Order of multiplication is chosen optimally:\n";
descr += "Y = (X1*X2)*X3', or Y = X1*(X2*X3').\n";
descr += "Thus, output Y has size R1xR3.\n";
descr += "\n";
descr += "Include -t (--transpose) to transpose X1 instead of X3.\n";
descr += "In this case: Y = X1'*X2*X3, and has size C1xC3.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ matmul3t X1 X2 X3 -o Y \n";
descr += "$ matmul3t X1 X2 X3 > Y \n";
descr += "$ cat X2 | matmul3t -t X1 - X3 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2,X3)");
struct arg_lit   *a_tr = arg_litn("t","transpose",0,1,"transpose X1 rather than X3 [default=false]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get tr
tr = (a_tr->count>0);

//Checks
if (i1.T!=i2.T || i1.T!=i3.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }
if (i3.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 3 (X3) found to be empty" << endl; return 1; }
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) must be a matrix" << endl; return 1; }
if (!i2.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) must be a matrix" << endl; return 1; }
if (!i3.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 3 (X3) must be a matrix" << endl; return 1; }
if (!major_compat(i1,i2) || !major_compat(i1,i3)) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same row/col major format unless vectors" << endl; return 1; }
if (tr && i1.R!=i2.R) { cerr << progstr+": " << __LINE__ << errstr << "R1 (nrows X1) must equal R2 (nrows X2) if transpose opt" << endl; return 1; }
if (tr && i2.C!=i3.R) { cerr << progstr+": " << __LINE__ << errstr << "C2 (ncols X2) must equal R3 (nrows X3) if transpose opt" << endl; return 1; }
if (!tr && i1.C!=i2.R) { cerr << progstr+": " << __LINE__ << errstr << "C1 (ncols X1) must equal R2 (nrows X2) if no transpose opt" << endl; return 1; }
if (!tr && i2.C!=i3.C) { cerr << progstr+": " << __LINE__ << errstr << "C2 (ncols X2) must equal C3 (ncols X3) if no transpose opt" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (tr) ? i1.C : i1.R;
o1.C = (tr) ? i3.C : i3.R;
o1.S = i1.S; o1.H = i1.H;

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
    if (codee::matmul3t_s(Y,X1,X2,X3,i1.R,i1.C,i3.R,i3.C,o1.iscolmajor(),tr))
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
    if (codee::matmul3t_c(Y,X1,X2,X3,i1.R,i1.C,i3.R,i3.C,o1.iscolmajor(),tr))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X1; delete[] X2; delete[] X3; delete[] Y;
}

//Finish
