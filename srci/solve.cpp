//Includes
#include "solve.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 2u, O = 1u;
int tr;

//Description
string descr;
descr += "Linear algebra function for 2 inputs.\n";
descr += "Least-squares solution of linear system: A*X = B.\n";
descr += "The system can be over- or under-determined. \n";
descr += "\n";
descr += "This is similar to the \\ (backslash) operator in Octave,\n";
descr += "with X = A\\B. \n";
descr += "\n";
descr += "Input 1 (A) is usually a matrix.\n";
descr += "Input 2 (B) can be vector or matrix.\n";
descr += "The output (X) is solution to minimize norm2(B-A*X).\n";
descr += "\n";
descr += "Include -t (--transpose) to use A' (Hermitian transpose of A).\n";
descr += "In this case, minimize: norm2(B-A'*X). \n";
descr += "\n";
descr += "Examples:\n";
descr += "$ solve A B -o X \n";
descr += "$ solve A B > X \n";
descr += "$ cat A | solve -t - B > X \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (A,B)");
struct arg_lit   *a_tr = arg_litn("t","transpose",0,1,"transpose A");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (X)");

//Get options

//Get tr
tr = (a_tr->count>0);

//Checks
if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (A) found to be empty" << endl; return 1; }
if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (B) found to be empty" << endl; return 1; }
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (A) must be a matrix" << endl; return 1; }
if (!i2.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (B) must be a vector or matrix" << endl; return 1; }
if (!major_compat(i1,i2)) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same row/col major format unless vectors" << endl; return 1; }
if (tr && i1.C!=i2.R) { cerr << progstr+": " << __LINE__ << errstr << "C1 (nrows A) must equal R2 (nrows B) if transpose of A" << endl; return 1; }
if (!tr && i1.R!=i2.R) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same num rows if no transpose of A" << endl; return 1; }
if (tr && i1.T>100) { cerr << progstr+": " << __LINE__ << errstr << "transpose option not working for complex data types" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (tr) ? i1.R : i1.C;
o1.C = i2.C;
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *A, *B;
    try { A = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (A)" << endl; return 1; }
    try { B = new float[i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (B)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(A),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (A)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(B),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (B)" << endl; return 1; }
    if (codee::solve_inplace_s(A,B,i1.R,i1.C,i2.C,o1.iscolmajor(),tr))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(B),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (X)" << endl; return 1; }
    }
    delete[] A; delete[] B;
}
else if (i1.T==101u)
{
    float *A, *B;
    try { A = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (A)" << endl; return 1; }
    try { B = new float[2u*i2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (B)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(A),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (A)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(B),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (B)" << endl; return 1; }
    if (codee::solve_inplace_c(A,B,i1.R,i1.C,i2.C,o1.iscolmajor(),tr))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(B),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (X)" << endl; return 1; }
    }
    delete[] A; delete[] B;
}

//Finish
