//Includes
#include "tril.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 1, O = 1;
int k;

//Description
string descr;
descr += "Matrix construct function.\n";
descr += "Makes lower-triangular Y from input matrix X.\n";
descr += "This zeros all elements above the kth diagonal.\n";
descr += "\n";
descr += "Use -k (--k) to specify the diagonal number [default=0].\n";
descr += "\n";
descr += "For k=0, Y is 0 above the main diagonal [default]. \n";
descr += "For k<0, Y is 0 above the kth sub-diagonal. \n";
descr += "For k>0, Y is 0 above the kth super-diagonal. \n";
descr += "\n";
descr += "For tensor X, each matrix slice is processed identically.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ tril X -o Y \n";
descr += "$ tril -k-1 X > Y \n";
descr += "$ cat X | tril -k2 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_k = arg_intn("k","k","<int>",0,1,"diagonal number [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get k
k = (a_k->count>0) ? a_k->ival[0] : 0;

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (k<=-int(i1.R) || k>=int(i1.C)) { cerr << progstr+": " << __LINE__ << errstr << "k must be in [1-R C-1]" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = i1.R; o1.C = i1.C; o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::tril_inplace_s(X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),k))
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
    float *X; //, *Y;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    //try { Y = new float[2u*o1.N()]; }
    //catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    //if (codee::tril_c(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),k))
    if (codee::tril_inplace_c(X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),k))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; //delete[] Y;
}

//Finish
// if (i1.T==1)
// {
//     float *X, *Y;
//     try { X = new float[i1.N()]; }
//     catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
//     try { Y = new float[o1.N()]; }
//     catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
//     try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
//     catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
//     if (codee::tril_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),k))
//     { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
//     if (wo1)
//     {
//         try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
//         catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
//     }
//     delete[] X; delete[] Y;
// }
