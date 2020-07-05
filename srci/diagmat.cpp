//Includes
#include "diagmat.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 1, O = 1;
int k;

//Description
string descr;
descr += "Matrix construct function.\n";
descr += "Puts vector X as the kth diagonal of matrix Y.\n";
descr += "\n";
descr += "Use -k (--k) to specify the diagonal number [default=0].\n";
descr += "\n";
descr += "For k=0, X becomes the main diagonal [default]. \n";
descr += "For k<0, X becomes the kth sub-diagonal. \n";
descr += "For k>0, X becomes the kth super-diagonal. \n";
descr += "\n";
descr += "Y has the minimal necessary size.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ diagmat X -o Y \n";
descr += "$ diagmat -k-1 X > Y \n";
descr += "$ cat X | diagmat -k2 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_k = arg_intn("k","k","<int>",0,1,"diagonal number [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get k
k = (a_k->count>0) ? a_k->ival[0] : 0;

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be a vector" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
if (k>0) { o1.R = i1.N(); o1.C = i1.N()+size_t(k); }
else { o1.R = i1.N()+size_t(-k); o1.C = i1.N(); }
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::diagmat_s(Y,X,o1.R,o1.C,o1.iscolmajor(),k))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}
else if (i1.T==101)
{
    float *X, *Y;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::diagmat_c(Y,X,o1.R,o1.C,o1.iscolmajor(),k))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
