//Includes
#include "diag.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 1, O = 1;
int k;

//Description
string descr;
descr += "Matsel (matrix select) function.\n";
descr += "Gets the kth diagonal of X as a column vector Y.\n";
descr += "\n";
descr += "Use -k (--k) to specify the diagonal [default=0].\n";
descr += "k must be in the range [1-R C-1].\n";
descr += "\n";
descr += "For k=0, the main diagonal is returned [default]. \n";
descr += "For k<0, the kth sub-diagonal is returned. \n";
descr += "For k>0, the kth super-diagonal is returned. \n";
descr += "\n";
descr += "Examples:\n";
descr += "$ diag X -o Y \n";
descr += "$ diag -k-1 X > Y \n";
descr += "$ cat X | diag -k2 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_k = arg_intn("k","k","<int>",0,1,"get kth diagonal [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get k
k = (a_k->count>0) ? a_k->ival[0] : 0;

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be a matrix" << endl; return 1; }
if (k>0 && k>=int(i1.C)) { cerr << progstr+": " << __LINE__ << errstr << "k must be int in [1-R C-1]" << endl; return 1; }
if (k<0 && -k>=int(i1.R)) { cerr << progstr+": " << __LINE__ << errstr << "k must be int in [1-R C-1]" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
if (k>=0) { o1.R = (int(i1.C)-k<int(i1.R)) ? i1.C-uint(k) : i1.R; }
else { o1.R = (int(i1.R)+k<int(i1.C)) ? i1.R-uint(-k) : i1.C; }
o1.C = 1u;
o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1)
{
    float *X; //*Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    //try { Y = new float[i1.N()]; }
    //catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    //if (codee::diag_s(Y,X,int(i1.R),int(i1.C),i1.iscolmajor(),k))
    if (codee::diag_inplace_s(X,int(i1.R),int(i1.C),i1.iscolmajor(),k))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; //delete[] Y;
}
else if (i1.T==101)
{
    float *X;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::diag_inplace_c(X,int(i1.R),int(i1.C),i1.iscolmajor(),k))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X;
}

//Finish

