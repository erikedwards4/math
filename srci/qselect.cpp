//Includes
#include "qselect.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 1u, O = 1u;
size_t dim, k, L;
int largest;

//Description
string descr;
descr += "Does quick-select for each vector in X along dim.\n";
descr += "This selects the kth smallest element in each vector.\n";
descr += "\n";
descr += "Use -k (--k) to specify k [default=0].\n";
descr += "This must be a non-negative int (using 0-based indexing!),\n";
descr += "in the range [0 L-1], where L is the length of vecs along dim.\n";
descr += "\n";
descr += "Include -l (--largest) to get the kth largest instead.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
descr += "Use -d0 to multiply along cols.\n";
descr += "Use -d1 to multiply along rows.\n";
descr += "Use -d2 to multiply along slices.\n";
descr += "Use -d3 to multiply along hyperslices.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ qselect -k5 X -o Y \n";
descr += "$ qselect -k8 X > Y \n";
descr += "$ qselect -k4 -l -d1 X > Y \n";
descr += "$ cat X | qselect -k29 -l > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
struct arg_int    *a_k = arg_intn("k","k","<uint>",0,1,"get kth smallest element [default=0]");
struct arg_lit    *a_l = arg_litn("l","largest",0,1,"get kth largest [default=False -> kth smallest]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isvec() ? i1.nonsingleton1() : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get k
if (a_k->count==0) { k = 0u; }
else if (a_k->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "k must be non-negative" << endl; return 1; }
else { k = size_t(a_k->ival[0]); }

//Get l
largest = (a_l->count>0);

//Checks
L = (dim==0u) ? i1.R : (dim==1u) ? i1.C : (dim==2u) ? i1.S : i1.H;
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (k>=L) { cerr << progstr+": " << __LINE__ << errstr << "k must be in [0 L-1]" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0u) ? 1u : i1.R;
o1.C = (dim==1u) ? 1u : i1.C;
o1.S = (dim==2u) ? 1u : i1.S;
o1.H = (dim==3u) ? 1u : i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::qselect_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,k,largest))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}
else if (i1.T==101u)
{
    float *X, *Y;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::qselect_c(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,k,largest))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
