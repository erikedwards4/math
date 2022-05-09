//Includes
#include "partsort.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 1u, O = 1u;
size_t dim, k, L;
int a;

//Description
string descr;
descr += "Vec2vec operation.\n";
descr += "Partial sort of each vector in X along dim using quickselect algorithm.\n";
descr += "The output is sorted up to the kth element for each vector.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
descr += "Use -d0 to sort along cols.\n";
descr += "Use -d1 to sort along rows.\n";
descr += "Use -d2 to sort along slices.\n";
descr += "Use -d3 to sort along hyperslices.\n";
descr += "\n";
descr += "Include -a (--ascend) to sort in ascending order [default=descend].\n";
descr += "\n";
descr += "Use -k (--k) to specify k, the element up to which to sort [default=L-1].\n";
descr += "where 0 <= k < L, and L is the length of vectors in X.If k==0. \n";
descr += "If k==0, then no sorting is done. \n";
descr += "\n";
descr += "Examples:\n";
descr += "$ partsort -k14 X -o Y \n";
descr += "$ partsort -k21 -d1 X > Y \n";
descr += "$ cat X | partsort -k5 -d2 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
struct arg_int    *a_k = arg_intn("k","k","<uint>",0,1,"sort up to kth element [default=L-1]");
struct arg_lit    *a_a = arg_litn("a","ascend",0,1,"sort ascending [default=descending]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isvec() ? i1.nonsingleton1() : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get k
L = (dim==0u) ? i1.R : (dim==1u) ? i1.C : (dim==2u) ? i1.S : i1.H;
if (a_k->count==0) { k = L - 1u; }
else if (a_k->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "k must be nonnegative" << endl; return 1; }
else { k = size_t(a_k->ival[0]); }
if (k>=L) { cerr << progstr+": " << __LINE__ << errstr << "k must be 0 <= k < L" << endl; return 1; }

//Get a
a = (a_a->count>0);

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = i1.R; o1.C = i1.C; o1.S = i1.S; o1.H = i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X; //, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    //try { Y = new float[o1.N()]; }
    //catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    //if (codee::partsort_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,a,k))
    if (codee::partsort_inplace_s(X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,a,k))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; //delete[] Y;
}
else if (i1.T==101u)
{
    float *X; //, *Y;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    //try { Y = new float[2u*o1.N()]; }
    //catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    //if (codee::partsort_c(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,a,k))
    if (codee::partsort_inplace_c(X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,a,k))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; //delete[] Y;
}

//Finish
