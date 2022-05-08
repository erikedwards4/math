//Includes
#include "qsorti.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 1u, O = 1u;
size_t dim;
int a;

//Description
string descr;
descr += "Vec2vec operation.\n";
descr += "Sorts elements of X along dim and returns the indices.\n";
descr += "These are in [0 L-1] and are the inverse of the ranks (see ranks).\n";
descr += "That is, the qsorti indices would recover X from the sorted X.\n";
descr += "\n";
descr += "This uses (my own implementation of) the quicksort algorithm.\n";
descr += "This is faster than insert_sorti if the vectors are long (e.g., L > 128).\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
descr += "Use -d0 to sort along cols.\n";
descr += "Use -d1 to sort along rows.\n";
descr += "Use -d2 to sort along slices.\n";
descr += "Use -d3 to sort along hyperslices.\n";
descr += "\n";
descr += "Include -a (--ascend) to sort in ascending order [default=descend].\n";
descr += "\n";
descr += "Although indices are ints, this returns float or double data type.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ qsorti X -o Y \n";
descr += "$ qsorti -d1 X > Y \n";
descr += "$ cat X | qsorti -d2 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
struct arg_lit    *a_a = arg_litn("a","ascend",0,1,"sort ascending [default=descending]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isvec() ? i1.nonsingleton1() : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get a
a = (a_a->count>0);

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }

//Set output header info
o1.F = i1.F;
o1.T = i1.isreal() ? i1.T : i1.T-100u;
o1.R = i1.R; o1.C = i1.C;
o1.S = i1.S; o1.H = i1.H;

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
    if (codee::qsorti_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,a))
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
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::qsorti_c(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,a))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
