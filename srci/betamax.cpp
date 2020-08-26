//Includes
#include "betamax.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u};
const size_t I = 1u, O = 1u;
size_t dim;
double base;

//Description
string descr;
descr += "Vec2vec function.\n";
descr += "\"Betamax\" function for each vector in X along dim.\n";
descr += "Output Y has the same size as X.\n";
descr += "\n";
descr += "This is just like softmax, except uses a different base (b).\n";
descr += "This is computed using exp(beta*X) instead of exp(X),\n";
descr += "where beta = log(b) and b = exp(beta). \n";
descr += "This is equivalent to using b^x instead of e^x.\n";
descr += "\n";
descr += "Use -b (--base) to specify the base [default=e].\n";
descr += "Be sure to enter b and not beta.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
descr += "Use -d0 to operate along cols.\n";
descr += "Use -d1 to operate along rows.\n";
descr += "Use -d2 to operate along slices.\n";
descr += "Use -d3 to operate along hyperslices.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ betamax -b2 X -o Y \n";
descr += "$ betamax -b3 -d1 X > Y \n";
descr += "$ cat X | betamax -b2.5 -d2 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_dbl    *a_b = arg_dbln("b","base","<dbl>",0,1,"base value [default=e]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get base
base = (a_b->count==0) ? exp(1.0) : a_b->dval[0];

//Get dim
if (a_d->count==0) { dim = i1.isvec() ? i1.nonsingleton1() : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

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
    if (codee::betamax_inplace_s(X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,float(base)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(X),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; //delete[] Y;
}

//Finish
