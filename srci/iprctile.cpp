//Includes
#include "iprctile.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u};
const size_t I = 1u, O = 1u;
size_t dim;
double p;

//Description
string descr;
descr += "Vec2scalar (reduction) operation.\n";
descr += "Gets index of the pth percentile for each vector in X along dim.\n";
descr += "This is the index with the value closest to the pth percentile.\n";
descr += "The output is an integer, but data type is float or double.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
descr += "Use -d0 to get indices along cols.\n";
descr += "Use -d1 to get indices along rows.\n";
descr += "Use -d2 to get indices along slices.\n";
descr += "Use -d3 to get indices along hyperslices.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ iprctile X -o Y \n";
descr += "$ iprctile X > Y \n";
descr += "$ iprctile -d1 X > Y \n";
descr += "$ cat X | iprctile > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_dbl    *a_p = arg_dbln("p","p","<dbl>",0,1,"percentile in [0 100] [default=50]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get p
p = (a_p->count==0) ? 50.0 : a_p->dval[0];
if (p<0.0 || p>100.0) { cerr << progstr+": " << __LINE__ << errstr << "p must be in [0 100]" << endl; return 1; }

//Get dim
if (a_d->count==0) { dim = i1.isvec() ? i1.nonsingleton1() : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }

//Set output header info
o1.F = i1.F;
o1.T = i1.isreal() ? i1.T : i1.T-100u;
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
    if (codee::iprctile_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,float(p)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
