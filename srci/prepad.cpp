//Includes
#include "prepad.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 1u, O = 1u;
size_t dim, P;
double val;

//Description
string descr;
descr += "Vec2vec operation.\n";
descr += "Prepads each vector in X with P elements equal to val.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
descr += "Use -d0 to get prepad along cols.\n";
descr += "Use -d1 to get prepad along rows.\n";
descr += "Use -d2 to get prepad along slices.\n";
descr += "Use -d3 to get prepad along hyperslices.\n";
descr += "\n";
descr += "Use -v (--val) to specify the fill value [default=0.0].\n";
descr += "\n";
descr += "Use -p (--P) to give the prepad length [default=0].\n";
descr += "The length of output Y along dim will be greater than X by P elements.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ prepad X -p10 -o Y \n";
descr += "$ prepad X -p10 -v1.5 > Y \n";
descr += "$ prepad X -p16 -d1 > Y \n";
descr += "$ cat X | prepad -p5 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
struct arg_int    *a_p = arg_intn("p","P","<uint>",0,1,"prepad length [default=0]");
struct arg_dbl   *a_v = arg_dbln("v","val","<dbl>",0,1,"fill value [default=0.0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isvec() ? i1.nonsingleton1() : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get P
if (a_p->count==0) { P = 0u; }
else if (a_p->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "P (prepad length) must be nonnegative" << endl; return 1; }
else { P = size_t(a_p->ival[0]); }

//Get val
val = (a_v->count==0) ? 0.0 : a_v->dval[0];

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) found to be empty" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0u) ? i1.R+P : i1.R;
o1.C = (dim==1u) ? i1.C+P : i1.C;
o1.S = (dim==2u) ? i1.S+P : i1.S;
o1.H = (dim==3u) ? i1.H+P : i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    if (codee::prepad_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,P,float(val)))
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
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X)" << endl; return 1; }
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X)" << endl; return 1; }
    if (codee::prepad_c(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,P,float(val)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
