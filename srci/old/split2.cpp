//Includes
#include "split2.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 1u, O = 2u;
size_t dim;

//Description
string descr;
descr += "Splits input X into 2 equal-sized outputs Y1, Y2.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
descr += "Use -d0 to work along cols --> Y1 and Y2 have R/2 rows.\n";
descr += "Use -d1 to work along rows --> Y1 and Y2 have C/2 cols.\n";
descr += "Use -d2 to work along slices --> Y1 and Y2 have S/2 slices.\n";
descr += "Use -d3 to work along hyperslices --> Y1 and Y2 have H/2 hyperslices.\n";
descr += "\n";
descr += "For dim=0, num rows X must be even (R%2==0).\n";
descr += "For dim=1, num cols X must be even (C%2==0).\n";
descr += "For dim=2, num slices X must be even (S%2==0).\n";
descr += "For dim=3, num hyperslices X must be even (H%2==0).\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ split2 X -o Y1 -o Y2 \n";
descr += "$ split2 X -o Y1 > Y2 \n";
descr += "$ split2 -d1 X -o Y1 -o Y2 \n";
descr += "$ cat X | split2 -o Y1 -o Y2 \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output files (Y1,Y2)");

//Get options

//Get dim
if (a_d->count==0) { dim = i1.isrowvec() ? 1u : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
if (dim==0u && i1.R%2) { cerr << progstr+": " << __LINE__ << errstr << "num rows X must be even for dim=0" << endl; return 1; }
if (dim==1u && i1.C%2) { cerr << progstr+": " << __LINE__ << errstr << "num cols X must be even for dim=1" << endl; return 1; }
if (dim==2u && i1.S%2) { cerr << progstr+": " << __LINE__ << errstr << "num slices X must be even for dim=2" << endl; return 1; }
if (dim==3u && i1.H%2) { cerr << progstr+": " << __LINE__ << errstr << "num hyperslices X must be even for dim=3" << endl; return 1; }

//Set output header info
o1.F = o2.F = i1.F;
o1.T = o2.T = i1.T;
o1.R = o2.R = (dim==0u) ? i1.R/2 : i1.R;
o1.C = o2.C = (dim==1u) ? i1.C/2 : i1.C;
o1.S = o2.S = (dim==2u) ? i1.S/2 : i1.S;
o1.H = o2.H = (dim==3u) ? i1.H/2 : i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X, *Y1, *Y2;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y1 = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (Y1)" << endl; return 1; }
    try { Y2 = new float[o2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (Y2)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::split2_s(Y1,Y2,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y1),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 1 (Y1)" << endl; return 1; }
    }
    if (wo2)
    {
        try { ofs2.write(reinterpret_cast<char*>(Y2),o2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 2 (Y2)" << endl; return 1; }
    }
    delete[] X; delete[] Y1; delete[] Y2;
}
else if (i1.T==101u)
{
    float *X, *Y1, *Y2;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y1 = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (Y1)" << endl; return 1; }
    try { Y2 = new float[2u*o2.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (Y2)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::split2_c(Y1,Y2,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y1),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 1 (Y1)" << endl; return 1; }
    }
    if (wo2)
    {
        try { ofs2.write(reinterpret_cast<char*>(Y2),o2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 2 (Y2)" << endl; return 1; }
    }
    delete[] X; delete[] Y1; delete[] Y2;
}

//Finish
