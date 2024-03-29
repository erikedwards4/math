//Includes
#include <cfloat>
#include "randn.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 0u, O = 1u;
double mn, std;

//Description
string descr;
descr += "Normally-distributed random floats.\n";
descr += "Outputs tensor of floats from a normal distribution with mean and stddev. \n";
descr += "\n";
descr += "This uses modified code from PCG random, but does not require it to be installed.\n";
descr += "\n";
descr += "For complex output, real/imag parts are separately set using the same params.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ randn -r2 -c3 -o Y \n";
descr += "$ randn -m4.1 -d2.5 -r2 -c3 -t1 > Y \n";
descr += "$ randn -m-1.5 -d3 -r2 -c3 -t102 > Y \n";

//Argtable
struct arg_dbl   *a_mn = arg_dbln("m","mean","<dbl>",0,1,"mean parameter [default=0.0]");
struct arg_dbl  *a_std = arg_dbln("d","stddev","<dbl>",0,1,"std dev parameter [default=1.0]");
struct arg_int   *a_nr = arg_intn("r","R","<uint>",0,1,"num rows in output [default=1]");
struct arg_int   *a_nc = arg_intn("c","C","<uint>",0,1,"num cols in output [default=1]");
struct arg_int   *a_ns = arg_intn("s","S","<uint>",0,1,"num slices in output [default=1]");
struct arg_int   *a_nh = arg_intn("y","H","<uint>",0,1,"num hyperslices in the output [default=1]");
struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=1]");
struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=147]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get o1.F
if (a_ofmt->count==0) { o1.F = 147u; }
else if (a_ofmt->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be nonnegative" << endl; return 1; }
else if (a_ofmt->ival[0]>255) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be < 256" << endl; return 1; }
else { o1.F = size_t(a_ofmt->ival[0]); }

//Get o1.T
if (a_otyp->count==0) { o1.T = 1u; }
else if (a_otyp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "output data type must be positive" << endl; return 1; }
else { o1.T = size_t(a_otyp->ival[0]); }
if ((o1.T==oktypes).sum()==0)
{
    cerr << progstr+": " << __LINE__ << errstr << "output data type must be in " << "{";
    for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1u]) ? "}" : ","); }
    cerr << endl; return 1;
}

//Get o1.R
if (a_nr->count==0) { o1.R = 1u; }
else if (a_nr->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "R (nrows) must be nonnegative" << endl; return 1; }
else { o1.R = size_t(a_nr->ival[0]); }

//Get o1.C
if (a_nc->count==0) { o1.C = 1u; }
else if (a_nc->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "C (ncols) must be nonnegative" << endl; return 1; }
else { o1.C = size_t(a_nc->ival[0]); }

//Get o1.S
if (a_ns->count==0) { o1.S = 1u; }
else if (a_ns->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "S (nslices) must be nonnegative" << endl; return 1; }
else { o1.S = size_t(a_ns->ival[0]); }

//Get o1.H
if (a_nh->count==0) { o1.H = 1u; }
else if (a_nh->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "H (nhyperslices) must be nonnegative" << endl; return 1; }
else { o1.H = size_t(a_nh->ival[0]); }

//Get mn
mn = (a_mn->count>0) ? a_mn->dval[0] : 0.0;
if (o1.T==1 && mn>=double(FLT_MAX)) { cerr << progstr+": " << __LINE__ << errstr << "m (mean) must be < " << double(FLT_MAX) << endl; return 1; }
if (o1.T==1 && mn<=-double(FLT_MAX)) { cerr << progstr+": " << __LINE__ << errstr << "m (mean) must be > " << -double(FLT_MAX) << endl; return 1; }

//Get std
std = (a_std->count>0) ? a_std->dval[0] : 1.0;
if (std<0.0) { cerr << progstr+": " << __LINE__ << errstr << "s (stddev) must be nonnegative " << endl; return 1; }
if (o1.T==1 && std>=double(FLT_MAX)) { cerr << progstr+": " << __LINE__ << errstr << "s (stddev) must be < " << double(FLT_MAX) << endl; return 1; }

//Checks

//Set output header

//Other prep

//Process
if (o1.T==1u)
{
    float *Y;
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (codee::randn_s(Y,(float)mn,(float)std,o1.N()))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(&Y[0]),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}
else if (o1.T==101u)
{
    float *Y;
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (codee::randn_c(Y,(float)mn,(float)std,o1.N()))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(&Y[0]),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}

//Finish
