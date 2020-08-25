//Includes
#include "e.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 0u, O = 1u;

//Description
string descr;
descr += "Generate function (0 inputs, 1 output).\n";
descr += "Output Y is a tensor (up to 4-D) filled with e (constant M_E).\n";
descr += "\n";
descr += "Use -t (--type) to specify output data type [default=1 -> float].\n";
descr += "Data type can also be 2 (double), 101 (float complex), 102 (double complex).\n";
descr += "\n";
descr += "Use -f (--fmt) to specify output file format [default=147 -> NumPy].\n";
descr += "File format can also be 1 (ArrayFire), 65 (Armadillo),\n";
descr += "101 (minimal row-major format), or 102 (minimal col-major format).\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ e -r2 -c3 -o Y \n";
descr += "$ e -r2 -c3 > Y \n";
descr += "$ e -r2 -c3 -s2 -t2 > Y \n";

//Argtable
struct arg_int   *a_nr = arg_intn("r","n_rows","<uint>",0,1,"num rows in output [default=1]");
struct arg_int   *a_nc = arg_intn("c","n_cols","<uint>",0,1,"num cols in output [default=1]");
struct arg_int   *a_ns = arg_intn("s","n_slices","<uint>",0,1,"num slices in output [default=1]");
struct arg_int   *a_nh = arg_intn("y","n_hyperslices","<uint>",0,1,"num hyperslices in output [default=1]");
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
else if (a_otyp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "data type must be positive int" << endl; return 1; }
else { o1.T = size_t(a_otyp->ival[0]); }
if ((o1.T==oktypes).sum()==0)
{
    cerr << progstr+": " << __LINE__ << errstr << "output data type must be in " << "{";
    for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1]) ? "}" : ","); }
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

//Checks
if (o1.H!=1u && o1.only_3D()) { cerr << progstr+": " << __LINE__ << errstr << "4D (hypercubes) not supported for Armadillo file format" << endl; return 1; }

//Set output header info

//Other prep

//Process
if (o1.T==1u)
{
    float *Y;
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (codee::e_s(Y,o1.N()))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}
else if (o1.T==101u)
{
    float *Y;
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (codee::e_c(Y,o1.N()))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}

//Finish
