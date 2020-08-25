//Includes
#include "logspace.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u,101u,102u};
const size_t I = 0u, O = 1u;
double a, b;
size_t dim;
size_t N;

//Description
string descr;
descr += "Generate function (0 inputs, 1 output).\n";
descr += "Output Y is a row vector with elements logarithmically spaced from 10^a to 10^b.\n";
descr += "\n";
descr += "Use -n (--N) to specify the length of Y [default=100].\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
descr += "Y is column vector for d=0, a row vector for d=1, etc.\n";
descr += "\n";
descr += "Use -t (--type) to specify output data type [default=1 (float)].\n";
descr += "Data type can also be 2 (double), 101 (complex float), or 102 (complex double).\n";
descr += "\n";
descr += "Use -f (--fmt) to specify output file format [default=147 -> NumPy].\n";
descr += "File format can also be 1 (ArrayFire), 65 (Armadillo),\n";
descr += "101 (minimal row-major format), or 102 (minimal col-major format).\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ logspace -n10 -a-1 -b2 -o Y \n";
descr += "$ logspace -n9 -a1 -b-3 > Y \n";
descr += "$ logspace -n11 -a-2.5 -b4.5 -t2 > Y \n";

//Argtable
struct arg_dbl    *a_a = arg_dbln("a","a","<dbl>",0,1,"start value [default=0.0]");
struct arg_dbl    *a_b = arg_dbln("b","b","<dbl>",0,1,"end value [default=1.0]");
struct arg_int    *a_n = arg_intn("n","N","<uint>",0,1,"length of output Y [default=100]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=1]");
struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=147]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get a
a = (a_a->count==0) ? 0.0 : a_a->dval[0];

//Get b
b = (a_b->count==0) ? 1.0 : a_b->dval[0];

//Get dim
if (a_d->count==0) { dim = 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get N
if (a_n->count==0) { N = 100u; }
else if (a_n->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "N (length of Y) must be positive" << endl; return 1; }
else { N = size_t(a_n->ival[0]); }

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

//Checks

//Set output header info
o1.R = (dim==0u) ? N : 1u;
o1.C = (dim==1u) ? N : 1u;
o1.S = (dim==2u) ? N : 1u;
o1.H = (dim==3u) ? N : 1u;

//Other prep

//Process
if (o1.T==1u)
{
    float *Y;
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (codee::logspace_s(Y,N,float(a),float(b)))
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
    if (codee::logspace_c(Y,N,float(a),float(b)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}

//Finish
