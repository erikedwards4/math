//Includes
#include "linspace.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 0, O = 1;
double a, b;

//Description
string descr;
descr += "Generate function (0 inputs, 1 output).\n";
descr += "Output Y is a row vector with elements linearly spaced from a to b.\n";
descr += "\n";
descr += "Use -n (--N) to specify the length of Y [default=100].\n";
descr += "\n";
descr += "Use -t (--type) to specify output data type [default=1 -> float].\n";
descr += "Data type can also be 2 (double).\n";
descr += "Internally, the Sieve of Eratosthanes uses integers.\n";
descr += "\n";
descr += "Use -f (--fmt) to specify output file format [default=147 -> NumPy].\n";
descr += "File format can also be 1 (ArrayFire), 65 (Armadillo),\n";
descr += "101 (minimal row-major format), or 102 (minimal col-major format).\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ linspace -n10 -a0.1 -b9.1 -o Y \n";
descr += "$ linspace -n9 -a1 -b-3 > Y \n";
descr += "$ linspace -n11 -a-2.5 -b4.5 -t2 > Y \n";

//Argtable
struct arg_dbl    *a_a = arg_dbln("a","a","<dbl>",0,1,"start value [default=0.0]");
struct arg_dbl    *a_b = arg_dbln("b","b","<dbl>",0,1,"end value [default=1.0]");
struct arg_int    *a_n = arg_intn("n","N","<uint>",0,1,"length of output Y [default=100]");
struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=1]");
struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=147]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get a
a = (a_a->count==0) ? 0.0 : a_a->dval[0];

//Get b
b = (a_b->count==0) ? 1.0 : a_b->dval[0];

//Get o1.C
if (a_n->count==0) { o1.C = 100u; }
else if (a_n->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "N must be positive" << endl; return 1; }
else { o1.C = uint32_t(a_n->ival[0]); }

//Get o1.F
if (a_ofmt->count==0) { o1.F = 147; }
else if (a_ofmt->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be nonnegative" << endl; return 1; }
else if (a_ofmt->ival[0]>255) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be < 256" << endl; return 1; }
else { o1.F = uint8_t(a_ofmt->ival[0]); }

//Get o1.T
if (a_otyp->count==0) { o1.T = 1; }
else if (a_otyp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "data type must be positive int" << endl; return 1; }
else { o1.T = uint8_t(a_otyp->ival[0]); }
if ((o1.T==oktypes).sum()==0)
{
    cerr << progstr+": " << __LINE__ << errstr << "output data type must be in " << "{";
    for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1]) ? "}" : ","); }
    cerr << endl; return 1;
}

//Checks

//Set output header info
o1.R = o1.S = o1.H = 1u;

//Other prep

//Process
if (o1.T==1)
{
    float *Y;
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (codee::linspace_s(Y,int(o1.C),float(a),float(b)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}
else if (o1.T==101)
{
    float *Y;
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (codee::linspace_c(Y,int(o1.C),float(a),float(b)))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}

//Finish

