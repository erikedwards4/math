//Includes
#include "primes.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u};
const size_t I = 0u, O = 1u;
size_t P, cnt;
size_t *X;

//Description
string descr;
descr += "Generate function (0 inputs, 1 output).\n";
descr += "Output Y is a row vector filled with prime numbers from 2 through P.\n";
descr += "\n";
descr += "Use -p (--P) to specify the largest integer to return [default=2].\n";
descr += "That is, all primes less than or equal to P will be returned.\n";
descr += "\n";
descr += "Use -t (--type) to specify output data type [default=1 -> float].\n";
descr += "Data type can also be 2 (double).\n";
descr += "Internally, the Sieve of Eratosthanes uses size_t integers.\n";
descr += "\n";
descr += "Use -f (--fmt) to specify output file format [default=147 -> NumPy].\n";
descr += "File format can also be 1 (ArrayFire), 65 (Armadillo),\n";
descr += "101 (minimal row-major format), or 102 (minimal col-major format).\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ primes -p100 -o Y \n";
descr += "$ primes -p100 > Y \n";
descr += "$ primes -p100 -t2 > Y \n";

//Argtable
struct arg_int    *a_p = arg_intn("p","P","<uint>",0,1,"return primes <= P [default=2]");
struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=1]");
struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=147]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get P
if (a_p->count==0) { P = 2u; }
else if (a_p->ival[0]<2) { cerr << progstr+": " << __LINE__ << errstr << "P must be greater than 1" << endl; return 1; }
else { P = size_t(a_p->ival[0]); }

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
    for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1u]) ? "}" : ","); }
    cerr << endl; return 1;
}

//Checks

//Set output header info
try { X = new size_t[1u+(P-1u)/2u]; }
catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for integer file (X)" << endl; return 1; }
if (codee::primes_u(X,&cnt,(P-1u)/2u))
{ cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
o1.C = cnt;
o1.R = o1.S = o1.H = 1u;

//Other prep

//Process
if (o1.T==1u)
{
    float *Y;
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    for (size_t c=0u; c<o1.C; ++c) { Y[c] = float(X[c]); }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}

//Finish
