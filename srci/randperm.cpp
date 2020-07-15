//Includes
#include "randperm.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2};
const size_t I = 0, O = 1;
size_t N, M;
size_t dim;

//Description
string descr;
descr += "Random function: random permutation.\n";
descr += "Outputs a row vector of M ints from 1 to N in random order.\n";
descr += "\n";
descr += "Use -n (--N) to specify the highest integer for the random shuffle [default=1]\n";
descr += "\n";
descr += "Use -m (--M) to specify the length of the output vector [default=N].\n";
descr += "If M is given, it must be less than N, and only the first M ints are output.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
descr += "Y is column vector for d=0, a row vector for d=1, etc.\n";
descr += "\n";
descr += "Output is ints, but data type is float or double.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ randperm -n9 -o Y \n";
descr += "$ randperm -n1000 > Y \n";
descr += "$ randperm -n50 -t2 > Y \n";

//Argtable
struct arg_int    *a_n = arg_intn("n","N","<uint>",0,1,"highest int in random shuffle [default=1]");
struct arg_int    *a_m = arg_intn("m","M","<uint>",0,1,"num elements in output [default=1]");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=1]");
struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=147]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get dim
if (a_d->count==0) { dim = 0; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Get o1.F
if (a_ofmt->count==0) { o1.F = 147; }
else if (a_ofmt->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be nonnegative" << endl; return 1; }
else if (a_ofmt->ival[0]>255) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be < 256" << endl; return 1; }
else { o1.F = uint8_t(a_ofmt->ival[0]); }

//Get o1.T
if (a_otyp->count==0) { o1.T = 1; }
else if (a_otyp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "output data type must be positive" << endl; return 1; }
else { o1.T = uint8_t(a_otyp->ival[0]); }
if ((o1.T==oktypes).sum()==0)
{
    cerr << progstr+": " << __LINE__ << errstr << "output data type must be in " << "{";
    for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1]) ? "}" : ","); }
    cerr << endl; return 1;
}

//Get N
if (a_n->count==0) { N = 1; }
else if (a_n->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "N must be nonnegative" << endl; return 1; }
else { N = size_t(a_n->ival[0]); }

//Get M
if (a_m->count==0) { M = N; }
else if (a_m->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "M (output length) must be nonnegative" << endl; return 1; }
else { M = size_t(a_m->ival[0]); }

//Checks
if (M>N) { cerr << progstr+": " << __LINE__ << errstr << "M must be less than or equal to N" << endl; return 1; }

//Set output header
o1.R = (dim==0) ? M : 1u;
o1.C = (dim==1) ? M : 1u;
o1.S = (dim==2) ? M : 1u;
o1.H = (dim==3) ? M : 1u;

//Other prep

//Process
if (o1.T==1)
{
    float *Y;
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    if (codee::randperm_s(Y,M,N))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] Y;
}

//Finish
