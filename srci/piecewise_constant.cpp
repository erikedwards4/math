//Includes
#include <random>
#include <complex>

//Declarations
const valarray<size_t> oktypes = {1u,2u};
const size_t I = 2u, O = 1u;
random_device rd;    //random device to seed Mersenne twister engine
mt19937 mt_eng(rd()); 

//Description
string descr;
descr += "Random function: piecewise_constant_distribution.\n";
descr += "Outputs tensor of floats from a piecewise_constant distribution, \n";
descr += "with bin edges B and probabilities P.\n";
descr += "\n";
descr += "B must be a row or column vector with length N,\n";
descr += "with a monotonically-increasing set of floats.\n";
descr += "\n";
descr += "P must be a row or column vector with length N-1.\n";
descr += "The elements of P must be positive, but do not have to be probabilities.\n";
descr += "They are treated as weights and normalized to sum to 1.\n";
descr += "\n";
descr += "B and P must have the same data type, and Y will have this same data type.\n";
descr += "\n";
descr += "The output (Y) has floats in [B[0] B[N-1]]. \n";
descr += "\n";
descr += "Examples:\n";
descr += "$ piecewise_constant B P -r2 -c3 -o Y \n";
descr += "$ piecewise_constant B P -r2 -c3 > Y \n";
descr += "$ cat P | piecewise_constant B -r2 -c3 > Y \n";
descr += "$ cat B | piecewise_constant - P -r2 -c3 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (B,P)");
struct arg_int   *a_nr = arg_intn("r","R","<uint>",0,1,"num rows in output [default=1]");
struct arg_int   *a_nc = arg_intn("c","C","<uint>",0,1,"num cols in output [default=1]");
struct arg_int   *a_ns = arg_intn("s","S","<uint>",0,1,"num slices in output [default=1]");
struct arg_int   *a_nh = arg_intn("y","H","<uint>",0,1,"num hyperslices in the output [default=1]");
struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default from B]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get o1.F
if (a_ofmt->count==0) { o1.F = i1.F; }
else if (a_ofmt->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be nonnegative" << endl; return 1; }
else if (a_ofmt->ival[0]>255) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be < 256" << endl; return 1; }
else { o1.F = size_t(a_ofmt->ival[0]); }

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
if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs 1 and 2 must have the same data type" << endl; return 1; }
if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (B) must be a vector" << endl; return 1; }
if (!i2.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (P) must be a vector" << endl; return 1; }
if (i2.N()!=i1.N()-1u) { cerr << progstr+": " << __LINE__ << errstr << "length(P) must equal length(B)-1" << endl; return 1; }

//Set output header
o1.T = i1.T;

//Other prep

//Process
if (o1.T==1u)
{
    valarray<float> B(i1.N()), P(i2.N()), Y(o1.N());
    try { ifs1.read(reinterpret_cast<char*>(&B[0]),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading data for input file 1 (B)" << endl; return 1; }
    try { ifs2.read(reinterpret_cast<char*>(&P[0]),i2.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading data for input file 2 (P)" << endl; return 1; }
    if ((P<0.0f).sum()>0) { cerr << progstr+": " << __LINE__ << errstr << "elements of P must be nonnegative" << endl; return 1; }
    piecewise_constant_distribution<float> distr(&B[0],&B[i1.N()],&P[0]);
    try { generate_n(begin(Y),o1.N(),[&distr,&mt_eng](){return distr(mt_eng);}); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem during generate" << endl; return 1; }
    if (Y.size()!=o1.N()) { cerr << progstr+": " << __LINE__ << errstr << "unexpected output size" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(&Y[0]),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
}

//Finish
