//Includes
#include <random>
#include <complex>
#include "pcg_random.hpp"

//Declarations
const valarray<size_t> oktypes = {1u,2u};
const size_t I = 1u, O = 1u;
valarray<double> P;  //Prepare different input data types for P
valarray<float> P1; valarray<long double> P3;
valarray<int8_t> P8; valarray<uint8_t> P9; valarray<int16_t> P16; valarray<uint8_t> P17;
valarray<int32_t> P32; valarray<int64_t> P64;

//Description
string descr;
descr += "Random function: discrete_distribution.\n";
descr += "Outputs tensor of uints from a discrete distribution, \n";
descr += "with N categories with probabilities P.\n";
descr += "\n";
descr += "P must be a row or column vector with length N.\n";
descr += "The elements of P must be positive, but do not have to be probabilities.\n";
descr += "They are treated as weights and normalized to sum to 1.0.\n";
descr += "\n";
descr += "The output (Y) has uints in [0 N-1], \n";
descr += "but can be requested to be of any real data type.\n";
descr += "Complex data types are not supported.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ discrete P -r2 -c3 -o Y \n";
descr += "$ discrete P -r2 -c3 > Y \n";
descr += "$ cat P | discrete -r2 -c3 > Y \n";
descr += "$ discrete P -r2 -c3 -t102 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (P)");
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
else if (a_otyp->ival[0]>103) { cerr << progstr+": " << __LINE__ << errstr << "output data type must be <= 103" << endl; return 1; }
else { o1.T = size_t(a_otyp->ival[0]); }
if ((o1.T==oktypes).sum()==0)
{
    cerr << progstr+": " << __LINE__ << errstr << "input data type must be in " << "{";
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
if (!i1.isvec()) { cerr << progstr+": " << __LINE__ << errstr << "input (P) must be a vector" << endl; return 1; }

//Set output header

//Other prep
//The top 2 lines would be for C++ random (no PCG)
//random_device rd;  //random device to seed Mersenne twister engine
//mt19937 mt_eng(rd());
pcg_extras::seed_seq_from<std::random_device> seed_source;
//these are a list of worthwhile generator engines (c is for cryptographic security, i.e. lowest predictability, but k is for equidistribution)
//pcg32 pcg_eng(seed_source);            //if need fast default generator
//pcg64 pcg_eng(seed_source);            //64-bit generator, 2^128 period, 2^127 streams
//pcg64_unique pcg_eng(seed_source);     //64-bit generator, 2^128 period, every instance has its own unique stream
//pcg32_k64 pcg_eng(seed_source);        //32-bit 64-dimensionally equidistributed generator, 2^2112 period, 2^63 streams (about the same state size and period as arc4random)
pcg64_k1024 pcg_eng(seed_source);      //64-bit 64-dimensionally equidistributed generator, 2^65664 period, 2^63 streams (larger period than the mersenne twister)
//pcg64_c1024 pcg_eng(seed_source);      //64-bit generator, 2^65664 period, 2^63 streams; uniform but not equidistributed; harder to predict than the above generator

//Get P
P.resize(i1.N(),0.0);
if (i1.T==1u)
{
    P1.resize(i1.N());
    try { ifs1.read(reinterpret_cast<char*>(&P1[0]),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (P)" << endl; return 1; }
    for (size_t n=0u; n<i1.N(); n++) { P[n] = double(P1[n]); }
}
else if (i1.T==2)
{
    try { ifs1.read(reinterpret_cast<char*>(&P[0]),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (P)" << endl; return 1; }
}
else
{
    cerr << progstr+": " << __LINE__ << errstr << "data type (" << int(i1.T) << ") not supported" << endl; return 1;
}
if ((P<0.0).sum()>0) { cerr << progstr+": " << __LINE__ << errstr << "elements of P must be non-negative" << endl; return 1; }

//Process
if (o1.T==1u)
{
    valarray<int> Yi(o1.N()); valarray<float> Y(o1.N());
    discrete_distribution<uint> distr(&P[0],&P[i1.N()]);
    try { generate_n(begin(Yi),o1.N(),[&distr,&pcg_eng](){return distr(pcg_eng);}); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem during generate" << endl; return 1; }
    for (size_t n=0u; n<o1.N(); ++n) { Y[n] = float(Yi[n]); }
    if (Y.size()!=o1.N()) { cerr << progstr+": " << __LINE__ << errstr << "unexpected output size" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(&Y[0]),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
}

//Finish
