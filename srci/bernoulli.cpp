//Includes
#include <random>

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 0, O = 1;
double p;
random_device rd;  //random device to seed Mersenne twister engine
mt19937 mt_eng(rd());

//Description
string descr;
descr += "Random function: bernoulli_distribution.\n";
descr += "Outputs tensor of ints in {0,1} from a Bernoulli distribution,\n";
descr += "with parameter p (probability of true outcome).\n";
descr += "\n";
descr += "The output data is in {0,1}, since boolean, but any output data type can be used.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ bernoulli -p0.9 -r2 -c3 -o Y \n";
descr += "$ bernoulli -p0.9 -r2 -c3 > Y \n";
descr += "$ bernoulli -p0.9 -r2 -c3 -t2 > Y \n";

//Argtable
struct arg_dbl    *a_p = arg_dbln("p","p","<dbl>",0,1,"probability of generating true [default=0.5]");
struct arg_int   *a_nr = arg_intn("r","R","<uint>",0,1,"num rows in output [default=1]");
struct arg_int   *a_nc = arg_intn("c","C","<uint>",0,1,"num cols in output [default=1]");
struct arg_int   *a_ns = arg_intn("s","S","<uint>",0,1,"num slices in output [default=1]");
struct arg_int   *a_nh = arg_intn("y","H","<uint>",0,1,"num hyperslices in the output [default=1]");
struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=1]");
struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=147]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get o1.F
if (a_ofmt->count==0) { o1.F = 102; }
else if (a_ofmt->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be nonnegative" << endl; return 1; }
else if (a_ofmt->ival[0]>255) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be < 256" << endl; return 1; }
else { o1.F = uint8_t(a_ofmt->ival[0]); }

//Get o1.T
if (a_otyp->count==0) { o1.T = 1; }
else if (a_otyp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "output data type must be positive" << endl; return 1; }
else { o1.T = uint8_t(a_otyp->ival[0]); }
if ((o1.T==oktypes).sum()==0)
{
    cerr << progstr+": " << __LINE__ << errstr << "input data type must be in " << "{";
    for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1]) ? "}" : ","); }
    cerr << endl; return 1;
}

//Get o1.R
if (a_nr->count==0) { o1.R = 1u; }
else if (a_nr->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "R (nrows) must be nonnegative" << endl; return 1; }
else { o1.R = uint32_t(a_nr->ival[0]); }

//Get o1.C
if (a_nc->count==0) { o1.C = 1u; }
else if (a_nc->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "C (ncols) must be nonnegative" << endl; return 1; }
else { o1.C = uint32_t(a_nc->ival[0]); }

//Get o1.S
if (a_ns->count==0) { o1.S = 1u; }
else if (a_ns->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "S (nslices) must be nonnegative" << endl; return 1; }
else { o1.S = uint32_t(a_ns->ival[0]); }

//Get o1.H
if (a_nh->count==0) { o1.H = 1u; }
else if (a_nh->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "H (nhyperslices) must be nonnegative" << endl; return 1; }
else { o1.H = uint32_t(a_nh->ival[0]); }

//Get p
p = (a_p->count>0) ? a_p->dval[0] : 0.5;
if (p<0.0) { cerr << progstr+": " << __LINE__ << errstr << "p must be >= 0.0" << endl; return 1; }
else if (p>1.0) { cerr << progstr+": " << __LINE__ << errstr << "p must be <= 1.0" << endl; return 1; }

//Checks

//Set output header

//Other prep
bernoulli_distribution distr(p);

//Process
if (o1.T==1)
{
    valarray<float> Y(o1.N());
    try { generate_n(begin(Y),o1.N(),[&distr,&mt_eng](){return distr(mt_eng);}); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem during generate" << endl; return 1; }
    if (Y.size()!=o1.N()) { cerr << progstr+": " << __LINE__ << errstr << "unexpected output size" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(&Y[0]),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
}
else if (o1.T==101)
{
    valarray<complex<float>> Y(o1.N());
    try { generate_n(begin(Y),o1.N(),[&distr,&mt_eng](){complex<float> y(distr(mt_eng),distr(mt_eng)); return y;}); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem during generate" << endl; return 1; }
    if (Y.size()!=o1.N()) { cerr << progstr+": " << __LINE__ << errstr << "unexpected output size" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(&Y[0]),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
}

//Finish

