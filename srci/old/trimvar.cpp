//Includes
#include "trimvar.c"

//Declarations
const valarray<size_t> oktypes = {1u,2u};
const size_t I = 1u, O = 1u;
size_t dim;
double p, q;
char b;

//Description
string descr;
descr += "Vec2scalar (reduction) operation.\n";
descr += "Gets trimmed (truncated) variance for each vector in X along dim,\n";
descr += "by trimming (excluding) the lower p% and upper q% of values.\n";
descr += "\n";
descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
descr += "Use -d0 to get trimmed var along cols.\n";
descr += "Use -d1 to get trimmed var along rows.\n";
descr += "Use -d2 to get trimmed var along slices.\n";
descr += "Use -d3 to get trimmed var along hyperslices.\n";
descr += "\n";
descr += "Use -p (--p) to give the lower percentage in [0 50) [default=10].\n";
descr += "Typical choices for p are 10 and 25.\n";
descr += "\n";
descr += "Use -q (--q) to give the upper percentage in [0 50) [default=p].\n";
descr += "Typically this is equal to p.\n";
descr += "\n";
descr += "The bottom p% and the top q% of values are excluded.\n";
descr += "These are not percentiles, just percentages of data to exclude,\n";
descr += "although they are approximately equal to the percentiles.\n";
descr += "\n";
descr += "Include -b (--biased) to use the biased denominator [default is unbiased].\n";
descr += "The biased denominator is N, and the unbiased denominator is N-1.\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ trimvar X -o Y \n";
descr += "$ trimvar -p25 X > Y \n";
descr += "$ trimvar -p5 -q0 -d1 X > Y \n";
descr += "$ cat X | trimvar > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
struct arg_dbl    *a_p = arg_dbln("p","p","<dbl>",0,1,"percentile in [0 50) [default=10]");
struct arg_dbl    *a_q = arg_dbln("q","q","<dbl>",0,1,"upper percentile in [0 50) [default=p]");
struct arg_lit    *a_b = arg_litn("b","biased",0,1,"use biased (N) denominator [default=N-1]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get p
p = (a_p->count==0) ? 10.0 : a_p->dval[0];
if (p<0.0 || p>=50.0) { cerr << progstr+": " << __LINE__ << errstr << "p must be in [0 50)" << endl; return 1; }

//Get q
q = (a_q->count==0) ? p : a_q->dval[0];
if (q<0.0 || q>=50.0) { cerr << progstr+": " << __LINE__ << errstr << "q must be in [0 50)" << endl; return 1; }

//Get b
b = (a_b->count>0);

//Get dim
if (a_d->count==0) { dim = i1.isrowvec() ? 1u : 0u; }
else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
else { dim = size_t(a_d->ival[0]); }
if (dim>3u) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = (dim==0u) ? 1u : i1.R;
o1.C = (dim==1u) ? 1u : i1.C;
o1.S = (dim==2u) ? 1u : i1.S;
o1.H = (dim==3u) ? 1u : i1.H;

//Other prep

//Process
if (i1.T==1u)
{
    float *X, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::trimvar_inplace_s(Y,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim,float(p),float(q),b))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish
