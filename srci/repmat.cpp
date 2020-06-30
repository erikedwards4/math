//Includes
#include "repmat.c"

//Declarations
const valarray<uint8_t> oktypes = {1,2,101,102};
const size_t I = 1, O = 1;
uint32_t Nr, Nc, Ns, Nh;

//Description
string descr;
descr += "Matrix construct function.\n";
descr += "Takes input tensor X and repeats it Nr x Nc x Ns x Nh times,\n";
descr += "such that the output Y has size R*Nr x C*Nc x S*Ns x H*Nh.\n";
descr += "\n";
descr += "Use -r (--Nr) to specify the number of row repetitions [default=1].\n";
descr += "Use -c (--Nc) to specify the number of col repetitions [default=1].\n";
descr += "Use -s (--Ns) to specify the number of slice repetitions [default=1].\n";
descr += "Use -y (--Nh) to specify the number of hyperslice repetitions [default=1].\n";
descr += "\n";
descr += "Examples:\n";
descr += "$ repmat X -r2 -c3 -o Y \n";
descr += "$ repmat X -s3 > Y \n";
descr += "$ cat X | repmat -r3 -c2 > Y \n";

//Argtable
struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
struct arg_int   *a_nr = arg_intn("r","Nr","<int>",0,1,"num row repetitions [default=1]");
struct arg_int   *a_nc = arg_intn("c","Nc","<int>",0,1,"num col repetitions [default=1]");
struct arg_int   *a_ns = arg_intn("s","Ns","<int>",0,1,"num slice repetitions [default=1]");
struct arg_int   *a_nh = arg_intn("y","Nh","<int>",0,1,"num hyperslice repetitions [default=1]");
struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");

//Get options

//Get Nr
if (a_nr->count==0) { Nr = 1u; }
else if (a_nr->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "Nr (num row reps) must be positive" << endl; return 1; }
else { Nr = uint32_t(a_nr->ival[0]); }

//Get Nc
if (a_nc->count==0) { Nc = 1u; }
else if (a_nc->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "Nc (num col reps) must be positive" << endl; return 1; }
else { Nc = uint32_t(a_nc->ival[0]); }

//Get Ns
if (a_ns->count==0) { Ns = 1u; }
else if (a_ns->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "Ns (num slice reps) must be positive" << endl; return 1; }
else { Ns = uint32_t(a_ns->ival[0]); }

//Get Nh
if (a_nh->count==0) { Nh = 1u; }
else if (a_nh->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "Nh (num hyperslice reps) must be positive" << endl; return 1; }
else { Nh = uint32_t(a_nh->ival[0]); }

//Checks
if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }

//Set output header info
o1.F = i1.F; o1.T = i1.T;
o1.R = i1.R * Nr; o1.C = i1.C * Nc;
o1.S = i1.S * Ns; o1.H = i1.H * Nh;

//Other prep

//Process
if (i1.T==1)
{
    float *X, *Y;
    try { X = new float[i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::repmat_s(Y,X,int(i1.R),int(i1.C),int(i1.S),int(i1.H),int(Nr),int(Nc),int(Ns),int(Nh),i1.iscolmajor()))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}
else if (i1.T==101)
{
    float *X, *Y;
    try { X = new float[2u*i1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
    try { Y = new float[2u*o1.N()]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
    try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
    if (codee::repmat_c(Y,X,int(i1.R),int(i1.C),int(i1.S),int(i1.H),int(Nr),int(Nc),int(Ns),int(Nh),i1.iscolmajor()))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    if (wo1)
    {
        try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
    }
    delete[] X; delete[] Y;
}

//Finish

