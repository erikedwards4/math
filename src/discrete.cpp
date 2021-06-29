//@author Erik Edwards
//@date 2018-present
//@license BSD 3-clause


#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <cstring>
#include <valarray>
#include <unordered_map>
#include <argtable2.h>
#include "../util/cmli.hpp"
#include <random>
#include <complex>
#include "pcg_random.hpp"

#ifdef I
#undef I
#endif


int main(int argc, char *argv[])
{
    using namespace std;


    //Declarations
    int ret = 0;
    const string errstr = ": \033[1;31merror:\033[0m ";
    const string warstr = ": \033[1;35mwarning:\033[0m ";
    const string progstr(__FILE__,string(__FILE__).find_last_of("/")+1,strlen(__FILE__)-string(__FILE__).find_last_of("/")-5);
    const valarray<size_t> oktypes = {1u,2u};
    const size_t I = 1u, O = 1u;
    ifstream ifs1; ofstream ofs1;
    int8_t stdi1, stdo1, wo1;
    ioinfo i1, o1;
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
    int nerrs;
    struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (P)");
    struct arg_int   *a_nr = arg_intn("r","R","<uint>",0,1,"num rows in output [default=1]");
    struct arg_int   *a_nc = arg_intn("c","C","<uint>",0,1,"num cols in output [default=1]");
    struct arg_int   *a_ns = arg_intn("s","S","<uint>",0,1,"num slices in output [default=1]");
    struct arg_int   *a_nh = arg_intn("y","H","<uint>",0,1,"num hyperslices in the output [default=1]");
    struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=1]");
    struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=147]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_fi, a_nr, a_nc, a_ns, a_nh, a_otyp, a_ofmt, a_fo, a_help, a_end};
    if (arg_nullcheck(argtable)!=0) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating argtable" << endl; return 1; }
    nerrs = arg_parse(argc, argv, argtable);
    if (a_help->count>0)
    {
        cout << "Usage: " << progstr; arg_print_syntax(stdout, argtable, "\n");
        cout << endl; arg_print_glossary(stdout, argtable, "  %-25s %s\n");
        cout << endl << descr; return 1;
    }
    if (nerrs>0) { arg_print_errors(stderr,a_end,(progstr+": "+to_string(__LINE__)+errstr).c_str()); return 1; }


    //Check stdin
    stdi1 = (a_fi->count==0 || strlen(a_fi->filename[0])==0 || strcmp(a_fi->filename[0],"-")==0);
    if (stdi1>0 && isatty(fileno(stdin))) { cerr << progstr+": " << __LINE__ << errstr << "no stdin detected" << endl; return 1; }


    //Check stdout
    if (a_fo->count>0) { stdo1 = (strlen(a_fo->filename[0])==0 || strcmp(a_fo->filename[0],"-")==0); }
    else { stdo1 = (!isatty(fileno(stdout))); }
    wo1 = (stdo1 || a_fo->count>0);


    //Open input
    if (stdi1) { ifs1.copyfmt(cin); ifs1.basic_ios<char>::rdbuf(cin.rdbuf()); } else { ifs1.open(a_fi->filename[0]); }
    if (!ifs1) { cerr << progstr+": " << __LINE__ << errstr << "problem opening input file" << endl; return 1; }


    //Read input header
    if (!read_input_header(ifs1,i1)) { cerr << progstr+": " << __LINE__ << errstr << "problem reading header for input file" << endl; return 1; }
    if ((i1.T==oktypes).sum()==0)
    {
        cerr << progstr+": " << __LINE__ << errstr << "input data type must be in " << "{";
        for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1u]) ? "}" : ","); }
        cerr << endl; return 1;
    }


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


    //Set output header info


    //Open output
    if (wo1)
    {
        if (stdo1) { ofs1.copyfmt(cout); ofs1.basic_ios<char>::rdbuf(cout.rdbuf()); } else { ofs1.open(a_fo->filename[0]); }
        if (!ofs1) { cerr << progstr+": " << __LINE__ << errstr << "problem opening output file 1" << endl; return 1; }
    }


    //Write output header
    if (wo1 && !write_output_header(ofs1,o1)) { cerr << progstr+": " << __LINE__ << errstr << "problem writing header for output file 1" << endl; return 1; }


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
    else if (o1.T==2)
    {
        valarray<int> Yi(o1.N()); valarray<double> Y(o1.N());
        discrete_distribution<uint> distr(&P[0],&P[i1.N()]);
        try { generate_n(begin(Yi),o1.N(),[&distr,&pcg_eng](){return distr(pcg_eng);}); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem during generate" << endl; return 1; }
        for (size_t n=0u; n<o1.N(); ++n) { Y[n] = double(Yi[n]); }
        if (Y.size()!=o1.N()) { cerr << progstr+": " << __LINE__ << errstr << "unexpected output size" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(&Y[0]),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
    }
    else
    {
        cerr << progstr+": " << __LINE__ << errstr << "data type not supported" << endl; return 1;
    }
    

    //Exit
    return ret;
}

