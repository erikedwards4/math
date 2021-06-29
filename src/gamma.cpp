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
    const valarray<size_t> oktypes = {1u,2u,101u,102u};
    const size_t O = 1u;
    ofstream ofs1;
    int8_t stdo1, wo1;
    ioinfo o1;
    double a, b;


    //Description
    string descr;
    descr += "Random function: gamma_distribution.\n";
    descr += "Outputs tensor of floats from a gamma distribution, \n";
    descr += "with shape param alpha>0 and scale param beta>0.\n";
    descr += "\n";
    descr += "For complex output, real/imag parts are set separately.\n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ gamma -a2.5 -b2 -r2 -c3 -o Y \n";
    descr += "$ gamma -a2.5 -b2 -r2 -c3 -t1 > Y \n";
    descr += "$ gamma -a2.5 -b2 -r2 -c3 -t102 > Y \n";


    //Argtable
    int nerrs;
    struct arg_dbl    *a_a = arg_dbln("a","alpha","<dbl>",0,1,"shape parameter [default=1]");
    struct arg_dbl    *a_b = arg_dbln("b","beta","<dbl>",0,1,"scale parameter [default=1]");
    struct arg_int   *a_nr = arg_intn("r","R","<uint>",0,1,"num rows in output [default=1]");
    struct arg_int   *a_nc = arg_intn("c","C","<uint>",0,1,"num cols in output [default=1]");
    struct arg_int   *a_ns = arg_intn("s","S","<uint>",0,1,"num slices in output [default=1]");
    struct arg_int   *a_nh = arg_intn("y","H","<uint>",0,1,"num hyperslices in the output [default=1]");
    struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=1]");
    struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=147]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_a, a_b, a_nr, a_nc, a_ns, a_nh, a_otyp, a_ofmt, a_fo, a_help, a_end};
    if (arg_nullcheck(argtable)!=0) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating argtable" << endl; return 1; }
    nerrs = arg_parse(argc, argv, argtable);
    if (a_help->count>0)
    {
        cout << "Usage: " << progstr; arg_print_syntax(stdout, argtable, "\n");
        cout << endl; arg_print_glossary(stdout, argtable, "  %-25s %s\n");
        cout << endl << descr; return 1;
    }
    if (nerrs>0) { arg_print_errors(stderr,a_end,(progstr+": "+to_string(__LINE__)+errstr).c_str()); return 1; }


    //Check stdout
    if (a_fo->count>0) { stdo1 = (strlen(a_fo->filename[0])==0 || strcmp(a_fo->filename[0],"-")==0); }
    else { stdo1 = (!isatty(fileno(stdout))); }
    wo1 = (stdo1 || a_fo->count>0);


    //Get options

    //Get o1.F
    if (a_ofmt->count==0) { o1.F = 147u; }
    else if (a_ofmt->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be nonnegative" << endl; return 1; }
    else if (a_ofmt->ival[0]>255) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be < 256" << endl; return 1; }
    else { o1.F = size_t(a_ofmt->ival[0]); }

    //Get o1.T
    if (a_otyp->count==0) { o1.T = 1u; }
    else if (a_otyp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "output data type must be positive" << endl; return 1; }
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

    //Get a
    a = (a_a->count>0) ? a_a->dval[0] : 1.0;
    if (a<0.0) { cerr << progstr+": " << __LINE__ << errstr << "a (alpha) must be nonnegative " << endl; return 1; }

    //Get b (b<0 gives negative output)
    b = (a_b->count>0) ? a_b->dval[0] : 1.0;


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
    

    //Process
    if (o1.T==1u)
    {
        valarray<float> Y(o1.N());
        gamma_distribution<float> distr(a,b);
        try { generate_n(begin(Y),o1.N(),[&distr,&pcg_eng](){return distr(pcg_eng);}); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem during generate" << endl; return 1; }
        if (Y.size()!=o1.N()) { cerr << progstr+": " << __LINE__ << errstr << "unexpected output size" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(&Y[0]),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
    }
    else if (o1.T==2)
    {
        valarray<double> Y(o1.N());
        gamma_distribution<double> distr(a,b);
        try { generate_n(begin(Y),o1.N(),[&distr,&pcg_eng](){return distr(pcg_eng);}); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem during generate" << endl; return 1; }
        if (Y.size()!=o1.N()) { cerr << progstr+": " << __LINE__ << errstr << "unexpected output size" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(&Y[0]),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
    }
    else if (o1.T==101u)
    {
        valarray<complex<float>> Y(o1.N());
        gamma_distribution<float> distr(a,b);
        try { generate_n(begin(Y),o1.N(),[&distr,&pcg_eng](){complex<float> y(distr(pcg_eng),distr(pcg_eng)); return y;}); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem during generate" << endl; return 1; }
        if (Y.size()!=o1.N()) { cerr << progstr+": " << __LINE__ << errstr << "unexpected output size" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(&Y[0]),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
    }
    else if (o1.T==102u)
    {
        valarray<complex<double>> Y(o1.N());
        gamma_distribution<double> distr(a,b);
        try { generate_n(begin(Y),o1.N(),[&distr,&pcg_eng](){complex<double> y(distr(pcg_eng),distr(pcg_eng)); return y;}); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem during generate" << endl; return 1; }
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

