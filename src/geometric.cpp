//@author Erik Edwards
//@date 2019-2020


#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <cstring>
#include <valarray>
#include <complex>
#include <unordered_map>
#include <argtable2.h>
#include "/home/erik/codee/util/cmli.hpp"
#include <random>

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
    const valarray<uint8_t> oktypes = {1,2,101,102};
    const size_t O = 1;
    ofstream ofs1;
    int8_t stdo1, wo1;
    ioinfo o1;
    double p;
    random_device rd;  //random device to seed Mersenne twister engine
    mt19937 mt_eng(rd());


    //Description
    string descr;
    descr += "Random function: geometric_distribution.\n";
    descr += "Outputs RxC matrix (or RxCxS cube) of ints from a geometric distribution\n";
    descr += "with parameter p (probability of success).\n";
    descr += "\n";
    descr += "Output is nonnegative ints, but any data type can be used.\n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ geometric -p0.2 -r2 -c3 -o Y \n";
    descr += "$ geometric -p0.2 -r2 -c3 > Y \n";
    descr += "$ geometric -p0.2 -r2 -c3 -t2 > Y \n";


    //Argtable
    int nerrs;
    struct arg_dbl    *a_p = arg_dbln("p","p","<dbl>",0,1,"p probability parameter [default=0.5]");
    struct arg_int   *a_nr = arg_intn("r","R","<uint>",0,1,"num rows in output [default=1]");
    struct arg_int   *a_nc = arg_intn("c","C","<uint>",0,1,"num cols in output [default=1]");
    struct arg_int   *a_ns = arg_intn("s","S","<uint>",0,1,"num slices in output [default=1]");
    struct arg_int   *a_nh = arg_intn("y","H","<uint>",0,1,"num hyperslices in the output [default=1]");
    struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=1]");
    struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=147]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_p, a_nr, a_nc, a_ns, a_nh, a_otyp, a_ofmt, a_fo, a_help, a_end};
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


    //Process
    if (o1.T==1)
    {
        valarray<float> Y(o1.N());
        geometric_distribution<uint> distr(p);
        try { generate_n(begin(Y),o1.N(),[&distr,&mt_eng](){return distr(mt_eng);}); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem during generate" << endl; return 1; }
        if (Y.size()!=o1.N()) { cerr << progstr+": " << __LINE__ << errstr << "unexpected output size" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(&Y[0]),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing data for output" << endl; return 1; }
        }
    }
    else if (o1.T==2)
    {
        valarray<double> Y(o1.N());
        geometric_distribution<uint> distr(p);
        try { generate_n(begin(Y),o1.N(),[&distr,&mt_eng](){return distr(mt_eng);}); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem during generate" << endl; return 1; }
        if (Y.size()!=o1.N()) { cerr << progstr+": " << __LINE__ << errstr << "unexpected output size" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(&Y[0]),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing data for output" << endl; return 1; }
        }
    }
    else if (o1.T==101)
    {
        valarray<complex<float>> Y(o1.N());
        geometric_distribution<uint> distr(p);
        try { generate_n(begin(Y),o1.N(),[&distr,&mt_eng](){complex<float> y(distr(mt_eng),distr(mt_eng)); return y;}); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem during generate" << endl; return 1; }
        if (Y.size()!=o1.N()) { cerr << progstr+": " << __LINE__ << errstr << "unexpected output size" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(&Y[0]),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing data for output" << endl; return 1; }
        }
    }
    else if (o1.T==102)
    {
        valarray<complex<double>> Y(o1.N());
        geometric_distribution<uint> distr(p);
        try { generate_n(begin(Y),o1.N(),[&distr,&mt_eng](){complex<double> y(distr(mt_eng),distr(mt_eng)); return y;}); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem during generate" << endl; return 1; }
        if (Y.size()!=o1.N()) { cerr << progstr+": " << __LINE__ << errstr << "unexpected output size" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(&Y[0]),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing data for output" << endl; return 1; }
        }
    }
    else
    {
        cerr << progstr+": " << __LINE__ << errstr << "data type not supported" << endl; return 1;
    }
    

    //Exit
    return ret;
}

