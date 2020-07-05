//@author Erik Edwards
//@date 2019-2020


#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <cstring>
#include <valarray>
#include <unordered_map>
#include <argtable2.h>
#include "/home/erik/codee/util/cmli.hpp"
#include "primes.c"

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
    const valarray<uint8_t> oktypes = {1,2};
    const size_t O = 1;
    ofstream ofs1;
    int8_t stdo1, wo1;
    ioinfo o1;
    size_t P, cnt;
    size_t *X;


    //Description
    string descr;
    descr += "Generate function (0 inputs, 1 output).\n";
    descr += "Output Y is a row vector filled with prime numbers from 2 up to P.\n";
    descr += "\n";
    descr += "Use -p (--P) to specify the largest integer below which to return primes [default=2].\n";
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
    int nerrs;
    struct arg_int    *a_p = arg_intn("p","P","<uint>",0,1,"return primes <= P [default=2]");
    struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=1]");
    struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=147]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_p, a_otyp, a_ofmt, a_fo, a_help, a_end};
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

    //Get P
    if (a_p->count==0) { P = 2u; }
    else if (a_p->ival[0]<2) { cerr << progstr+": " << __LINE__ << errstr << "P must be greater than 1" << endl; return 1; }
    else { P = size_t(a_p->ival[0]); }

    //Get o1.F
    if (a_ofmt->count==0) { o1.F = 147; }
    else if (a_ofmt->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be nonnegative" << endl; return 1; }
    else if (a_ofmt->ival[0]>255) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be < 256" << endl; return 1; }
    else { o1.F = uint8_t(a_ofmt->ival[0]); }

    //Get o1.T
    if (a_otyp->count==0) { o1.T = 1; }
    else if (a_otyp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "data type must be positive int" << endl; return 1; }
    else { o1.T = uint8_t(a_otyp->ival[0]); }
    if ((o1.T==oktypes).sum()==0)
    {
        cerr << progstr+": " << __LINE__ << errstr << "output data type must be in " << "{";
        for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1]) ? "}" : ","); }
        cerr << endl; return 1;
    }


    //Set output header info
    try { X = new size_t[(P-1)/2]; }
    catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for integer file (X)" << endl; return 1; }
    if (codee::primes_i(X,&cnt,(P-1)/2))
    { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
    o1.C = size_t(cnt);
    o1.R = o1.S = o1.H = 1u;


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
        float *Y;
        try { Y = new float[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        for (size_t c=0u; c<o1.C; c++) { Y[c] = float(X[c]); }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] Y;
    }
    else if (o1.T==2)
    {
        double *Y;
        try { Y = new double[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        for (size_t c=0u; c<o1.C; c++) { Y[c] = double(X[c]); }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] Y;
    }
    else
    {
        cerr << progstr+": " << __LINE__ << errstr << "data type not supported" << endl; return 1;
    }
    

    //Exit
    return ret;
}

