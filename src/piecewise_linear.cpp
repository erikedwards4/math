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
#include <cfloat>
#include <random>
#include <complex>

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
    const size_t I = 2, O = 1;
    ifstream ifs1, ifs2; ofstream ofs1;
    int8_t stdi1, stdi2, stdo1, wo1;
    ioinfo i1, i2, o1;
    random_device rd;    //random device to seed Mersenne twister engine
    mt19937 mt_eng(rd()); 


    //Description
    string descr;
    descr += "Random function: piecewise_linear_distribution.\n";
    descr += "Outputs tensor of floats from a piecewise_linear distribution, \n";
    descr += "with bin edges B and probabilities P.\n";
    descr += "\n";
    descr += "B must be a row or column vector with length N,\n";
    descr += "with a monotonically-increasing set of floats.\n";
    descr += "\n";
    descr += "P must be a row or column vector with length N.\n";
    descr += "The elements of P must be positive, but do not have to be probabilities.\n";
    descr += "They are treated as weights and normalized to sum to 1.\n";
    descr += "\n";
    descr += "The probability at bin edge B[n] will be exactly P[n],\n";
    descr += "with linear change of probability in between.\n";
    descr += "\n";
    descr += "B and P must have the same data type, and Y will have this same data type.\n";
    descr += "\n";
    descr += "The output (Y) has floats in [B[0] B[N-1]]. \n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ piecewise_linear B P -r2 -c3 -o Y \n";
    descr += "$ piecewise_linear B P -r2 -c3 > Y \n";
    descr += "$ cat P | piecewise_linear B -r2 -c3 > Y \n";
    descr += "$ cat B | piecewise_linear - P -r2 -c3 > Y \n";


    //Argtable
    int nerrs;
    struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (B,P)");
    struct arg_int   *a_nr = arg_intn("r","R","<uint>",0,1,"num rows in output [default=1]");
    struct arg_int   *a_nc = arg_intn("c","C","<uint>",0,1,"num cols in output [default=1]");
    struct arg_int   *a_ns = arg_intn("s","S","<uint>",0,1,"num slices in output [default=1]");
    struct arg_int   *a_nh = arg_intn("y","H","<uint>",0,1,"num hyperslices in the output [default=1]");
    struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default from B]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_fi, a_nr, a_nc, a_ns, a_nh, a_ofmt, a_fo, a_help, a_end};
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
    stdi2 = (a_fi->count<=1 || strlen(a_fi->filename[1])==0 || strcmp(a_fi->filename[1],"-")==0);
    if (stdi1+stdi2>1) { cerr << progstr+": " << __LINE__ << errstr << "can only use stdin for one input" << endl; return 1; }
    if (stdi1+stdi2>0 && isatty(fileno(stdin))) { cerr << progstr+": " << __LINE__ << errstr << "no stdin detected" << endl; return 1; }


    //Check stdout
    if (a_fo->count>0) { stdo1 = (strlen(a_fo->filename[0])==0 || strcmp(a_fo->filename[0],"-")==0); }
    else { stdo1 = (!isatty(fileno(stdout))); }
    wo1 = (stdo1 || a_fo->count>0);


    //Open inputs
    if (stdi1) { ifs1.copyfmt(cin); ifs1.basic_ios<char>::rdbuf(cin.rdbuf()); } else { ifs1.open(a_fi->filename[0]); }
    if (!ifs1) { cerr << progstr+": " << __LINE__ << errstr << "problem opening input file 1" << endl; return 1; }
    if (stdi2) { ifs2.copyfmt(cin); ifs2.basic_ios<char>::rdbuf(cin.rdbuf()); } else { ifs2.open(a_fi->filename[1]); }
    if (!ifs2) { cerr << progstr+": " << __LINE__ << errstr << "problem opening input file 2" << endl; return 1; }


    //Read input headers
    if (!read_input_header(ifs1,i1)) { cerr << progstr+": " << __LINE__ << errstr << "problem reading header for input file 1" << endl; return 1; }
    if (!read_input_header(ifs2,i2)) { cerr << progstr+": " << __LINE__ << errstr << "problem reading header for input file 2" << endl; return 1; }
    if ((i1.T==oktypes).sum()==0 || (i2.T==oktypes).sum()==0)
    {
        cerr << progstr+": " << __LINE__ << errstr << "input data type must be in " << "{";
        for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1]) ? "}" : ","); }
        cerr << endl; return 1;
    }


    //Get options

    //Get o1.F
    if (a_ofmt->count==0) { o1.F = i1.F; }
    else if (a_ofmt->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be nonnegative" << endl; return 1; }
    else if (a_ofmt->ival[0]>255) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be < 256" << endl; return 1; }
    else { o1.F = uint8_t(a_ofmt->ival[0]); }

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


    //Set output header info
    o1.T = i1.T;


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
        valarray<float> B(i1.N()), P(i2.N()), Y(o1.N());
        try { ifs1.read(reinterpret_cast<char*>(&B[0]),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading data for input file 1 (B)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(&P[0]),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading data for input file 2 (P)" << endl; return 1; }
        if ((P<0.0f).sum()>0) { cerr << progstr+": " << __LINE__ << errstr << "elements of P must be nonnegative" << endl; return 1; }
        piecewise_linear_distribution<float> distr(&B[0],&B[i1.N()],&P[0]);
        try { generate_n(begin(Y),o1.N(),[&distr,&mt_eng](){return distr(mt_eng);}); }
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
        valarray<double> B(i1.N()), P(i2.N()), Y(o1.N());
        try { ifs1.read(reinterpret_cast<char*>(&B[0]),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading data for input file 1 (B)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(&P[0]),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading data for input file 2 (P)" << endl; return 1; }
        if ((P<0.0).sum()>0) { cerr << progstr+": " << __LINE__ << errstr << "elements of P must be nonnegative" << endl; return 1; }
        piecewise_linear_distribution<double> distr(&B[0],&B[i1.N()],&P[0]);
        try { generate_n(begin(Y),o1.N(),[&distr,&mt_eng](){return distr(mt_eng);}); }
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

