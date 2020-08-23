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
#include "split3.c"

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
    const size_t I = 1u, O = 3u;
    ifstream ifs1; ofstream ofs1, ofs2, ofs3;
    int8_t stdi1, stdo1, stdo2, stdo3, wo1, wo2, wo3;
    ioinfo i1, o1, o2, o3;
    size_t dim;


    //Description
    string descr;
    descr += "Splits input X into 3 equal-sized outputs Y1, Y2, Y3.\n";
    descr += "\n";
    descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
    descr += "Use -d0 to work along cols --> Y1, Y2, Y3 have R/3 rows.\n";
    descr += "Use -d1 to work along rows --> Y1, Y2, Y3 have C/3 cols.\n";
    descr += "Use -d2 to work along slices --> Y1, Y2, Y3 have S/3 slices.\n";
    descr += "Use -d3 to work along hyperslices --> Y1, Y2, Y3 have H/3 hyperslices.\n";
    descr += "\n";
    descr += "For dim=0, num rows X must be multiple of 3 (R%3==0).\n";
    descr += "For dim=1, num cols X must be multiple of 3 (C%3==0).\n";
    descr += "For dim=2, num slices X must be multiple of 3 (S%3==0).\n";
    descr += "For dim=3, num hyperslices X must be multiple of 3 (H%3==0).\n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ split3 X -o Y1 -o Y2 -o Y3 \n";
    descr += "$ split3 X -o Y1 -o Y2 > Y3 \n";
    descr += "$ split3 -d1 X -o Y1 -o Y2 -o Y3 \n";
    descr += "$ cat X | split3 -o Y1 -o Y2 -o Y3 \n";


    //Argtable
    int nerrs;
    struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
    struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output files (Y1,Y2,Y3)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_fi, a_d, a_fo, a_help, a_end};
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
    if (a_fo->count>1) { stdo2 = (strlen(a_fo->filename[1])==0 || strcmp(a_fo->filename[1],"-")==0); }
    else { stdo2 = (!isatty(fileno(stdout)) && a_fo->count==1 && stdo1==0); }
    if (a_fo->count>2) { stdo3 = (strlen(a_fo->filename[2])==0 || strcmp(a_fo->filename[2],"-")==0); }
    else { stdo3 = (!isatty(fileno(stdout)) && a_fo->count==2 && stdo1+stdo2==0); }
    if (stdo1+stdo2+stdo3>1) { cerr << progstr+": " << __LINE__ << errstr << "can only use stdout for one output" << endl; return 1; }
    wo1 = (stdo1 || a_fo->count>0); wo2 = (stdo2 || a_fo->count>1); wo3 = (stdo3 || a_fo->count>2);


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

    //Get dim
    if (a_d->count==0) { dim = 0u; }
    else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
    else { dim = size_t(a_d->ival[0]); }
    if (dim>3) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1u,2u,3}" << endl; return 1; }


    //Checks
    if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
    if (dim==0 && i1.R%3) { cerr << progstr+": " << __LINE__ << errstr << "num rows X must be a multiple of 3 for dim=0" << endl; return 1; }
    if (dim==1 && i1.C%3) { cerr << progstr+": " << __LINE__ << errstr << "num cols X must be a multiple of 3 for dim=1" << endl; return 1; }
    if (dim==2 && i1.S%3) { cerr << progstr+": " << __LINE__ << errstr << "num slices X must be a multiple of 3 for dim=2" << endl; return 1; }
    if (dim==3 && i1.H%3) { cerr << progstr+": " << __LINE__ << errstr << "num hyperslices X must be a multiple of 3 for dim=3" << endl; return 1; }


    //Set output header infos
    o1.F = o2.F = o3.F = i1.F;
    o1.T = o2.T = o3.T = i1.T;
    o1.R = o2.R = o3.R = (dim==0) ? i1.R/3 : i1.R;
    o1.C = o2.C = o3.C = (dim==1) ? i1.C/3 : i1.C;
    o1.S = o2.S = o3.S = (dim==2) ? i1.S/3 : i1.S;
    o1.H = o2.H = o3.H = (dim==3) ? i1.H/3 : i1.H;


    //Open outputs
    if (wo1)
    {
        if (stdo1) { ofs1.copyfmt(cout); ofs1.basic_ios<char>::rdbuf(cout.rdbuf()); } else { ofs1.open(a_fo->filename[0]); }
        if (!ofs1) { cerr << progstr+": " << __LINE__ << errstr << "problem opening output file 1" << endl; return 1; }
    }
    if (wo2)
    {
        if (stdo2) { ofs2.copyfmt(cout); ofs2.basic_ios<char>::rdbuf(cout.rdbuf()); } else { ofs2.open(a_fo->filename[1]); }
        if (!ofs2) { cerr << progstr+": " << __LINE__ << errstr << "problem opening output file 2" << endl; return 1; }
    }
    if (wo3)
    {
        if (stdo3) { ofs3.copyfmt(cout); ofs3.basic_ios<char>::rdbuf(cout.rdbuf()); } else { ofs3.open(a_fo->filename[2]); }
        if (!ofs3) { cerr << progstr+": " << __LINE__ << errstr << "problem opening output file 3" << endl; return 1; }
    }


    //Write output headers
    if (wo1 && !write_output_header(ofs1,o1)) { cerr << progstr+": " << __LINE__ << errstr << "problem writing header for output file 1" << endl; return 1; }
    if (wo2 && !write_output_header(ofs2,o2)) { cerr << progstr+": " << __LINE__ << errstr << "problem writing header for output file 2" << endl; return 1; }
    if (wo3 && !write_output_header(ofs3,o3)) { cerr << progstr+": " << __LINE__ << errstr << "problem writing header for output file 3" << endl; return 1; }


    //Other prep


    //Process
    if (i1.T==1u)
    {
        float *X, *Y1, *Y2, *Y3;
        try { X = new float[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y1 = new float[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (Y1)" << endl; return 1; }
        try { Y2 = new float[o2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (Y2)" << endl; return 1; }
        try { Y3 = new float[o3.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 3 (Y3)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (codee::split3_s(Y1,Y2,Y3,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y1),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 1 (Y1)" << endl; return 1; }
        }
        if (wo2)
        {
            try { ofs2.write(reinterpret_cast<char*>(Y2),o2.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 2 (Y2)" << endl; return 1; }
        }
        if (wo3)
        {
            try { ofs3.write(reinterpret_cast<char*>(Y3),o3.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 3 (Y3)" << endl; return 1; }
        }
        delete[] X; delete[] Y1; delete[] Y2; delete[] Y3;
    }
    else if (i1.T==2)
    {
        double *X, *Y1, *Y2, *Y3;
        try { X = new double[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y1 = new double[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (Y1)" << endl; return 1; }
        try { Y2 = new double[o2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (Y2)" << endl; return 1; }
        try { Y3 = new double[o3.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 3 (Y3)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (codee::split3_d(Y1,Y2,Y3,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y1),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 1 (Y1)" << endl; return 1; }
        }
        if (wo2)
        {
            try { ofs2.write(reinterpret_cast<char*>(Y2),o2.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 2 (Y2)" << endl; return 1; }
        }
        if (wo3)
        {
            try { ofs3.write(reinterpret_cast<char*>(Y3),o3.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 3 (Y3)" << endl; return 1; }
        }
        delete[] X; delete[] Y1; delete[] Y2; delete[] Y3;
    }
    else if (i1.T==101u)
    {
        float *X, *Y1, *Y2, *Y3;
        try { X = new float[2u*i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y1 = new float[2u*o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (Y1)" << endl; return 1; }
        try { Y2 = new float[2u*o2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (Y2)" << endl; return 1; }
        try { Y3 = new float[2u*o3.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 3 (Y3)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (codee::split3_c(Y1,Y2,Y3,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y1),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 1 (Y1)" << endl; return 1; }
        }
        if (wo2)
        {
            try { ofs2.write(reinterpret_cast<char*>(Y2),o2.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 2 (Y2)" << endl; return 1; }
        }
        if (wo3)
        {
            try { ofs3.write(reinterpret_cast<char*>(Y3),o3.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 3 (Y3)" << endl; return 1; }
        }
        delete[] X; delete[] Y1; delete[] Y2; delete[] Y3;
    }
    else if (i1.T==102u)
    {
        double *X, *Y1, *Y2, *Y3;
        try { X = new double[2u*i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y1 = new double[2u*o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (Y1)" << endl; return 1; }
        try { Y2 = new double[2u*o2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (Y2)" << endl; return 1; }
        try { Y3 = new double[2u*o3.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 3 (Y3)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (codee::split3_z(Y1,Y2,Y3,X,i1.R,i1.C,i1.S,i1.H,i1.iscolmajor(),dim))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y1),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 1 (Y1)" << endl; return 1; }
        }
        if (wo2)
        {
            try { ofs2.write(reinterpret_cast<char*>(Y2),o2.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 2 (Y2)" << endl; return 1; }
        }
        if (wo3)
        {
            try { ofs3.write(reinterpret_cast<char*>(Y3),o3.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 3 (Y3)" << endl; return 1; }
        }
        delete[] X; delete[] Y1; delete[] Y2; delete[] Y3;
    }
    else
    {
        cerr << progstr+": " << __LINE__ << errstr << "data type not supported" << endl; return 1;
    }
    

    //Exit
    return ret;
}

