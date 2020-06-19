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
#include "hypot.c"

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
    const size_t I = 2, O = 1;
    ifstream ifs1, ifs2; ofstream ofs1;
    int8_t stdi1, stdi2, stdo1, wo1;
    ioinfo i1, i2, o1;


    //Description
    string descr;
    descr += "Elementwise function for 2 inputs with broadcasting.\n";
    descr += "Gets hypotenuse for each corresponding element of X1, X2:\n";
    descr += "y = sqrt(x1^2 + x2^2) \n";
    descr += "\n";
    descr += "For complex input, output is real-valued: \n";
    descr += "y = sqrt(|x1|^2 + |x2|^2) \n";
    descr += "\n";
    descr += "X1 and X2 must have the same size or broadcast-compatible sizes.\n";
    descr += "Output (Y) has size max(R1,R2) x max(C1,C2).\n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ hypot X1 X2 -o Y \n";
    descr += "$ hypot X1 X2 > Y \n";
    descr += "$ cat X1 | hypot - X2 > Y \n";


    //Argtable
    int nerrs;
    struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2)");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_fi, a_fo, a_help, a_end};
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


    //Checks
    if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
    if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }
    if (!i2.isvec() && i2.iscolmajor()!=i1.iscolmajor()) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same row/col major format" << endl; return 1; }
    if (i2.T!=i1.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
    if (!bcast_compat(i1,i2)) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have same size or broadcast-compatible sizes" << endl; return 1; }


    //Set output header info
    o1.F = i1.F; o1.T = (i1.isreal()) ? i1.T : i1.T-100;
    o1.R = (i1.R>i2.R) ? i1.R : i2.R;
    o1.C = (i1.C>i2.C) ? i1.C : i2.C;
    o1.S = (i1.S>i2.S) ? i1.S : i2.S;
    o1.H = (i1.H>i2.H) ? i1.H : i2.H;


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
    if (i1.T==1)
    {
        float *X1, *X2;
        try { X1 = new float[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
        try { X2 = new float[i2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
        if (i1.N()==o1.N())
        {
            if (codee::hypot_inplace_s(X1,X2,int(i1.R),int(i1.C),int(i1.S),int(i1.H),int(i2.R),int(i2.C),int(i2.S),int(i2.H),o1.iscolmajor()))
            { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
            if (wo1)
            {
                try { ofs1.write(reinterpret_cast<char*>(X1),o1.nbytes()); }
                catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
            }
        }
        else
        {
            float *Y;
            try { Y = new float[o1.N()]; }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
            if (codee::hypot_s(Y,X1,X2,int(i1.R),int(i1.C),int(i1.S),int(i1.H),int(i2.R),int(i2.C),int(i2.S),int(i2.H),o1.iscolmajor()))
            { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
            if (wo1)
            {
                try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
                catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
            }
            delete[] Y;
        }
        delete[] X1; delete[] X2;
    }
    else if (i1.T==2)
    {
        double *X1, *X2;
        try { X1 = new double[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
        try { X2 = new double[i2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
        if (i1.N()==o1.N())
        {
            if (codee::hypot_inplace_d(X1,X2,int(i1.R),int(i1.C),int(i1.S),int(i1.H),int(i2.R),int(i2.C),int(i2.S),int(i2.H),o1.iscolmajor()))
            { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
            if (wo1)
            {
                try { ofs1.write(reinterpret_cast<char*>(X1),o1.nbytes()); }
                catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
            }
        }
        else
        {
            double *Y;
            try { Y = new double[o1.N()]; }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
            if (codee::hypot_d(Y,X1,X2,int(i1.R),int(i1.C),int(i1.S),int(i1.H),int(i2.R),int(i2.C),int(i2.S),int(i2.H),o1.iscolmajor()))
            { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
            if (wo1)
            {
                try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
                catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
            }
            delete[] Y;
        }
        delete[] X1; delete[] X2;
    }
    else if (i1.T==101)
    {
        float *X1, *X2;
        try { X1 = new float[2u*i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
        try { X2 = new float[2u*i2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
        if (i1.N()==o1.N())
        {
            if (codee::hypot_inplace_c(X1,X2,int(i1.R),int(i1.C),int(i1.S),int(i1.H),int(i2.R),int(i2.C),int(i2.S),int(i2.H),o1.iscolmajor()))
            { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
            if (wo1)
            {
                try { ofs1.write(reinterpret_cast<char*>(X1),o1.nbytes()); }
                catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
            }
        }
        else
        {
            float *Y;
            try { Y = new float[o1.N()]; }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
            if (codee::hypot_c(Y,X1,X2,int(i1.R),int(i1.C),int(i1.S),int(i1.H),int(i2.R),int(i2.C),int(i2.S),int(i2.H),o1.iscolmajor()))
            { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
            if (wo1)
            {
                try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
                catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
            }
            delete[] Y;
        }
        delete[] X1; delete[] X2;
    }
    else if (i1.T==102)
    {
        double *X1, *X2;
        try { X1 = new double[2u*i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
        try { X2 = new double[2u*i2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
        if (i1.N()==o1.N())
        {
            if (codee::hypot_inplace_z(X1,X2,int(i1.R),int(i1.C),int(i1.S),int(i1.H),int(i2.R),int(i2.C),int(i2.S),int(i2.H),o1.iscolmajor()))
            { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
            if (wo1)
            {
                try { ofs1.write(reinterpret_cast<char*>(X1),o1.nbytes()); }
                catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
            }
        }
        else
        {
            double *Y;
            try { Y = new double[o1.N()]; }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
            if (codee::hypot_z(Y,X1,X2,int(i1.R),int(i1.C),int(i1.S),int(i1.H),int(i2.R),int(i2.C),int(i2.S),int(i2.H),o1.iscolmajor()))
            { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
            if (wo1)
            {
                try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
                catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
            }
            delete[] Y;
        }
        delete[] X1; delete[] X2;
    }
    else
    {
        cerr << progstr+": " << __LINE__ << errstr << "data type not supported" << endl; return 1;
    }
    

    //Exit
    return ret;
}

