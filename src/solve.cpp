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
#include "solve.c"

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
    const size_t I = 2u, O = 1u;
    ifstream ifs1, ifs2; ofstream ofs1;
    int8_t stdi1, stdi2, stdo1, wo1;
    ioinfo i1, i2, o1;
    char tr;


    //Description
    string descr;
    descr += "Linear algebra function for 2 inputs.\n";
    descr += "Least-squares solution of linear system: A*X = B.\n";
    descr += "The system can be over- or under-determined. \n";
    descr += "\n";
    descr += "This is similar to the \\ (backslash) operator in Octave,\n";
    descr += "with X = A\\B. \n";
    descr += "\n";
    descr += "Input 1 (A) is usually a matrix.\n";
    descr += "Input 2 (B) can be vector or matrix.\n";
    descr += "The output (X) is solution to minimize norm2(B-A*X).\n";
    descr += "\n";
    descr += "Include -t (--transpose) to use A' (Hermitian transpose of A).\n";
    descr += "In this case, minimize: norm2(B-A'*X). \n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ solve A B -o X \n";
    descr += "$ solve A B > X \n";
    descr += "$ cat A | solve -t - B > X \n";


    //Argtable
    int nerrs;
    struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (A,B)");
    struct arg_lit   *a_tr = arg_litn("t","transpose",0,1,"transpose A");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (X)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_fi, a_tr, a_fo, a_help, a_end};
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
        for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1u]) ? "}" : ","); }
        cerr << endl; return 1;
    }


    //Get options

    //Get tr
    tr = (a_tr->count>0);


    //Checks
    if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
    if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (A) found to be empty" << endl; return 1; }
    if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (B) found to be empty" << endl; return 1; }
    if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (A) must be a matrix" << endl; return 1; }
    if (!i2.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (B) must be a vector or matrix" << endl; return 1; }
    if (!major_compat(i1,i2)) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same row/col major format unless vectors" << endl; return 1; }
    if (tr && i1.C!=i2.R) { cerr << progstr+": " << __LINE__ << errstr << "C1 (nrows A) must equal R2 (nrows B) if transpose of A" << endl; return 1; }
    if (!tr && i1.R!=i2.R) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same num rows if no transpose of A" << endl; return 1; }
    if (tr && i1.T>100) { cerr << progstr+": " << __LINE__ << errstr << "transpose option not working for complex data types" << endl; return 1; }


    //Set output header info
    o1.F = i1.F; o1.T = i1.T;
    o1.R = (tr) ? i1.R : i1.C;
    o1.C = i2.C;
    o1.S = i1.S; o1.H = i1.H;


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
    if (i1.T==1u)
    {
        float *A, *B;
        try { A = new float[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (A)" << endl; return 1; }
        try { B = new float[i2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (B)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(A),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (A)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(B),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (B)" << endl; return 1; }
        if (codee::solve_inplace_s(A,B,i1.R,i1.C,i2.C,o1.iscolmajor(),tr))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(B),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (X)" << endl; return 1; }
        }
        delete[] A; delete[] B;
    }
    else if (i1.T==2)
    {
        double *A, *B;
        try { A = new double[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (A)" << endl; return 1; }
        try { B = new double[i2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (B)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(A),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (A)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(B),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (B)" << endl; return 1; }
        if (codee::solve_inplace_d(A,B,i1.R,i1.C,i2.C,o1.iscolmajor(),tr))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(B),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (X)" << endl; return 1; }
        }
        delete[] A; delete[] B;
    }
    else if (i1.T==101u)
    {
        float *A, *B;
        try { A = new float[2u*i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (A)" << endl; return 1; }
        try { B = new float[2u*i2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (B)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(A),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (A)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(B),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (B)" << endl; return 1; }
        if (codee::solve_inplace_c(A,B,i1.R,i1.C,i2.C,o1.iscolmajor(),tr))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(B),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (X)" << endl; return 1; }
        }
        delete[] A; delete[] B;
    }
    else if (i1.T==102u)
    {
        double *A, *B;
        try { A = new double[2u*i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (A)" << endl; return 1; }
        try { B = new double[2u*i2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (B)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(A),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (A)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(B),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (B)" << endl; return 1; }
        if (codee::solve_inplace_z(A,B,i1.R,i1.C,i2.C,o1.iscolmajor(),tr))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(B),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (X)" << endl; return 1; }
        }
        delete[] A; delete[] B;
    }
    else
    {
        cerr << progstr+": " << __LINE__ << errstr << "data type not supported" << endl; return 1;
    }
    

    //Exit
    return ret;
}

