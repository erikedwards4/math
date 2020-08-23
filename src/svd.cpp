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
#include "svd.c"

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
    size_t Kmax, K;


    //Description
    string descr;
    descr += "Linear algebra function.\n";
    descr += "Singular decomposition of matrix X.\n";
    descr += "\n";
    descr += "The outputs U and V are the left and right singular vectors, respectively.\n";
    descr += "The output S is a vector of singular values, such that: X = U*diagmat(S)*V'.\n";
    descr += "\n";
    descr += "Use -k (--K) to give the number of singular vals and vecs to keep [default=all].\n";
    descr += "The K largest singular values and corresponding vectors are returned.\n";
    descr += "The default for K is K_all = min(R,C). \n";
    descr += "\n";
    descr += "List each output filename following a -o opt, in the order U, S, V.\n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ svd X -o U -o S -o V \n";
    descr += "$ svd X -o U -o - -o V > S \n";
    descr += "$ cat X | svd -o U -o S -o V \n";


    //Argtable
    int nerrs;
    struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
    struct arg_int    *a_k = arg_intn("k","K","<uint>",0,1,"num svdencomponents to keep [default=all]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output files (U,S,V)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_fi, a_k, a_fo, a_help, a_end};
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

    //Get K
    Kmax = (i1.R<i1.C) ? i1.R : i1.C;
    if (a_k->count==0) { K = Kmax; }
    else if (a_k->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "K must be positive" << endl; return 1; }
    else { K = size_t(a_k->ival[0]); }
    if (K>Kmax) { cerr << progstr+": " << __LINE__ << errstr << "K must be <= min(R,C)" << endl; return 1; }


    //Checks
    if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }
    if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) must be a matrix" << endl; return 1; }


    //Set output header infos
    o1.F = o2.F = o3.F = i1.F;
    o1.T = o3.T = i1.T;
    o2.T = (i1.isreal()) ? i1.T : i1.T-100u;
    o1.R = i1.R;
    o1.C = K;
    o2.R = K;
    o2.C = 1u;
    o3.R = K;
    o3.C = i1.C;
    o1.S = o2.S = o3.S = i1.S;
    o1.H = o2.H = o3.H = i1.H;


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
        float *X, *U, *S, *Vt;
        try { X = new float[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { U = new float[o1.R*Kmax]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (U)" << endl; return 1; }
        try { S = new float[Kmax]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (S)" << endl; return 1; }
        try { Vt = new float[Kmax*o3.C]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 3 (Vt)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (codee::svd_inplace_s(U,S,Vt,X,i1.R,i1.C,i1.iscolmajor(),K))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(U),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 1 (U)" << endl; return 1; }
        }
        if (wo2)
        {
            try { ofs2.write(reinterpret_cast<char*>(S),o2.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 2 (S)" << endl; return 1; }
        }
        if (wo3)
        {
            try { ofs3.write(reinterpret_cast<char*>(Vt),o3.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 3 (Vt)" << endl; return 1; }
        }
        delete[] X; delete[] U; delete[] S; delete[] Vt;
    }
    else if (i1.T==2)
    {
        double *X, *U, *S, *Vt;
        try { X = new double[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { U = new double[o1.R*Kmax]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (U)" << endl; return 1; }
        try { S = new double[Kmax]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (S)" << endl; return 1; }
        try { Vt = new double[Kmax*o3.C]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 3 (Vt)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (codee::svd_inplace_d(U,S,Vt,X,i1.R,i1.C,i1.iscolmajor(),K))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(U),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 1 (U)" << endl; return 1; }
        }
        if (wo2)
        {
            try { ofs2.write(reinterpret_cast<char*>(S),o2.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 2 (S)" << endl; return 1; }
        }
        if (wo3)
        {
            try { ofs3.write(reinterpret_cast<char*>(Vt),o3.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 3 (Vt)" << endl; return 1; }
        }
        delete[] X; delete[] U; delete[] S; delete[] Vt;
    }
    else if (i1.T==101u)
    {
        float *X, *U, *S, *Vt;
        try { X = new float[2u*i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { U = new float[2u*o1.R*Kmax]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (U)" << endl; return 1; }
        try { S = new float[Kmax]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (S)" << endl; return 1; }
        try { Vt = new float[2u*Kmax*o3.C]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 3 (Vt)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (codee::svd_inplace_c(U,S,Vt,X,i1.R,i1.C,i1.iscolmajor(),K))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(U),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 1 (U)" << endl; return 1; }
        }
        if (wo2)
        {
            try { ofs2.write(reinterpret_cast<char*>(S),o2.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 2 (S)" << endl; return 1; }
        }
        if (wo3)
        {
            try { ofs3.write(reinterpret_cast<char*>(Vt),o3.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 3 (Vt)" << endl; return 1; }
        }
        delete[] X; delete[] U; delete[] S; delete[] Vt;
    }
    else if (i1.T==102u)
    {
        double *X, *U, *S, *Vt;
        try { X = new double[2u*i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { U = new double[2u*o1.R*Kmax]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (U)" << endl; return 1; }
        try { S = new double[Kmax]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (S)" << endl; return 1; }
        try { Vt = new double[2u*Kmax*o3.C]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 3 (Vt)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (codee::svd_inplace_z(U,S,Vt,X,i1.R,i1.C,i1.iscolmajor(),K))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(U),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 1 (U)" << endl; return 1; }
        }
        if (wo2)
        {
            try { ofs2.write(reinterpret_cast<char*>(S),o2.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 2 (S)" << endl; return 1; }
        }
        if (wo3)
        {
            try { ofs3.write(reinterpret_cast<char*>(Vt),o3.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file 3 (Vt)" << endl; return 1; }
        }
        delete[] X; delete[] U; delete[] S; delete[] Vt;
    }
    else
    {
        cerr << progstr+": " << __LINE__ << errstr << "data type not supported" << endl; return 1;
    }
    

    //Exit
    return ret;
}

