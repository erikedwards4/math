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
#include "dist_cos.c"

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
    const size_t I = 2u, O = 1u;
    ifstream ifs1, ifs2; ofstream ofs1;
    int8_t stdi1, stdi2, stdo1, wo1;
    ioinfo i1, i2, o1;
    size_t dim;


    //Description
    string descr;
    descr += "Vecs2scalar function for 2 inputs with broadcasting.\n";
    descr += "Gets cosine distance for each pair of vectors in X1, X2 along dim.\n";
    descr += "Thus, output Y has length 1 along dim.\n";
    descr += "\n";
    descr += "Use -d (--dim) to give the dimension (axis) [default=0].\n";
    descr += "Use -d0 to work along cols (column vecs).\n";
    descr += "Use -d1 to work along rows (row vecs).\n";
    descr += "Use -d2 to work along slices (tube vecs).\n";
    descr += "Use -d3 to work along hyperslices (hypertube vecs).\n";
    descr += "\n";
    descr += "X1 and X2 must have the same size and data type, or\n";
    descr += "either X1 or X2 can be a single vector of appropriate length.\n";
    descr += "In that case, the single vector will be broadcast to the other input.\n";
    descr += "\n";
    descr += "This is equal to 1 minus the cosine similarity, so is in [0 2].\n";
    descr += "Although it is not a proper distance metric, it is often used.\n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ dist_cos X1 X2 -o Y \n";
    descr += "$ dist_cos -d1 X1 X2 > Y \n";
    descr += "$ cat X1 | dist_cos -d2 - X2 > Y \n";


    //Argtable
    int nerrs;
    struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input files (X1,X2)");
    struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension [default=0]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
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

    //Get dim
    if (a_d->count==0) { dim = 0u; }
    else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
    else { dim = size_t(a_d->ival[0]); }
    if (dim>3) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1u,2u,3}" << endl; return 1; }


    //Checks
    if (i1.T!=i2.T) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same data type" << endl; return 1; }
    if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X1) found to be empty" << endl; return 1; }
    if (i2.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 2 (X2) found to be empty" << endl; return 1; }
    if (!bcast_compat(i1,i2)) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have same size or broadcast-compatible sizes" << endl; return 1; }
    if (!major_compat(i1,i2)) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have the same row/col major format unless vectors" << endl; return 1; }
    size_t L1 = (dim==0) ? i1.R : (dim==1) ? i1.C : (dim==2) ? i1.S : i1.H;
    size_t L2 = (dim==0) ? i2.R : (dim==1) ? i2.C : (dim==2) ? i2.S : i2.H;
    if (L1!=L2) { cerr << progstr+": " << __LINE__ << errstr << "inputs must have same length along dim" << endl; return 1; }


    //Set output header info
    o1.F = i1.F; o1.T = i1.T;
    o1.R = (dim==0) ? 1 : (i1.R>i2.R) ? i1.R : i2.R;
    o1.C = (dim==1) ? 1 : (i1.C>i2.C) ? i1.C : i2.C;
    o1.S = (dim==2) ? 1 : (i1.S>i2.S) ? i1.S : i2.S;
    o1.H = (dim==3) ? 1 : (i1.H>i2.H) ? i1.H : i2.H;


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
        float *X1, *X2, *Y;
        try { X1 = new float[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
        try { X2 = new float[i2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
        try { Y = new float[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
        if (codee::dist_cos_s(Y,X1,X2,i1.R,i1.C,i1.S,i1.H,i2.R,i2.C,i2.S,i2.H,o1.iscolmajor(),dim))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] X1; delete[] X2; delete[] Y;
    }
    else if (i1.T==2)
    {
        double *X1, *X2, *Y;
        try { X1 = new double[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 1 (X1)" << endl; return 1; }
        try { X2 = new double[i2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file 2 (X2)" << endl; return 1; }
        try { Y = new double[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X1),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 1 (X1)" << endl; return 1; }
        try { ifs2.read(reinterpret_cast<char*>(X2),i2.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file 2 (X2)" << endl; return 1; }
        if (codee::dist_cos_d(Y,X1,X2,i1.R,i1.C,i1.S,i1.H,i2.R,i2.C,i2.S,i2.H,o1.iscolmajor(),dim))
        { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] X1; delete[] X2; delete[] Y;
    }
    else
    {
        cerr << progstr+": " << __LINE__ << errstr << "data type not supported" << endl; return 1;
    }
    

    //Exit
    return ret;
}

