//Gets kth diagonal of input X as column vector Y

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int diag_s (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int k);
int diag_d (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int k);
int diag_c (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int k);
int diag_z (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int k);


int diag_s (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int k)
{
    const size_t K = (iscolmajor) ? R+1u : C+1u;

    if (k>=0 && (int)C>k)
    {
        const size_t L = (C-(size_t)k<R) ? C-(size_t)k : R;
        const int S = (iscolmajor) ? k*(int)R : k;
        
        X += S;
        for (size_t l=0u; l<L; ++l, X+=K, ++Y) { *Y = *X; }
    }
    else if (k<=0 && (int)R>-k)
    {
        const size_t L = ((int)R+k<(int)C) ? (size_t)((int)R+k) : C;
        const int S = (iscolmajor) ? -k : -k*(int)C;
        
        X += S;
        for (size_t l=0u; l<L; ++l, X+=K, ++Y) { *Y = *X; }
    }
    else
    {
        fprintf(stderr,"error in diag_s: k out of range [1-R C-1]\n"); return 1;
    }

    // if (iscolmajor)
    // {
    //     if (k>0) { cblas_scopy((int)L,&X[k*(int)R],(int)R+1,Y,1); }
    //     else { cblas_scopy((int)L,&X[-k],(int)R+1,Y,1); }
    // }
    // else
    // {
    //     if (k>0) { cblas_scopy((int)L,&X[k],(int)C+1,Y,1); }
    //     else { cblas_scopy((int)L,&X[-k*(int)C],(int)C+1,Y,1); }
    // }

    return 0;
}


int diag_d (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int k)
{
    const size_t K = (iscolmajor) ? R+1u : C+1u;

    if (k>=0 && (int)C>k)
    {
        const size_t L = (C-(size_t)k<R) ? C-(size_t)k : R;
        const int S = (iscolmajor) ? k*(int)R : k;
        
        X += S;
        for (size_t l=0u; l<L; ++l, X+=K, ++Y) { *Y = *X; }
    }
    else if (k<=0 && (int)R>-k)
    {
        const size_t L = ((int)R+k<(int)C) ? (size_t)((int)R+k) : C;
        const int S = (iscolmajor) ? -k : -k*(int)C;
        
        X += S;
        for (size_t l=0u; l<L; ++l, X+=K, ++Y) { *Y = *X; }
    }
    else
    {
        fprintf(stderr,"error in diag_d: k out of range [1-R C-1]\n"); return 1;
    }

    return 0;
}


int diag_c (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int k)
{
    const size_t K = (iscolmajor) ? 2u*R+1u : 2u*C+1u;

    if (k>=0 && (int)C>k)
    {
        const size_t L = (C-(size_t)k<R) ? C-(size_t)k : R;
        const int S = (iscolmajor) ? k*(int)R : k;
        
        X += 2*S;
        for (size_t l=0u; l<L; ++l, X+=K, ++Y) { *Y = *X; *++Y = *++X; }
    }
    else if (k<=0 && (int)R>-k)
    {
        const size_t L = ((int)R+k<(int)C) ? (size_t)((int)R+k) : C;
        const int S = (iscolmajor) ? -k : -k*(int)C;
        
        X += 2*S;
        for (size_t l=0u; l<L; ++l, X+=K, ++Y) { *Y = *X; *++Y = *++X; }
    }
    else
    {
        fprintf(stderr,"error in diag_c: k out of range [1-R C-1]\n"); return 1;
    }

    return 0;
}


int diag_z (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int k)
{
    const size_t K = (iscolmajor) ? 2u*R+1u : 2u*C+1u;

    if (k>=0 && (int)C>k)
    {
        const size_t L = (C-(size_t)k<R) ? C-(size_t)k : R;
        const int S = (iscolmajor) ? k*(int)R : k;
        
        X += 2*S;
        for (size_t l=0u; l<L; ++l, X+=K, ++Y) { *Y = *X; *++Y = *++X; }
    }
    else if (k<=0 && (int)R>-k)
    {
        const size_t L = ((int)R+k<(int)C) ? (size_t)((int)R+k) : C;
        const int S = (iscolmajor) ? -k : -k*(int)C;
        
        X += 2*S;
        for (size_t l=0u; l<L; ++l, X+=K, ++Y) { *Y = *X; *++Y = *++X; }
    }
    else
    {
        fprintf(stderr,"error in diag_z: k out of range [1-R C-1]\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
