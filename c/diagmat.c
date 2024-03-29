//Puts vector X on kth diagonal of matrix Y.

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int diagmat_s (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int k)
{
    const size_t L = (k>0) ? R : C;
    const size_t K = (iscolmajor) ? R+1u : C+1u;
    const int S = (iscolmajor) ? (k<0) ? -k : k*(int)R : (k<0) ? -k*(int)C : k;

    for (size_t n=R*C; n>0u; --n, ++Y) { *Y = 0.0f; }
    Y -= (int)(R*C) - S;
    for (size_t l=L; l>0u; --l, ++X, Y+=K) { *Y = *X; }
    
    return 0;
}


int diagmat_d (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int k)
{
    const size_t L = (k>0) ? R : C;
    const size_t K = (iscolmajor) ? R+1u : C+1u;
    const int S = (iscolmajor) ? (k<0) ? -k : k*(int)R : (k<0) ? -k*(int)C : k;

    for (size_t n=R*C; n>0u; --n, ++Y) { *Y = 0.0; }
    Y -= (int)(R*C) - S;
    for (size_t l=L; l>0u; --l, ++X, Y+=K) { *Y = *X; }

    return 0;
}


int diagmat_c (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int k)
{
    const size_t L = (k>0) ? R : C;
    const size_t K = (iscolmajor) ? R+1u : C+1u;
    const int S = (iscolmajor) ? (k<0) ? -2*k : 2*k*(int)R : (k<0) ? -2*k*(int)C : 2*k;

    for (size_t n=0u; n<2u*R*C; ++n, ++Y) { *Y = 0.0f; }
    Y -= 2*(int)(R*C) - S;
    for (size_t l=L; l>0u; --l, ++X, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }

    return 0;
}


int diagmat_z (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int k)
{
    const size_t L = (k>0) ? R : C;
    const size_t K = (iscolmajor) ? R+1u : C+1u;
    const int S = (iscolmajor) ? (k<0) ? -2*k : 2*k*(int)R : (k<0) ? -2*k*(int)C : 2*k;

    for (size_t n=0u; n<2u*R*C; ++n, ++Y) { *Y = 0.0; }
    Y -= 2*(int)(R*C) - S;
    for (size_t l=L; l>0u; --l, ++X, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
