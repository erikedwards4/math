//Gets count of nonzero values for each row or col of X according to dim.
//This is also the Hamming norm (or L0 norm) of each vector in X.

//For complex case, real and imag parts calculated separately.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cnt_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int cnt_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int cnt_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int cnt_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);

int cnt_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int cnt_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int cnt_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int cnt_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int cnt_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cnt_s: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t G = N/(B*L);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (size_t g=0; g<G; g++, X+=B*(L-J))
    {
        for (size_t b=0; b<B; b++, X+=J, Y++)
        {
            *Y = 0.0f;
            for (size_t l=0; l<L; l++) { *Y += (float)(*(X+l*K)!=0.0f); }
        }
    }

    return 0;
}


int cnt_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cnt_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t G = N/(B*L);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (size_t g=0; g<G; g++, X+=B*(L-J))
    {
        for (size_t b=0; b<B; b++, X+=J, Y++)
        {
            *Y = 0.0;
            for (size_t l=0; l<L; l++) { *Y += (double)(*(X+l*K)!=0.0); }
        }
    }

    return 0;
}


int cnt_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cnt_c: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t G = N/(B*L);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (size_t g=0; g<G; g++, X+=2*B*(L-J))
    {
        for (size_t b=0; b<B; b++, X+=2*J, Y++)
        {
            *Y = 0.0f;
            for (size_t l=0; l<2*L; l+=2) { *Y += (float)(*(X+l*K)!=0.0f || *(X+l*K+1)!=0.0f); }
        }
    }

    return 0;
}


int cnt_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cnt_z: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t G = N/(B*L);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (size_t g=0; g<G; g++, X+=2*B*(L-J))
    {
        for (size_t b=0; b<B; b++, X+=2*J, Y++)
        {
            *Y = 0.0;
            for (size_t l=0; l<2*L; l+=2) { *Y += (double)(*(X+l*K)!=0.0 || *(X+l*K+1)!=0.0); }
        }
    }

    return 0;
}


int cnt_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cnt_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t G = N/(B*L);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    for (size_t g=0, n=0, v=0; g<G; g++, n+=B*(L-J))
    {
        for (size_t b=0; b<B; b++, n+=J, v++)
        {
            float sm = 0.0f;
            for (size_t l=0; l<L; l++) { sm += (float)(X[n+l*K]!=0.0f); }
            X[v] = sm;
        }
    }

    return 0;
}


int cnt_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cnt_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t G = N/(B*L);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    for (size_t g=0, n=0, v=0; g<G; g++, n+=B*(L-J))
    {
        for (size_t b=0; b<B; b++, n+=J, v++)
        {
            double sm = 0.0;
            for (size_t l=0; l<L; l++) { sm += (double)(X[n+l*K]!=0.0); }
            X[v] = sm;
        }
    }

    return 0;
}


int cnt_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cnt_inplace_c: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t G = N/(B*L);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    for (size_t g=0, n=0, v=0; g<G; g++, n+=2*B*(L-J))
    {
        for (size_t b=0; b<B; b++, n+=2*J, v++)
        {
            float sm = 0.0f;
            for (size_t l=0; l<2*L; l+=2) { sm += (float)(X[n+l*K]!=0.0f || X[n+l*K+1]!=0.0f); }
            X[v] = sm;
        }
    }

    return 0;
}


int cnt_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cnt_inplace_z: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t G = N/(B*L);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    for (size_t g=0, n=0, v=0; g<G; g++, n+=2*B*(L-J))
    {
        for (size_t b=0; b<B; b++, n+=2*J, v++)
        {
            double sm = 0.0;
            for (size_t l=0; l<2*L; l+=2) { sm += (double)(X[n+l*K]!=0.0 || X[n+l*K+1]!=0.0); }
            X[v] = sm;
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
