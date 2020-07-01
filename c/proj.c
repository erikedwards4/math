//Gets proj part of complex-valued input X.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <complex.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int proj_c (float *Y, const float *X, const int N);
int proj_z (double *Y, const double *X, const int N);

int proj_inplace_c (float *X, const int N);
int proj_inplace_z (double *X, const int N);


int proj_c (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in proj_s: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex float y;

    for (int n=0; n<2*N; n+=2)
    {
        y = cprojf(X[n]+1.0if*X[n+1]);
        memcpy(&Y[n],(float *)&y,2*sizeof(float));
    }

    return 0;
}


int proj_z (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in proj_z: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex double y;

    for (int n=0; n<2*N; n+=2)
    {
        y = cproj(X[n]+1.0i*X[n+1]);
        memcpy(&Y[n],(double *)&y,2*sizeof(double));
    }
    
    return 0;
}


int proj_inplace_c (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in proj_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex float y;

    for (int n=0; n<2*N; n+=2)
    {
        y = cprojf(X[n]+1.0if*X[n+1]);
        memcpy(&X[n],(float *)&y,2*sizeof(float));
    }
    
    return 0;
}


int proj_inplace_z (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in proj_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex double y;

    for (int n=0; n<2*N; n+=2)
    {
        y = cproj(X[n]+1.0i*X[n+1]);
        memcpy(&X[n],(double *)&y,2*sizeof(double));
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
