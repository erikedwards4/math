//Gets inverse hyberbolic tangent (atanh) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int atanh_s (float *Y, const float *X, const int N);
int atanh_d (double *Y, const double *X, const int N);
int atanh_c (float *Y, const float *X, const int N);
int atanh_z (double *Y, const double *X, const int N);

int atanh_inplace_s (float *X, const int N);
int atanh_inplace_d (double *X, const int N);
int atanh_inplace_c (float *X, const int N);
int atanh_inplace_z (double *X, const int N);


int atanh_s (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atanh_s: N (num elements X) must be nonnegative\n"); return 1; }

    //for (int n=0; n<N; n++) { Y[n] = atanhf(X[n]); }
    for (int n=0; n<N; n++) { Y[n] = 0.5f * logf((1.0f+X[n])/(1.0f-X[n])); }

    return 0;
}


int atanh_d (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atanh_d: N (num elements X) must be nonnegative\n"); return 1; }

    //for (int n=0; n<N; n++) { Y[n] = atanh(X[n]); }
    for (int n=0; n<N; n++) { Y[n] = 0.5 * log((1.0+X[n])/(1.0-X[n])); }
    
    return 0;
}


int atanh_c (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atanh_c: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex float x, y;

    for (int n2=0; n2<2*N; n2+=2)
    {
        //y = catanhf(X[n2]+1.0if*X[n2+1]);
        x = X[n2] + 1.0if*X[n2+1];
        y = 0.5f * clogf((1.0f+x)/(1.0f-x));
        memcpy(&Y[n2],(float *)&y,2*sizeof(float));
    }
    
    return 0;
}


int atanh_z (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atanh_z: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex double x, y;

    for (int n2=0; n2<2*N; n2+=2)
    {
        //y = catanh(X[n2]+1.0i*X[n2+1]);
        x = X[n2] + 1.0i*X[n2+1];
        y = 0.5 * clog((1.0+x)/(1.0-x));
        memcpy(&Y[n2],(double *)&y,2*sizeof(double));
    }
    
    return 0;
}


int atanh_inplace_s (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atanh_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    //for (int n=0; n<N; n++) { X[n] = atanhf(X[n]); }
    for (int n=0; n<N; n++) { X[n] = 0.5f * logf((1.0f+X[n])/(1.0f-X[n])); }

    return 0;
}


int atanh_inplace_d (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atanh_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    //for (int n=0; n<N; n++) { X[n] = atanh(X[n]); }
    for (int n=0; n<N; n++) { X[n] = 0.5 * log((1.0+X[n])/(1.0-X[n])); }
    
    return 0;
}


int atanh_inplace_c (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atanh_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex float x;

    for (int n2=0; n2<2*N; n2+=2)
    {
        //x = catanhf(X[n2]+1.0if*X[n2+1]);
        x = X[n2] + 1.0if*X[n2+1];
        x = 0.5f * clogf((1.0f+x)/(1.0f-x));
        memcpy(&X[n2],(float *)&x,2*sizeof(float));
    }
    
    return 0;
}


int atanh_inplace_z (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atanh_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }
    
    _Complex double x;

    for (int n2=0; n2<2*N; n2+=2)
    {
        //x = catanh(X[n2]+1.0i*X[n2+1]);
        x = X[n2] + 1.0i*X[n2+1];
        x = 0.5 * clog((1.0+x)/(1.0-x));
        memcpy(&X[n2],(double *)&x,2*sizeof(double));
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
