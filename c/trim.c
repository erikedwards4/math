//Vec2vec operation.
//Trims each vector in X along dim.
//This is the operation used in trimmed mean, trimmed var, etc.
//The output Y is returned in sorted order.

//The bottom p% and the top q% of values are excluded.
//These are not percentiles, just percentages of data to exclude,
//although they are approximately equal to the percentiles.

//The inplace version still outputs Y, but modifies X during processing.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int trim_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q);
int trim_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q);

int trim_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q);
int trim_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q);


int trim_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q)
{
    if (dim>3) { fprintf(stderr,"error in trim_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>50.0f) { fprintf(stderr,"error in trim_s: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const float p1 = (p/100.0f)*(L-1), p2 = (1.0f-q/100.0f)*(L-1);
    size_t i1 = (size_t)floorf(p1), i2 = (size_t)ceilf(p2);
    const size_t Ly = i2-i1+1;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in trim_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        cblas_scopy((int)L,X,1,X1,1);
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trim_s: problem with LAPACKE function\n"); }
        cblas_scopy((int)Ly,&X1[i1],1,Y,1);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L, Y+=Ly)
            {
                cblas_scopy((int)L,X,1,X1,1);
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trim_s: problem with LAPACKE function\n"); }
                cblas_scopy((int)Ly,&X1[i1],1,Y,1);
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(Ly-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    cblas_scopy((int)L,X,(int)K,X1,1);
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trim_s: problem with LAPACKE function\n"); }
                    cblas_scopy((int)Ly,&X1[i1],1,Y,(int)K);
                }
            }
        }
    }

    free(X1);
    return 0;
}


int trim_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q)
{
    if (dim>3) { fprintf(stderr,"error in trim_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in trim_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const double p1 = (p/100.0)*(L-1), p2 = (1.0-q/100.0)*(L-1);
    size_t i1 = (size_t)floor(p1), i2 = (size_t)ceil(p2);
    const size_t Ly = i2-i1+1;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in trim_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        cblas_dcopy((int)L,X,1,X1,1);
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trim_d: problem with LAPACKE function\n"); }
        cblas_dcopy((int)Ly,&X1[i1],1,Y,1);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L, Y+=Ly)
            {
                cblas_dcopy((int)L,X,1,X1,1);
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trim_d: problem with LAPACKE function\n"); }
                cblas_dcopy((int)Ly,&X1[i1],1,Y,1);
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(Ly-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    cblas_dcopy((int)L,X,(int)K,X1,1);
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trim_d: problem with LAPACKE function\n"); }
                    cblas_dcopy((int)Ly,&X1[i1],1,Y,(int)K);
                }
            }
        }
    }

    free(X1);
    return 0;
}


int trim_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q)
{
    if (dim>3) { fprintf(stderr,"error in trim_inplace_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>50.0f) { fprintf(stderr,"error in trim_inplace_s: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const float p1 = (p/100.0f)*(L-1), p2 = (1.0f-q/100.0f)*(L-1);
    size_t i1 = (size_t)floorf(p1), i2 = (size_t)ceilf(p2);
    const size_t Ly = i2-i1+1;

    if (N==0) {}
    else if (L==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trim_s: problem with LAPACKE function\n"); }
        cblas_scopy((int)Ly,&X[i1],1,Y,1);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L-i1, Y+=Ly)
            {
                if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trim_s: problem with LAPACKE function\n"); }
                X += i1;
                cblas_scopy((int)Ly,X,1,Y,1);
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in trim_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(Ly-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    cblas_scopy((int)L,X,(int)K,X1,1);
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trim_s: problem with LAPACKE function\n"); }
                    cblas_scopy((int)Ly,&X1[i1],1,Y,(int)K);
                }
            }
            free(X1);
        }
    }

    return 0;
}


int trim_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q)
{
    if (dim>3) { fprintf(stderr,"error in trim_inplace_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in trim_inplace_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const double p1 = (p/100.0)*(L-1), p2 = (1.0-q/100.0)*(L-1);
    size_t i1 = (size_t)floor(p1), i2 = (size_t)ceil(p2);
    const size_t Ly = i2-i1+1;

    if (N==0) {}
    else if (L==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trim_d: problem with LAPACKE function\n"); }
        cblas_dcopy((int)Ly,&X[i1],1,Y,1);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L-i1, Y+=Ly)
            {
                if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trim_d: problem with LAPACKE function\n"); }
                X += i1;
                cblas_dcopy((int)Ly,X,1,Y,1);
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in trim_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(Ly-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    cblas_dcopy((int)L,X,(int)K,X1,1);
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trim_d: problem with LAPACKE function\n"); }
                    cblas_dcopy((int)Ly,&X1[i1],1,Y,(int)K);
                }
            }
            free(X1);
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
