//Vec2scalar (reduction) operation.
//Gets pth percentile for each vector in X along dim.

//The in-place versions still return the same Y, but modify X during processing.
//However, it turns out to be almost the identical speed for matrices.

#include <stdio.h>
#include <math.h>
// #include <lapacke.h>
#include "codee_math.h"
#include "kselect.c"
#include "extremum.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int prctile_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p)
{
    if (dim>3u) { fprintf(stderr,"error in prctile_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>100.0f) { fprintf(stderr,"error in prctile_s: p must be in [0 100]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Prep interpolation
    const float p1 = (p/100.0f)*(float)(L-1u);
    const size_t i1 = (p<100.0f) ? (size_t)floorf(p1) : L-2u;
    const size_t i2 = i1 + 1u;
    const float w2 = (p<100.0f) ? p1-floorf(p1) : 1.0f;
    const float w1 = 1.0f - w2;
    float x1, x2;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in prctile_s: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        x2 = kselect_s(X1,L,i2,1);
        // x1 = kselect_s(X1,i2,i1,1);
        x1 = extremum_s(X,i2,0);
        *Y = w1*x1 + w2*x2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, ++Y)
            {
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                x2 = kselect_s(X1,L,i2,1);
                x1 = extremum_s(X,i2,0);
                *Y = w1*x1 + w2*x2;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    x2 = kselect_s(X1,L,i2,1);
                    x1 = extremum_s(X1,i2,0);
                    *Y = w1*x1 + w2*x2;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int prctile_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p)
{
    if (dim>3u) { fprintf(stderr,"error in prctile_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>100.0) { fprintf(stderr,"error in prctile_d: p must be in [0 100]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Prep interpolation
    const double p1 = (p/100.0)*(double)(L-1u);
    const size_t i1 = (p<100.0) ? (size_t)floor(p1) : L-2u;
    const size_t i2 = i1 + 1u;
    const double w2 = (p<100.0) ? p1-floor(p1) : 1.0;
    const double w1 = 1.0 - w2;
    double x1, x2;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in prctile_d: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        x2 = kselect_d(X1,L,i2,1);
        x1 = extremum_d(X,i2,0);
        *Y = w1*x1 + w2*x2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, ++Y)
            {
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                x2 = kselect_d(X1,L,i2,1);
                x1 = extremum_d(X,i2,0);
                *Y = w1*x1 + w2*x2;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    x2 = kselect_d(X1,L,i2,1);
                    x1 = extremum_d(X1,i2,0);
                    *Y = w1*x1 + w2*x2;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int prctile_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p)
{
    if (dim>3u) { fprintf(stderr,"error in prctile_inplace_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>100.0f) { fprintf(stderr,"error in prctile_inplace_s: p must be in [0 100]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Prep interpolation
    const float p1 = (p/100.0f)*(float)(L-1u);
    const size_t i1 = (p<100.0f) ? (size_t)floorf(p1) : L-2u;
    const size_t i2 = i1 + 1u;
    const float w2 = (p<100.0f) ? p1-floorf(p1) : 1.0f;
    const float w1 = 1.0f - w2;
    float x1, x2;

    // struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        x2 = kselect_s(X,L,i2,1);
        //x1 = kselect_s(X,i1,i1,1);
        x1 = extremum_s(X,i2,0);
        *Y = w1*x1 + w2*x2;
        //if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in prctile_inplace_s: problem with LAPACKE function\n"); }
        //quicksort_s(X,L,1);
        //partsort_s(X,L,i1+1u,1);
        //X += i1;
        //*Y = w1**X + w2**(X+1);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=L, ++Y)
            {
                x2 = kselect_s(X,L,i2,1);
                x1 = extremum_s(X,i2,0);
                *Y = w1*x1 + w2*x2;
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in prctile_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    x2 = kselect_s(X1,L,i2,1);
                    x1 = extremum_s(X1,i2,0);
                    *Y = w1*x1 + w2*x2;
                }
            }
            free(X1);
        }
    }
    
    // clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int prctile_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p)
{
    if (dim>3u) { fprintf(stderr,"error in prctile_inplace_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>100.0) { fprintf(stderr,"error in prctile_inplace_d: p must be in [0 100]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Prep interpolation
    const double p1 = (p/100.0)*(double)(L-1u);
    const size_t i1 = (p<100.0) ? (size_t)floor(p1) : L-2u;
    const size_t i2 = i1 + 1u;
    const double w2 = (p<100.0) ? p1-floor(p1) : 1.0;
    const double w1 = 1.0 - w2;
    double x1, x2;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        x2 = kselect_d(X,L,i2,1);
        x1 = extremum_d(X,i2,0);
        *Y = w1*x1 + w2*x2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=L, ++Y)
            {
                x2 = kselect_d(X,L,i2,1);
                x1 = extremum_d(X,i2,0);
                *Y = w1*x1 + w2*x2;
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in prctile_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    x2 = kselect_d(X1,L,i2,1);
                    x1 = extremum_d(X1,i2,0);
                    *Y = w1*x1 + w2*x2;
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
