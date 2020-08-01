//Vec2scalar (reduction) operation.
//Gets root-mean-square (RMS) for each vector in X along dim.
//This is the square root of the mean of squares for each vector.
//This is the generalized mean with p=2.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int rms_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int rms_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int rms_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int rms_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int rms_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in rms_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f / sqrtf(L);

    if (N==0) {}
    else if (L==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
        //*Y = 0.0f; for (size_t l=0; l<L; ++l) { *Y += X[l]*X[l]; }
        //*Y = sqrtf(*Y/L); //this is faster for small L, but slower (and less numerically accurate) for large L
        *Y = cblas_snrm2((int)L,X,1) * den;
        //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L, ++Y)
            {
                *Y++ = cblas_snrm2((int)L,X,1) * den;
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X) { *Y++ = *X**X; }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X) { *Y++ += *X**X; }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y = sqrtf(*Y/L); }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X)
                {
                    *Y++ = cblas_snrm2((int)L,X,(int)K);
                }
            }
        }
    }

    return 0;
}


int rms_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in rms_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0 / sqrt(L);

    if (N==0) {}
    else if (L==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        *Y = cblas_dnrm2((int)L,X,1) * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L, ++Y)
            {
                *Y++ = cblas_dnrm2((int)L,X,1) * den;
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X) { *Y++ = *X**X; }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X) { *Y++ += *X**X; }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y = sqrt(*Y/L); }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X)
                {
                    *Y++ = cblas_dnrm2((int)L,X,(int)K);
                }
            }
        }
    }

    return 0;
}


int rms_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in rms_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f / sqrtf(L);

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; n+=2) { *Y++ = sqrtf(X[n]*X[n]+X[n+1]*X[n+1]); }
    }
    else if (L==N)
    {
        *Y = cblas_scnrm2((int)L,X,1) * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=2*L)
            {
                *Y++ = cblas_scnrm2((int)L,X,1) * den;
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, X+=2) { *Y++ = *X**X + *(X+1)**(X+1); }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, X+=2) { *Y++ += *X**X + *(X+1)**(X+1); }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y = sqrtf(*Y); }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X+=2)
                {
                    *Y++ = cblas_scnrm2((int)L,X,(int)K) * den;
                }
            }
        }
    }

    return 0;
}


int rms_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in rms_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0 / sqrt(L);

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; n+=2) { *Y++ = sqrt(X[n]*X[n]+X[n+1]*X[n+1]); }
    }
    else if (L==N)
    {
        *Y = cblas_dznrm2((int)L,X,1) * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=2*L)
            {
                *Y++ = cblas_dznrm2((int)L,X,1) * den;
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, X+=2) { *Y++ = *X**X + *(X+1)**(X+1); }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, X+=2) { *Y++ += *X**X + *(X+1)**(X+1); }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y = sqrt(*Y); }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X+=2)
                {
                    *Y++ = cblas_dznrm2((int)L,X,(int)K) * den;
                }
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
