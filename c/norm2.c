//Vec2scalar (reduction) operation.
//Gets L2-norm (Euclidean) for each vector in X along dim.
//This is the root-sum-square for each vector in X.
//Note that this is not the Frobenius matrix norm of X (see matnorm.c).

//For complex case, output is real.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int norm2_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int norm2_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int norm2_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int norm2_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int norm2_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in norm2_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X[n]*X[n]; }
    }
    else if (L==N)
    {
        *Y = cblas_snrm2((int)L,X,1);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, X+=L)
            {
                *Y++ = cblas_snrm2((int)L,X,1);
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; v++, X++) { *Y++ = *X**X; }
            Y -= V;
            for (size_t l=1; l<L; l++, Y-=V)
            {
                for (size_t v=0; v<V; v++, X++) { *Y++ += *X**X; }
            }
            for (size_t v=0; v<V; v++, Y++) { *Y = sqrtf(*Y); }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, X++)
                {
                    *Y++ = cblas_snrm2((int)L,X,(int)K);
                }
            }
        }
    }

    return 0;
}


int norm2_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in norm2_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X[n]*X[n]; }
    }
    else if (L==N)
    {
        *Y = cblas_dnrm2((int)L,X,1);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, X+=L)
            {
                *Y++ = cblas_dnrm2((int)L,X,1);
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; v++, X++) { *Y++ = *X**X; }
            Y -= V;
            for (size_t l=1; l<L; l++, Y-=V)
            {
                for (size_t v=0; v<V; v++, X++) { *Y++ += *X**X; }
            }
            for (size_t v=0; v<V; v++, Y++) { *Y = sqrt(*Y); }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, X++)
                {
                    *Y++ = cblas_dnrm2((int)L,X,(int)K);
                }
            }
        }
    }

    return 0;
}


int norm2_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in norm2_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; n++, X+=2) { *Y++ = *X**X + *(X+1)**(X+1); }
    }
    else if (L==N)
    {
        *Y = cblas_scnrm2((int)L,X,1);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, X+=2*L)
            {
                *Y++ = cblas_scnrm2((int)L,X,1);
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; v++, X+=2) { *Y++ = *X**X + *(X+1)**(X+1); }
            Y -= V;
            for (size_t l=1; l<L; l++, Y-=V)
            {
                for (size_t v=0; v<V; v++, X+=2) { *Y++ += *X**X + *(X+1)**(X+1); }
            }
            for (size_t v=0; v<V; v++, Y++) { *Y = sqrtf(*Y); }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; b++, X+=2)
                {
                    *Y++ = cblas_scnrm2((int)L,X,(int)K);
                }
            }
        }
    }

    return 0;
}


int norm2_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in norm2_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; n++, X+=2) { *Y++ = *X**X + *(X+1)**(X+1); }
    }
    else if (L==N)
    {
        *Y = cblas_dznrm2((int)L,X,1);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, X+=2*L)
            {
                *Y++ = cblas_dznrm2((int)L,X,1);
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; v++, X+=2) { *Y++ = *X**X + *(X+1)**(X+1); }
            Y -= V;
            for (size_t l=1; l<L; l++, Y-=V)
            {
                for (size_t v=0; v<V; v++, X+=2) { *Y++ += *X**X + *(X+1)**(X+1); }
            }
            for (size_t v=0; v<V; v++, Y++) { *Y = sqrt(*Y); }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; b++, X+=2)
                {
                    *Y++ = cblas_dznrm2((int)L,X,(int)K);
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
