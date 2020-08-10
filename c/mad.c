//Vec2scalar (reduction) operation.
//Gets MAD (median absolute deviation from the median) for each vector in X along dim.

//The in-place versions still return the same Y, but modify X during processing.
//However, it turns out to be almost the identical speed for matrices.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int mad_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mad_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);

int mad_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mad_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int mad_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mad_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t i2 = L/2;
    float med;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in mad_s: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_s: problem with LAPACKE function\n"); }
        X1 += i2;
        med = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
        X1 -= i2;
        for (size_t l=0; l<L; ++l, ++X1) { *X1 = fabsf(*X1-med); }
        X1 -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_s: problem with LAPACKE function\n"); }
        X1 += i2;
        *Y = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
        X1 -= i2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X1-=i2, ++Y)
            {
                for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_s: problem with LAPACKE function\n"); }
                X1 += i2;
                med = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=0; l<L; ++l, ++X1) { *X1 = fabsf(*X1-med); }
                X1 -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_s: problem with LAPACKE function\n"); }
                X1 += i2;
                *Y = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, X1-=i2, ++Y)
                {
                    for (size_t l=0; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_s: problem with LAPACKE function\n"); }
                    X1 += i2;
                    med = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=0; l<L; ++l, ++X1) { *X1 = fabsf(*X1-med); }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_s: problem with LAPACKE function\n"); }
                    X1 += i2;
                    *Y = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
                }
            }
        }
    }

    free(X1);
    return 0;
}


int mad_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mad_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t i2 = L/2;
    double med;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in mad_d: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_d: problem with LAPACKE function\n"); }
        X1 += i2;
        med = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
        X1 -= i2;
        for (size_t l=0; l<L; ++l, ++X1) { *X1 = fabs(*X1-med); }
        X1 -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_d: problem with LAPACKE function\n"); }
        X1 += i2;
        *Y = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
        X1 -= i2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X1-=i2, ++Y)
            {
                for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_d: problem with LAPACKE function\n"); }
                X1 += i2;
                med = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=0; l<L; ++l, ++X1) { *X1 = fabs(*X1-med); }
                X1 -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_d: problem with LAPACKE function\n"); }
                X1 += i2;
                *Y = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, X1-=i2, ++Y)
                {
                    for (size_t l=0; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_d: problem with LAPACKE function\n"); }
                    X1 += i2;
                    med = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=0; l<L; ++l, ++X1) { *X1 = fabs(*X1-med); }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_d: problem with LAPACKE function\n"); }
                    X1 += i2;
                    *Y = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
                }
            }
        }
    }

    free(X1);
    return 0;
}


int mad_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mad_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t i2 = L/2;
    float med;
    
    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in mad_inplace_s: problem with LAPACKE function\n"); }
        X += i2;
        med = (L%2) ? *X : 0.5f*(*X + *(X-1));
        X -= i2;
        for (size_t l=0; l<L; ++l, ++X) { *X = fabsf(*X-med); }
        X -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in mad_inplace_s: problem with LAPACKE function\n"); }
        X += i2;
        *Y = (L%2) ? *X : 0.5f*(*X + *(X-1));
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L-i2, ++Y)
            {
                if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in mad_inplace_s: problem with LAPACKE function\n"); }
                X += i2; med = (L%2) ? *X : 0.5f*(*X + *(X-1)); X -= i2;
                for (size_t l=0; l<L; ++l, ++X) { *X = fabsf(*X-med); }
                X -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in mad_inplace_s: problem with LAPACKE function\n"); }
                X += i2;
                *Y = (L%2) ? *X : 0.5f*(*X + *(X-1));
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in mad_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, X1-=i2, ++Y)
                {
                    for (size_t l=0; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_inplace_s: problem with LAPACKE function\n"); }
                    X1 += i2; med = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1)); X1 -= i2;
                    for (size_t l=0; l<L; ++l, ++X1) { *X1 = fabsf(*X1-med); }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_inplace_s: problem with LAPACKE function\n"); }
                    X1 += i2;
                    *Y = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
                }
            }
            free(X1);
        }
    }

    return 0;
}


int mad_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mad_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t i2 = L/2;
    double med;
    
    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in mad_inplace_d: problem with LAPACKE function\n"); }
        X += i2;
        med = (L%2) ? *X : 0.5*(*X + *(X-1));
        X -= i2;
        for (size_t l=0; l<L; ++l, ++X) { *X = fabs(*X-med); }
        X -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in mad_inplace_d: problem with LAPACKE function\n"); }
        X += i2;
        *Y = (L%2) ? *X : 0.5*(*X + *(X-1));
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L-i2, ++Y)
            {
                if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in mad_inplace_d: problem with LAPACKE function\n"); }
                X += i2; med = (L%2) ? *X : 0.5*(*X + *(X-1)); X -= i2;
                for (size_t l=0; l<L; ++l, ++X) { *X = fabs(*X-med); }
                X -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in mad_inplace_d: problem with LAPACKE function\n"); }
                X += i2;
                *Y = (L%2) ? *X : 0.5*(*X + *(X-1));
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in mad_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, X1-=i2, ++Y)
                {
                    for (size_t l=0; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_inplace_d: problem with LAPACKE function\n"); }
                    X1 += i2; med = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1)); X1 -= i2;
                    for (size_t l=0; l<L; ++l, ++X1) { *X1 = fabs(*X1-med); }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_inplace_d: problem with LAPACKE function\n"); }
                    X1 += i2;
                    *Y = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
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
