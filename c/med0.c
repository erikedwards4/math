//Zeros (subtracts) the median of each vector in X along dim.
//For each vector, estimates the median and then subtracts it from each element.
//This operates in-place.

#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int med0_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int med0_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int med0_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in med0_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t i2 = L/2;
    float med;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in med0_s: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X) { *X = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in med0_s: problem with LAPACKE function\n"); }
        X1 += i2;
        med = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
        X1 -= i2;
        for (size_t l=0; l<L; ++l) { --X; *X -= med; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L; X -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in med0_s: problem with LAPACKE function\n"); }
                X1 += i2;
                med = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=0; l<L; ++l, ++X) { *X -= med; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X)
                {
                    for (size_t l=0; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in med0_s: problem with LAPACKE function\n"); }
                    X1 += i2;
                    med = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=0; l<L; ++l) { X-=K; *X -= med; }
                }
            }
        }
    }

    free(X1);
    return 0;
}


int med0_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in med0_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t i2 = L/2;
    double med;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in med0_d: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X) { *X = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in med0_d: problem with LAPACKE function\n"); }
        X1 += i2;
        med = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
        X1 -= i2;
        for (size_t l=0; l<L; ++l) { --X; *X -= med; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L; X -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in med0_d: problem with LAPACKE function\n"); }
                X1 += i2;
                med = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=0; l<L; ++l, ++X) { *X -= med; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X)
                {
                    for (size_t l=0; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in med0_d: problem with LAPACKE function\n"); }
                    X1 += i2;
                    med = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=0; l<L; ++l) { X-=K; *X -= med; }
                }
            }
        }
    }

    free(X1);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
