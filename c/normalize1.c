//Normalizes each vector in X along dim to have unit L1 norm.
//This operates in-place.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int normalize1_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int normalize1_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int normalize1_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int normalize1_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int normalize1_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in normalize1_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float nrm = 0.0f;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X) { *X = (*X<0.0f) ? -1.0f : 1.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X) { nrm += fabsf(*X); }
        for (size_t l=0; l<L; ++l) { *--X /= nrm; }
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
                nrm = 0.0f;
                for (size_t l=0; l<L; ++l, ++X) { nrm += fabsf(*X); }
                X -= L;
                for (size_t l=0; l<L; ++l, ++X) { *X /= nrm; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X)
                {
                    nrm = 0.0f;
                    for (size_t l=0; l<L; ++l, X+=K) { nrm += fabsf(*X); }
                    for (size_t l=0; l<L; ++l) { X-=K; *X /= nrm; }
                }
            }
        }
    }

    return 0;
}


int normalize1_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in normalize1_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double nrm = 0.0;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X) { *X = (*X<0.0) ? -1.0 : 1.0; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X) { nrm += fabs(*X); }
        for (size_t l=0; l<L; ++l) { *--X /= nrm; }
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
                nrm = 0.0;
                for (size_t l=0; l<L; ++l, ++X) { nrm += fabs(*X); }
                X -= L;
                for (size_t l=0; l<L; ++l, ++X) { *X /= nrm; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X)
                {
                    nrm = 0.0;
                    for (size_t l=0; l<L; ++l, X+=K) { nrm += fabs(*X); }
                    for (size_t l=0; l<L; ++l) { X-=K; *X /= nrm; }
                }
            }
        }
    }

    return 0;
}


int normalize1_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in normalize1_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float nrm = 0.0f;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X)
        {
            nrm = sqrtf(*X**X + *(X+1)**(X+1));
            *X /= nrm; *++X /= nrm;
        }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, X+=2) { nrm += sqrtf(*X**X + *(X+1)**(X+1)); }
        for (size_t l=0; l<2*L; ++l) { *--X /= nrm; }
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
                nrm = 0.0f;
                for (size_t l=0; l<L; ++l, X+=2) { nrm += sqrtf(*X**X + *(X+1)**(X+1)); }
                X -= 2*L;
                for (size_t l=0; l<2*L; ++l, ++X) { *X /= nrm; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X+=2)
                {
                    nrm = 0.0f;
                    for (size_t l=0; l<L; ++l, X+=2*K) { nrm += sqrtf(*X**X + *(X+1)**(X+1)); }
                    for (size_t l=0; l<L; ++l) { X-=2*K; *X /= nrm; *(X+1) /= nrm; }
                }
            }
        }
    }

    return 0;
}


int normalize1_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in normalize1_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double nrm = 0.0;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X)
        {
            nrm = sqrt(*X**X + *(X+1)**(X+1));
            *X /= nrm; *++X /= nrm;
        }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, X+=2) { nrm += sqrt(*X**X + *(X+1)**(X+1)); }
        for (size_t l=0; l<2*L; ++l) { *--X /= nrm; }
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
                nrm = 0.0;
                for (size_t l=0; l<L; ++l, X+=2) { nrm += sqrt(*X**X + *(X+1)**(X+1)); }
                X -= 2*L;
                for (size_t l=0; l<2*L; ++l, ++X) { *X /= nrm; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X+=2)
                {
                    nrm = 0.0;
                    for (size_t l=0; l<L; ++l, X+=2*K) { nrm += sqrt(*X**X + *(X+1)**(X+1)); }
                    for (size_t l=0; l<L; ++l) { X-=2*K; *X /= nrm; *(X+1) /= nrm; }
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
