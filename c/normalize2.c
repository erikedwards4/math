//Normalizes each vector in X along dim to have unit L2 norm.
//This operates in-place.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int normalize2_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int normalize2_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int normalize2_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int normalize2_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int normalize2_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in normalize2_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float nrm = 0.0f;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = (*X<0.0f) ? -1.0f : 1.0f; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X) { nrm += *X * *X; }
        nrm = sqrtf(nrm);
        for (size_t l=L; l>0u; --l) { *--X /= nrm; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                nrm = 0.0f;
                for (size_t l=L; l>0u; --l, ++X) { nrm += *X * *X; }
                nrm = sqrtf(nrm);
                X -= L;
                for (size_t l=L; l>0u; --l, ++X) { *X /= nrm; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    nrm = 0.0f;
                    for (size_t l=L; l>0u; --l, X+=K) { nrm += *X * *X; }
                    nrm = sqrtf(nrm);
                    for (size_t l=L; l>0u; --l) { X-=K; *X /= nrm; }
                }
            }
        }
    }

    return 0;
}


int normalize2_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in normalize2_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double nrm = 0.0;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = (*X<0.0) ? -1.0 : 1.0; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X) { nrm += *X * *X; }
        nrm = sqrt(nrm);
        for (size_t l=L; l>0u; --l) { *--X /= nrm; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                nrm = 0.0;
                for (size_t l=L; l>0u; --l, ++X) { nrm += *X * *X; }
                nrm = sqrt(nrm);
                X -= L;
                for (size_t l=L; l>0u; --l, ++X) { *X /= nrm; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    nrm = 0.0;
                    for (size_t l=L; l>0u; --l, X+=K) { nrm += *X * *X; }
                    nrm = sqrt(nrm);
                    for (size_t l=L; l>0u; --l) { X-=K; *X /= nrm; }
                }
            }
        }
    }

    return 0;
}


int normalize2_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in normalize2_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float nrm = 0.0f;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X)
        {
            nrm = sqrtf(*X**X + *(X+1)**(X+1));
            *X /= nrm; *++X /= nrm;
        }
    }
    else if (L==N)
    {
        for (size_t l=2u*L; l>0u; --l, ++X) { nrm += *X * *X; }
        nrm = sqrtf(nrm);
        for (size_t l=2u*L; l>0u; --l) { *--X /= nrm; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                nrm = 0.0f;
                for (size_t l=2u*L; l>0u; --l, ++X) { nrm += *X * *X; }
                nrm = sqrtf(nrm);
                X -= 2u*L;
                for (size_t l=2u*L; l>0u; --l, ++X) { *X /= nrm; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X+=2)
                {
                    nrm = 0.0f;
                    for (size_t l=L; l>0u; --l, X+=2u*K) { nrm += *X**X + *(X+1)**(X+1); }
                    nrm = sqrtf(nrm);
                    for (size_t l=L; l>0u; --l) { X-=2u*K; *X /= nrm; *(X+1) /= nrm; }
                }
            }
        }
    }

    return 0;
}


int normalize2_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in normalize2_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double nrm = 0.0;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X)
        {
            nrm = sqrt(*X**X + *(X+1)**(X+1));
            *X /= nrm; *++X /= nrm;
        }
    }
    else if (L==N)
    {
        for (size_t l=2u*L; l>0u; --l, ++X) { nrm += *X * *X; }
        nrm = sqrt(nrm);
        for (size_t l=2u*L; l>0u; --l) { *--X /= nrm; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                nrm = 0.0;
                for (size_t l=2u*L; l>0u; --l, ++X) { nrm += *X * *X; }
                nrm = sqrt(nrm);
                X -= 2u*L;
                for (size_t l=2u*L; l>0u; --l, ++X) { *X /= nrm; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X+=2)
                {
                    nrm = 0.0;
                    for (size_t l=L; l>0u; --l, X+=2u*K) { nrm += *X**X + *(X+1)**(X+1); }
                    nrm = sqrt(nrm);
                    for (size_t l=L; l>0u; --l) { X-=2u*K; *X /= nrm; *(X+1) /= nrm; }
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
