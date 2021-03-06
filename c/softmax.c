//Vec2vec operation
//Softmax of each vector in X.
//This is a common layer-wise activation function in NNs.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int softmax_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int softmax_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);

int softmax_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int softmax_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int softmax_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in softmax_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float sm = 0.0f;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = 1.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = expf(*X); sm += *Y; }
        for (size_t l=0u; l<L; ++l) { *--Y /= sm; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                sm = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = expf(*X); sm += *Y; }
                Y -= L;
                for (size_t l=0u; l<L; ++l, ++Y) { *Y /= sm; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1), Y+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=K, Y+=K) { *Y = expf(*X); sm += *Y; }
                    for (size_t l=0u; l<L; ++l) { Y-=K; *Y /= sm; }
                }
            }
        }
    }

    return 0;
}


int softmax_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in softmax_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double sm = 0.0;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = 1.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = exp(*X); sm += *Y; }
        for (size_t l=0u; l<L; ++l) { *--Y /= sm; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                sm = 0.0;
                for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = exp(*X); sm += *Y; }
                Y -= L;
                for (size_t l=0u; l<L; ++l, ++Y) { *Y /= sm; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1), Y+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=K, Y+=K) { *Y = exp(*X); sm += *Y; }
                    for (size_t l=0u; l<L; ++l) { Y-=K; *Y /= sm; }
                }
            }
        }
    }

    return 0;
}


int softmax_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in softmax_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float sm = 0.0f;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X) { *X = 1.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { *X = expf(*X); sm += *X; }
        for (size_t l=0u; l<L; ++l) { *--X /= sm; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                sm = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X) { *X = expf(*X); sm += *X; }
                X -= L;
                for (size_t l=0u; l<L; ++l, ++X) { *X /= sm; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, ++X)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=K) { *X = expf(*X); sm += *X; }
                    for (size_t l=0u; l<L; ++l) { X-=K; *X /= sm; }
                }
            }
        }
    }

    return 0;
}


int softmax_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in softmax_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double sm = 0.0;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X) { *X = 1.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { *X = exp(*X); sm += *X; }
        for (size_t l=0u; l<L; ++l) { *--X /= sm; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                sm = 0.0;
                for (size_t l=0u; l<L; ++l, ++X) { *X = exp(*X); sm += *X; }
                X -= L;
                for (size_t l=0u; l<L; ++l, ++X) { *X /= sm; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, ++X)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=K) { *X = exp(*X); sm += *X; }
                    for (size_t l=0u; l<L; ++l) { X-=K; *X /= sm; }
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
