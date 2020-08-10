//Vec2scalar (reduction) operation.
//Gets the p-norm for each vector in X along dim.
//This is the Lp norm of each vector in X.
//For each vector, y = sum(|x|^p)^1/p.
//For complex case, output is real.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int normp_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p);
int normp_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p);
int normp_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p);
int normp_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p);


int normp_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p)
{
    if (dim>3) { fprintf(stderr,"error in normp_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ip = 1.0f / p;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        float sm = 0.0f;
        for (size_t l=0; l<L; ++l, ++X) { sm += powf(fabsf(*X),p); }
        *Y = powf(sm,ip);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            float sm;
            for (size_t v=0; v<V; ++v, ++Y)
            {
                sm = 0.0f;
                for (size_t l=0; l<L; ++l, ++X) { sm += powf(fabsf(*X),p); }
                *Y = powf(sm,ip);
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = powf(fabsf(*X),p); }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += powf(fabsf(*X),p); }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y = powf(*Y,ip); }
        }
        else
        {
            float sm;
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=0; l<L; ++l, X+=K) { sm += powf(fabsf(*X),p); }
                    *Y = powf(sm,ip);
                }
            }
        }
    }

    return 0;
}


int normp_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p)
{
    if (dim>3) { fprintf(stderr,"error in normp_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ip = 1.0 / p;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        double sm = 0.0;
        for (size_t l=0; l<L; ++l, ++X) { sm += pow(fabs(*X),p); }
        *Y = pow(sm,ip);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            double sm;
            for (size_t v=0; v<V; ++v, ++Y)
            {
                sm = 0.0;
                for (size_t l=0; l<L; ++l, ++X) { sm += pow(fabs(*X),p); }
                *Y = pow(sm,ip);
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = pow(fabs(*X),p); }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += pow(fabs(*X),p); }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y = pow(*Y,ip); }
        }
        else
        {
            double sm;
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=0; l<L; ++l, X+=K) { sm += pow(fabs(*X),p); }
                    *Y = pow(sm,ip);
                }
            }
        }
    }

    return 0;
}


int normp_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p)
{
    if (dim>3) { fprintf(stderr,"error in normp_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ip = 1.0f / p;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, X+=2, ++Y) { *Y = sqrtf(*X**X + *(X+1)**(X+1)); }
    }
    else if (L==N)
    {
        float sm = 0.0f;
        for (size_t l=0; l<L; ++l, X+=2) { sm += powf(sqrtf(*X**X + *(X+1)**(X+1)),p); }
        *Y = powf(sm,ip);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            float sm;
            for (size_t v=0; v<V; ++v, ++Y)
            {
                sm = 0.0f;
                for (size_t l=0; l<L; ++l, X+=2) { sm += powf(sqrtf(*X**X + *(X+1)**(X+1)),p); }
                *Y = powf(sm,ip);
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, X+=2, ++Y) { *Y = powf(sqrtf(*X**X + *(X+1)**(X+1)),p); }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, X+=2, ++Y) { *Y += powf(sqrtf(*X**X + *(X+1)**(X+1)),p); }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y = powf(*Y,ip); }
        }
        else
        {
            float sm;
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=0; l<L; ++l, X+=2*K) { sm += powf(sqrtf(*X**X + *(X+1)**(X+1)),p); }
                    *Y = powf(sm,ip);
                }
            }
        }
    }

    return 0;
}


int normp_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p)
{
    if (dim>3) { fprintf(stderr,"error in normp_z: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ip = 1.0 / p;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, X+=2, ++Y) { *Y = sqrt(*X**X + *(X+1)**(X+1)); }
    }
    else if (L==N)
    {
        double sm = 0.0;
        for (size_t l=0; l<L; ++l, X+=2) { sm += pow(sqrt(*X**X + *(X+1)**(X+1)),p); }
        *Y = pow(sm,ip);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            double sm;
            for (size_t v=0; v<V; ++v, ++Y)
            {
                sm = 0.0;
                for (size_t l=0; l<L; ++l, X+=2) { sm += pow(sqrt(*X**X + *(X+1)**(X+1)),p); }
                *Y = pow(sm,ip);
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, X+=2, ++Y) { *Y = pow(sqrt(*X**X + *(X+1)**(X+1)),p); }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, X+=2, ++Y) { *Y += pow(sqrt(*X**X + *(X+1)**(X+1)),p); }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y = pow(*Y,ip); }
        }
        else
        {
            double sm;
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=0; l<L; ++l, X+=2*K) { sm += pow(sqrt(*X**X + *(X+1)**(X+1)),p); }
                    *Y = pow(sm,ip);
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
