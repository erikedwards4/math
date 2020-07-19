//Gets L1 norm (sum of absolute-values) for each vector in X along dim.
//This is the 'taxicab' (L1) norm of each vector in X.

//For complex case, output is real. This does not use cblas_scasum or cblas_dzasum,
//because those sum |xr|+|xi|, which is not the definition of the L1 norm.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int norm1_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int norm1_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int norm1_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int norm1_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int norm1_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in norm1_s: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (L==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = fabsf(X[n]); }
    }
    else if (L==N) { *Y = cblas_sasum((int)L,X,1); }
    else if (SH==1)
    {
        const size_t V = N/L;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t v=0; v<V; v++, X+=L)
            {
                *Y++ = cblas_sasum((int)L,X,1);
            }
        }
        else
        {
            // for (size_t v=0; v<V; v++, X++)
            // {
            //     *Y++ = cblas_sasum((int)L,X,(int)V);
            // }
            for (size_t v=0; v<V; v++, X++) { *Y++ = fabsf(*X); }
            Y -= V;
            for (size_t l=1; l<L; l++, Y-=V)
            {
                for (size_t v=0; v<V; v++, X++) { *Y++ += fabsf(*X); }
            }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, X+=L)
            {
                *Y++ = cblas_sasum((int)L,X,1);
            }
        }
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=J)
            {
                *Y++ = cblas_sasum((int)L,X,(int)K);
            }
        }
    }

    return 0;
}


int norm1_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in norm1_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (L==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = fabs(X[n]); }
    }
    else if (L==N) { *Y = cblas_dasum((int)L,X,1); }
    else if (SH==1)
    {
        const size_t V = N/L;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t v=0; v<V; v++, X+=L)
            {
                *Y++ = cblas_dasum((int)L,X,1);
            }
        }
        else
        {
            for (size_t v=0; v<V; v++, X++) { *Y++ = fabs(*X); }
            Y -= V;
            for (size_t l=1; l<L; l++, Y-=V)
            {
                for (size_t v=0; v<V; v++, X++) { *Y++ += fabs(*X); }
            }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, X+=L)
            {
                *Y++ = cblas_dasum((int)L,X,1);
            }
        }
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=J)
            {
                *Y++ = cblas_dasum((int)L,X,(int)K);
            }
        }
    }

    return 0;
}


int norm1_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in norm1_c: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (L==1)
    {
        for (size_t n=0; n<N; n+=2) { *Y++ = sqrtf(X[n]*X[n]+X[n+1]*X[n+1]); }
    }
    else if (L==N)
    {
        *Y = sqrtf(X[0]*X[0]+X[1]*X[1]);
        for (size_t n=1; n<N; n+=2) { *Y += sqrtf(X[n]*X[n]+X[n+1]*X[n+1]); }
    }
    else if (SH==1)
    {
        const size_t V = N/L;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t v=0; v<V; v++, Y++)
            {
                *Y = sqrtf(*X**X+*(X+1)**(X+1)); X += 2;
                for (size_t l=1; l<L; l++, X+=2) { *Y += sqrtf(*X**X+*(X+1)**(X+1)); }
            }
        }
        else
        {
            for (size_t v=0; v<V; v++, X+=2) { *Y++ = sqrtf(*X**X+*(X+1)**(X+1)); }
            Y -= V;
            for (size_t l=1; l<L; l++, Y-=V)
            {
                for (size_t v=0; v<V; v++, X+=2) { *Y++ += sqrtf(*X**X+*(X+1)**(X+1)); }
            }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y++)
            {
                *Y = sqrtf(*X**X+*(X+1)**(X+1)); X += 2;
                for (size_t s=1; s<S; s++, X+=2) { *Y += sqrtf(*X**X+*(X+1)**(X+1)); }
            }
        }
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, X-=2*(K*L-J), Y++)
            {
                *Y = sqrtf(*X**X+*(X+1)**(X+1)); X += 2*K;
                for (size_t l=1; l<L; l++, X+=2*K) { *Y += sqrtf(*X**X+*(X+1)**(X+1)); }
            }
        }
    }

    return 0;
}


int norm1_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in norm1_z: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (L==1)
    {
        for (size_t n=0; n<N; n+=2) { *Y++ = sqrt(X[n]*X[n]+X[n+1]*X[n+1]); }
    }
    else if (L==N)
    {
        *Y = sqrt(X[0]*X[0]+X[1]*X[1]);
        for (size_t n=1; n<N; n+=2) { *Y += sqrt(X[n]*X[n]+X[n+1]*X[n+1]); }
    }
    else if (SH==1)
    {
        const size_t V = N/L;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t v=0; v<V; v++, Y++)
            {
                *Y = sqrt(*X**X+*(X+1)**(X+1)); X += 2;
                for (size_t l=1; l<L; l++, X+=2) { *Y += sqrt(*X**X+*(X+1)**(X+1)); }
            }
        }
        else
        {
            for (size_t v=0; v<V; v++, X+=2) { *Y++ = sqrt(*X**X+*(X+1)**(X+1)); }
            Y -= V;
            for (size_t l=1; l<L; l++, Y-=V)
            {
                for (size_t v=0; v<V; v++, X+=2) { *Y++ += sqrt(*X**X+*(X+1)**(X+1)); }
            }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y++)
            {
                *Y = sqrt(*X**X+*(X+1)**(X+1)); X += 2;
                for (size_t s=1; s<S; s++, X+=2) { *Y += sqrt(*X**X+*(X+1)**(X+1)); }
            }
        }
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, X-=2*(K*L-J), Y++)
            {
                *Y = sqrt(*X**X+*(X+1)**(X+1)); X += 2*K;
                for (size_t l=1; l<L; l++, X+=2*K) { *Y += sqrt(*X**X+*(X+1)**(X+1)); }
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
