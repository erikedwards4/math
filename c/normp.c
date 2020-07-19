//Gets the p-norm for each vector in X along dim.
//This is the Lp norm of each vector in X.
//For each vector, y = sum(|x|^p)^1/p.
//For complex case, output is real.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

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
    else if (L==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        *Y = powf(fabsf(X[0]),p);
        for (size_t l=1; l<L; l++) { *Y += powf(fabsf(X[l]),p); }
        *Y = powf(*Y,ip);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, Y++)
            {
                *Y = powf(fabsf(*X),p); X++;
                for (size_t l=1; l<L; l++, X++) { *Y += powf(fabsf(*X),p); }
                *Y = powf(*Y,ip);
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; v++, X++) { *Y++ = powf(fabsf(*X),p); }
            Y -= V;
            for (size_t l=1; l<L; l++, Y-=V)
            {
                for (size_t v=0; v<V; v++, X++) { *Y++ += powf(fabsf(*X),p); }
            }
            for (size_t v=0; v<V; v++, Y++) { *Y = powf(*Y,ip); }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, X-=K*L-1, Y++)
                {
                    *Y = powf(fabsf(*X),p); X += K;
                    for (size_t l=1; l<L; l++, X+=K) { *Y += powf(fabsf(*X),p); }
                    *Y = powf(*Y,ip);
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
    else if (L==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        *Y = pow(fabs(X[0]),p);
        for (size_t l=1; l<L; l++) { *Y += pow(fabs(X[l]),p); }
        *Y = pow(*Y,ip);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, Y++)
            {
                *Y = pow(fabs(*X),p); X++;
                for (size_t l=1; l<L; l++, X++) { *Y += pow(fabs(*X),p); }
                *Y = pow(*Y,ip);
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; v++, X++) { *Y++ = pow(fabs(*X),p); }
            Y -= V;
            for (size_t l=1; l<L; l++, Y-=V)
            {
                for (size_t v=0; v<V; v++, X++) { *Y++ += pow(fabs(*X),p); }
            }
            for (size_t v=0; v<V; v++, Y++) { *Y = pow(*Y,ip); }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, X-=K*L-1, Y++)
                {
                    *Y = pow(fabs(*X),p); X += K;
                    for (size_t l=1; l<L; l++, X+=K) { *Y += pow(fabs(*X),p); }
                    *Y = pow(*Y,ip);
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
        for (size_t n=0; n<N; n+=2) { *Y++ = sqrtf(X[n]*X[n]+X[n+1]*X[n+1]); }
    }
    else if (L==N)
    {
        *Y = powf(sqrtf(X[0]*X[0]+X[1]*X[1]),p);
        for (size_t l=1; l<L; l+=2) { *Y += powf(sqrtf(X[l]*X[l]+X[l+1]*X[l+1]),p); }
        *Y = powf(*Y,ip);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, Y++)
            {
                *Y = powf(sqrtf(*X**X+*(X+1)**(X+1)),p); X += 2;
                for (size_t l=1; l<L; l++, X+=2) { *Y += powf(sqrtf(*X**X+*(X+1)**(X+1)),p); }
                *Y = powf(*Y,ip);
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; v++, X+=2) { *Y++ = powf(sqrtf(*X**X+*(X+1)**(X+1)),p); }
            Y -= V;
            for (size_t l=1; l<L; l++, Y-=V)
            {
                for (size_t v=0; v<V; v++, X+=2) { *Y++ += powf(sqrtf(*X**X+*(X+1)**(X+1)),p); }
            }
            for (size_t v=0; v<V; v++, Y++) { *Y = powf(*Y,ip); }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; b++, X-=2*K*L-2, Y++)
                {
                    *Y = powf(sqrtf(*X**X+*(X+1)**(X+1)),p); X += 2*K;
                    for (size_t l=1; l<L; l++, X+=2*K) { *Y += powf(sqrtf(*X**X+*(X+1)**(X+1)),p); }
                    *Y = powf(*Y,ip);
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
        for (size_t n=0; n<N; n+=2) { *Y++ = sqrt(X[n]*X[n]+X[n+1]*X[n+1]); }
    }
    else if (L==N)
    {
        *Y = pow(sqrt(X[0]*X[0]+X[1]*X[1]),p);
        for (size_t l=1; l<L; l+=2) { *Y += pow(sqrt(X[l]*X[l]+X[l+1]*X[l+1]),p); }
        *Y = pow(*Y,ip);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, Y++)
            {
                *Y = pow(sqrt(*X**X+*(X+1)**(X+1)),p); X += 2;
                for (size_t l=1; l<L; l++, X+=2) { *Y += pow(sqrt(*X**X+*(X+1)**(X+1)),p); }
                *Y = pow(*Y,ip);
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; v++, X+=2) { *Y++ = pow(sqrt(*X**X+*(X+1)**(X+1)),p); }
            Y -= V;
            for (size_t l=1; l<L; l++, Y-=V)
            {
                for (size_t v=0; v<V; v++, X+=2) { *Y++ += pow(sqrt(*X**X+*(X+1)**(X+1)),p); }
            }
            for (size_t v=0; v<V; v++, Y++) { *Y = pow(*Y,ip); }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; b++, X-=2*K*L-2, Y++)
                {
                    *Y = pow(sqrt(*X**X+*(X+1)**(X+1)),p); X += 2*K;
                    for (size_t l=1; l<L; l++, X+=2*K) { *Y += pow(sqrt(*X**X+*(X+1)**(X+1)),p); }
                    *Y = pow(*Y,ip);
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
