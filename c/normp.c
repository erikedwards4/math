//Gets the p-norm along dim of X.
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

int normp_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor, const float p);
int normp_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor, const double p);
int normp_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor, const float p);
int normp_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor, const double p);


int normp_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor, const float p)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ip = 1.0f / p;

    if (N1==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = fabsf(X[n]); }
    }
    else if (N1==N)
    {
        *Y = powf(fabsf(X[0]),p);
        for (size_t n=1; n<N; n++) { *Y += powf(fabsf(X[n]),p); }
        *Y = powf(*Y,ip);
    }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, Y++)
            {
                *Y = powf(fabsf(*X),p); X++;
                for (size_t n1=1; n1<N1; n1++, X++) { *Y += powf(fabsf(*X),p); }
                *Y = powf(*Y,ip);
            }
        }
        else
        {
            for (size_t n2=0; n2<N2; n2++, X++) { *Y++ = powf(fabsf(*X),p); }
            Y -= N2;
            for (size_t n1=1; n1<N1; n1++, Y-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X++) { *Y++ += powf(fabsf(*X),p); }
            }
            for (size_t n2=0; n2<N2; n2++, Y++) { *Y = powf(*Y,ip); }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y++)
            {
                *Y = powf(fabsf(*X),p); X++;
                for (size_t s=1; s<S; s++, X++) { *Y += powf(fabsf(*X),p); }
                *Y = powf(*Y,ip);
            }
        }
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X-=K*N1-J, Y++)
            {
                *Y = powf(fabsf(*X),p); X += K;
                for (size_t n1=1; n1<N1; n1++, X+=K) { *Y += powf(fabsf(*X),p); }
                *Y = powf(*Y,ip);
            }
        }
    }

    return 0;
}


int normp_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor, const double p)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ip = 1.0 / p;

    if (N1==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = fabs(X[n]); }
    }
    else if (N1==N)
    {
        *Y = pow(fabs(X[0]),p);
        for (size_t n=1; n<N; n++) { *Y += pow(fabs(X[n]),p); }
        *Y = pow(*Y,ip);
    }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, Y++)
            {
                *Y = pow(fabs(*X),p); X++;
                for (size_t n1=1; n1<N1; n1++, X++) { *Y += pow(fabs(*X),p); }
                *Y = pow(*Y,ip);
            }
        }
        else
        {
            for (size_t n2=0; n2<N2; n2++, X++) { *Y++ = pow(fabs(*X),p); }
            Y -= N2;
            for (size_t n1=1; n1<N1; n1++, Y-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X++) { *Y++ += pow(fabs(*X),p); }
            }
            for (size_t n2=0; n2<N2; n2++, Y++) { *Y = pow(*Y,ip); }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y++)
            {
                *Y = pow(fabs(*X),p); X++;
                for (size_t s=1; s<S; s++, X++) { *Y += pow(fabs(*X),p); }
                *Y = pow(*Y,ip);
            }
        }
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X-=K*N1-J, Y++)
            {
                *Y = pow(fabs(*X),p); X += K;
                for (size_t n1=1; n1<N1; n1++, X+=K) { *Y += pow(fabs(*X),p); }
                *Y = pow(*Y,ip);
            }
        }
    }

    return 0;
}


int normp_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor, const float p)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ip = 1.0f / p;

    if (N1==1)
    {
        for (size_t n=0; n<N; n+=2) { *Y++ = sqrtf(X[n]*X[n]+X[n+1]*X[n+1]); }
    }
    else if (N1==N)
    {
        *Y = powf(sqrtf(X[0]*X[0]+X[1]*X[1]),p);
        for (size_t n=1; n<N; n+=2) { *Y += powf(sqrtf(X[n]*X[n]+X[n+1]*X[n+1]),p); }
        *Y = powf(*Y,ip);
    }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, Y++)
            {
                *Y = powf(sqrtf(*X**X+*(X+1)**(X+1)),p); X += 2;
                for (size_t n1=1; n1<N1; n1++, X+=2) { *Y += powf(sqrtf(*X**X+*(X+1)**(X+1)),p); }
                *Y = powf(*Y,ip);
            }
        }
        else
        {
            for (size_t n2=0; n2<N2; n2++, X+=2) { *Y++ = powf(sqrtf(*X**X+*(X+1)**(X+1)),p); }
            Y -= N2;
            for (size_t n1=1; n1<N1; n1++, Y-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X+=2) { *Y++ += powf(sqrtf(*X**X+*(X+1)**(X+1)),p); }
            }
            for (size_t n2=0; n2<N2; n2++, Y++) { *Y = powf(*Y,ip); }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y++)
            {
                *Y = powf(sqrtf(*X**X+*(X+1)**(X+1)),p); X += 2;
                for (size_t s=1; s<S; s++, X+=2) { *Y += powf(sqrtf(*X**X+*(X+1)**(X+1)),p); }
                *Y = powf(*Y,ip);
            }
        }
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X-=2*(K*N1-J), Y++)
            {
                *Y = powf(sqrtf(*X**X+*(X+1)**(X+1)),p); X += 2*K;
                for (size_t n1=1; n1<N1; n1++, X+=2*K) { *Y += powf(sqrtf(*X**X+*(X+1)**(X+1)),p); }
                *Y = powf(*Y,ip);
            }
        }
    }

    return 0;
}


int normp_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor, const double p)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ip = 1.0 / p;

    if (N1==1)
    {
        for (size_t n=0; n<N; n+=2) { *Y++ = sqrt(X[n]*X[n]+X[n+1]*X[n+1]); }
    }
    else if (N1==N)
    {
        *Y = pow(sqrt(X[0]*X[0]+X[1]*X[1]),p);
        for (size_t n=1; n<N; n+=2) { *Y += pow(sqrt(X[n]*X[n]+X[n+1]*X[n+1]),p); }
        *Y = pow(*Y,ip);
    }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, Y++)
            {
                *Y = pow(sqrt(*X**X+*(X+1)**(X+1)),p); X += 2;
                for (size_t n1=1; n1<N1; n1++, X+=2) { *Y += pow(sqrt(*X**X+*(X+1)**(X+1)),p); }
                *Y = pow(*Y,ip);
            }
        }
        else
        {
            for (size_t n2=0; n2<N2; n2++, X+=2) { *Y++ = pow(sqrt(*X**X+*(X+1)**(X+1)),p); }
            Y -= N2;
            for (size_t n1=1; n1<N1; n1++, Y-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X+=2) { *Y++ += pow(sqrt(*X**X+*(X+1)**(X+1)),p); }
            }
            for (size_t n2=0; n2<N2; n2++, Y++) { *Y = pow(*Y,ip); }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y++)
            {
                *Y = pow(sqrt(*X**X+*(X+1)**(X+1)),p); X += 2;
                for (size_t s=1; s<S; s++, X+=2) { *Y += pow(sqrt(*X**X+*(X+1)**(X+1)),p); }
                *Y = pow(*Y,ip);
            }
        }
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X-=2*(K*N1-J), Y++)
            {
                *Y = pow(sqrt(*X**X+*(X+1)**(X+1)),p); X += 2*K;
                for (size_t n1=1; n1<N1; n1++, X+=2*K) { *Y += pow(sqrt(*X**X+*(X+1)**(X+1)),p); }
                *Y = pow(*Y,ip);
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
