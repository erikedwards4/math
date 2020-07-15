//Gets L2-norm (root-sum-square) of each row or col of X according to dim.
//This is the Euclidean (L2) norm of each vector in X.
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

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X[n]*X[n]; }
    }
    else if (N1==N) { *Y = cblas_snrm2((int)N1,X,1); }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, X+=N1)
            {
                *Y++ = cblas_snrm2((int)N1,X,1);
            }
        }
        else
        {
            // for (size_t n2=0; n2<N2; n2++, X++)
            // {
            //     *Y++ = cblas_snrm2((int)N1,X,(int)N2);
            // }
            for (size_t n2=0; n2<N2; n2++, X++) { *Y++ = *X**X; }
            Y -= N2;
            for (size_t n1=1; n1<N1; n1++, Y-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X++) { *Y++ += *X**X; }
            }
            for (size_t n2=0; n2<N2; n2++, Y++) { *Y = sqrtf(*Y); }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, X+=N1)
            {
                *Y++ = cblas_snrm2((int)N1,X,1);
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
            for (size_t m=0; m<M; m++, X+=J)
            {
                *Y++ = cblas_snrm2((int)N1,X,(int)K);
            }
        }
    }

    return 0;
}


int norm2_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in norm2_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X[n]*X[n]; }
    }
    else if (N1==N) { *Y = cblas_dnrm2((int)N1,X,1); }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, X+=N1)
            {
                *Y++ = cblas_dnrm2((int)N1,X,1);
            }
        }
        else
        {
            for (size_t n2=0; n2<N2; n2++, X++) { *Y++ = *X**X; }
            Y -= N2;
            for (size_t n1=1; n1<N1; n1++, Y-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X++) { *Y++ += *X**X; }
            }
            for (size_t n2=0; n2<N2; n2++, Y++) { *Y = sqrt(*Y); }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, X+=N1)
            {
                *Y++ = cblas_dnrm2((int)N1,X,1);
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
            for (size_t m=0; m<M; m++, X+=J)
            {
                *Y++ = cblas_dnrm2((int)N1,X,(int)K);
            }
        }
    }

    return 0;
}


int norm2_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in norm2_c: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1)
    {
        for (size_t n=0; n<N; n+=2) { Y[n/2] = X[n]*X[n] + X[n+1]*X[n+1]; }
    }
    else if (N1==N) { *Y = cblas_scnrm2((int)N1,X,1); }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, X+=2*N1)
            {
                *Y++ = cblas_scnrm2((int)N1,X,1);
            }
        }
        else
        {
            for (size_t n2=0; n2<N2; n2++, X+=2) { *Y++ = *X**X + *(X+1)**(X+1); }
            Y -= N2;
            for (size_t n1=1; n1<N1; n1++, Y-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X+=2) { *Y++ += *X**X + *(X+1)**(X+1); }
            }
            for (size_t n2=0; n2<N2; n2++, Y++) { *Y = sqrtf(*Y); }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, X+=2*N1)
            {
                *Y++ = cblas_scnrm2((int)N1,X,1);
            }
        }
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t l=0; l<L; l++, X+=2*M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=2*J)
            {
                *Y++ = cblas_scnrm2((int)N1,X,(int)K);
            }
        }
    }

    return 0;
}


int norm2_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in norm2_z: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1)
    {
        for (size_t n=0; n<N; n+=2) { Y[n/2] = X[n]*X[n] + X[n+1]*X[n+1]; }
    }
    else if (N1==N) { *Y = cblas_dznrm2((int)N1,X,1); }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, X+=2*N1)
            {
                *Y++ = cblas_dznrm2((int)N1,X,1);
            }
        }
        else
        {
            for (size_t n2=0; n2<N2; n2++, X+=2) { *Y++ = *X**X + *(X+1)**(X+1); }
            Y -= N2;
            for (size_t n1=1; n1<N1; n1++, Y-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X+=2) { *Y++ += *X**X + *(X+1)**(X+1); }
            }
            for (size_t n2=0; n2<N2; n2++, Y++) { *Y = sqrt(*Y); }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, X+=2*N1)
            {
                *Y++ = cblas_dznrm2((int)N1,X,1);
            }
        }
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t l=0; l<L; l++, X+=2*M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=2*J)
            {
                *Y++ = cblas_dznrm2((int)N1,X,(int)K);
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
