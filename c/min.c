//Gets minimum of values for each row or col of X according to dim.
//For complex case, finds min absolute value and outputs the complex number.
//For complex case, this is the proper absolute value; see amin for the other one.

//The in-place version was always slower, since requires rewind to start of X,
//so those are no longer included for any of the Stats functions.
//Instead, in-place will mean for Stats that X is allowed to be modified.

//Unfortunately, compiler doesn't find cblas_i?amin (but finds cblas_i?amax),
//so I do a direct solution for complex-valued case for now.


#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int min_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int min_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int min_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int min_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int min_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in min_s: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        *Y = *X;
        //for (size_t n=1; n<N; n++, X++) { if (*X<*Y) { *Y = *X; } }
        for (size_t n=1; n<N; n++) { if (X[n]<*Y) { *Y = X[n]; } }
    }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, Y++)
            {
                *Y = *X++;
                for (size_t n1=1; n1<N1; n1++, X++) { if (*X<*Y) { *Y = *X; } }
            }
        }
        else
        {
            for (size_t n1=0; n1<N1; n1++, Y-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X++, Y++)
                {
                    if (n1==0 || *X<*Y) { *Y = *X; }
                }
            }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y++)
            {
                *Y = *X++;
                for (size_t n1=1; n1<N1; n1++, X++) { if (*X<*Y) { *Y = *X; } }
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
                *Y = *X; X += K;
                for (size_t n1=1; n1<N1; n1++, X+=K) { if (*X<*Y) { *Y = *X; } }
            }
        }
    }

    return 0;
}


int min_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in min_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        *Y = *X;
        for (size_t n=1; n<N; n++) { if (X[n]<*Y) { *Y = X[n]; } }
    }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, Y++)
            {
                *Y = *X++;
                for (size_t n1=1; n1<N1; n1++, X++) { if (*X<*Y) { *Y = *X; } }
            }
        }
        else
        {
            for (size_t n1=0; n1<N1; n1++, Y-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X++, Y++)
                {
                    if (n1==0 || *X<*Y) { *Y = *X; }
                }
            }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y++)
            {
                *Y = *X++;
                for (size_t n1=1; n1<N1; n1++, X++) { if (*X<*Y) { *Y = *X; } }
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
                *Y = *X; X += K;
                for (size_t n1=1; n1<N1; n1++, X+=K) { if (*X<*Y) { *Y = *X; } }
            }
        }
    }
    
    return 0;
}


int min_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in min_c: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float xx, mn;

    if (N1==1) { cblas_ccopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        mn = *X**X + *(X+1)**(X+1);
        *Y = *X; *(Y+1) = *(X+1); X += 2;
        for (size_t n=1; n<N; n++, X+=2)
        {
            xx = *X**X + *(X+1)**(X+1);
            if (xx<mn) { mn = xx; *Y = *X; *(Y+1) = *(X+1); }
        }
    }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, Y+=2)
            {
                mn = *X**X + *(X+1)**(X+1);
                *Y = *X; *(Y+1) = *(X+1); X += 2;
                for (size_t n1=1; n1<N1; n1++, X+=2)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (xx<mn) { mn = xx; *Y = *X; *(Y+1) = *(X+1); }
                }
            }
        }
        else
        {
            float *mns;
            if (!(mns=(float *)malloc(N2*sizeof(float)))) { fprintf(stderr,"error in min_c: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t n1=0; n1<N1; n1++, Y-=2*N2, mns-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X+=2, Y+=2, mns++)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (n1==0 || xx<*mns) { *mns = xx; *Y = *X; *(Y+1) = *(X+1); }
                }
            }
            free(mns);
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
            for (size_t m=0; m<M; m++, X-=2*(K*N1-J), Y+=2)
            {
                mn = *X**X + *(X+1)**(X+1);
                *Y = *X; *(Y+1) = *(X+1); X += 2*K;
                for (size_t n1=1; n1<N1; n1++, X+=2*K)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (xx<mn) { mn = xx; *Y = *X; *(Y+1) = *(X+1); }
                }
            }
        }
    }

    return 0;
}


int min_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in min_z: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double xx, mn;

    if (N1==1) { cblas_zcopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        mn = *X**X + *(X+1)**(X+1);
        *Y = *X; *(Y+1) = *(X+1); X += 2;
        for (size_t n=1; n<N; n++, X+=2)
        {
            xx = *X**X + *(X+1)**(X+1);
            if (xx<mn) { mn = xx; *Y = *X; *(Y+1) = *(X+1); }
        }
    }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, Y+=2)
            {
                mn = *X**X + *(X+1)**(X+1);
                *Y = *X; *(Y+1) = *(X+1); X += 2;
                for (size_t n1=1; n1<N1; n1++, X+=2)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (xx<mn) { mn = xx; *Y = *X; *(Y+1) = *(X+1); }
                }
            }
        }
        else
        {
            double *mns;
            if (!(mns=(double *)malloc(N2*sizeof(double)))) { fprintf(stderr,"error in min_z: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t n1=0; n1<N1; n1++, Y-=2*N2, mns-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X+=2, Y+=2, mns++)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (n1==0 || xx<*mns) { *mns = xx; *Y = *X; *(Y+1) = *(X+1); }
                }
            }
            free(mns);
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
            for (size_t m=0; m<M; m++, X-=2*(K*N1-J), Y+=2)
            {
                mn = *X**X + *(X+1)**(X+1);
                *Y = *X; *(Y+1) = *(X+1);
                for (size_t n1=2; n1<2*N1; n1+=2)
                {
                    xx = *(X+n1*K)**(X+n1*K) + *(X+n1*K+1)**(X+n1*K+1);
                    if (xx<mn) { mn = xx; *Y = *(X+n1*K); *(Y+1) = *(X+n1*K+1); }
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
