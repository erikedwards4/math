//Gets maximum of absolute values for each row or col of X according to dim.
//This is also the -Inf-norm (or min-norm) of each vector in X.

//For complex case, finds max absolute value and outputs the complex number.
//For complex case, this is |Xr|+|Xi|; see max for the usual sqrt(Xr*Xr+Xi*Xi).


#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int amin_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int amin_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int amin_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int amin_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int amin_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in amin_s: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    size_t i;

    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    if (N1==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        i = cblas_isamin((int)N,X,1);
        *Y = *(X+i);
    }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, Y++)
            {
                i = cblas_isamin((int)N1,X,1);
                X += i; *Y = *X; X += N1-i;
            }
        }
        else
        {
            float *mns;
            if (!(mns=(float *)calloc(N2,sizeof(float)))) { fprintf(stderr,"error in amin_s: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t n1=0; n1<N1; n1++, Y-=N2, mns-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X++, Y++, mns++)
                {
                    if (*X**X<*mns) { *Y = *X; *mns = *X**X; }
                }
            }
            free(mns);
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y++)
            {
                i = cblas_isamin((int)N1,X,1);
                X += i; *Y = *X; X += N1-i;
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
            for (size_t m=0; m<M; m++, Y++)
            {
                i = cblas_isamin((int)N1,X,(int)K);
                X += i*K; *Y = *X;
                X -= (int)((N1-i)*K)-(int)J;
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int amin_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in amin_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    size_t i;

    if (N1==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        i = cblas_idamin((int)N,X,1);
        *Y = *(X+i);
    }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, Y++)
            {
                i = cblas_idamin((int)N1,X,1);
                X += i; *Y = *X; X += N1-i;
            }
        }
        else
        {
            double *mns;
            if (!(mns=(double *)calloc(N2,sizeof(double)))) { fprintf(stderr,"error in amin_d: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t n1=0; n1<N1; n1++, Y-=N2, mns-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X++, Y++, mns++)
                {
                    if (*X**X<*mns) { *Y = *X; *mns = *X**X; }
                }
            }
            free(mns);
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y++)
            {
                i = cblas_idamin((int)N1,X,1);
                X += i; *Y = *X; X += N1-i;
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
            for (size_t m=0; m<M; m++, Y++)
            {
                i = cblas_idamin((int)N1,X,(int)K);
                X += i*K; *Y = *X;
                X -= (int)((N1-i)*K)-(int)J;
            }
        }
    }

    return 0;
}


int amin_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in amin_c: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    size_t i;

    if (N1==1) { cblas_ccopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        i = cblas_icamin((int)N,X,1);
        X += 2*i; *Y = *X; *(Y+1) = *(X+1);
    }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, Y+=2)
            {
                i = cblas_icamin((int)N1,X,1);
                X += 2*i; *Y = *X; *(Y+1) = *(X+1);
                X += 2*(N1-i);
            }
        }
        else
        {
            float *mns;
            if (!(mns=(float *)calloc(N2,sizeof(float)))) { fprintf(stderr,"error in amin_c: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t n1=0; n1<N1; n1++, Y-=2*N2, mns-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X+=2, Y+=2, mns++)
                {
                    if (sqrt(*X**X)+sqrt(*(X+1)**(X+1))<*mns)
                    {
                        *Y = *X; *(Y+1) = *(X+1);
                        *mns = sqrt(*X**X)+sqrt(*(X+1)**(X+1));
                    }
                }
            }
            free(mns);
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y+=2)
            {
                i = cblas_icamin((int)N1,X,1);
                X += 2*i; *Y = *X; *(Y+1) = *(X+1);
                X += 2*(N1-i);
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
            for (size_t m=0; m<M; m++, Y+=2)
            {
                i = cblas_icamin((int)N1,X,(int)K);
                X += 2*i*K; *Y = *X; *(Y+1) = *(X+1);
                X -= 2*((int)((N1-i)*K)-(int)J);
            }
        }
    }

    return 0;
}


int amin_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in amin_z: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    size_t i;

    if (N1==1) { cblas_zcopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        i = cblas_izamin((int)N,X,1);
        X += 2*i; *Y = *X; *(Y+1) = *(X+1);
    }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++, Y+=2)
            {
                i = cblas_izamin((int)N1,X,1);
                X += 2*i; *Y = *X; *(Y+1) = *(X+1);
                X += 2*(N1-i);
            }
        }
        else
        {
            double *mns;
            if (!(mns=(double *)calloc(N2,sizeof(double)))) { fprintf(stderr,"error in amin_z: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t n1=0; n1<N1; n1++, Y-=2*N2, mns-=N2)
            {
                for (size_t n2=0; n2<N2; n2++, X+=2, Y+=2, mns++)
                {
                    if (sqrt(*X**X)+sqrt(*(X+1)**(X+1))<*mns)
                    {
                        *Y = *X; *(Y+1) = *(X+1);
                        *mns = sqrt(*X**X)+sqrt(*(X+1)**(X+1));
                    }
                }
            }
            free(mns);
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y+=2)
            {
                i = cblas_izamin((int)N1,X,1);
                X += 2*i; *Y = *X; *(Y+1) = *(X+1);
                X += 2*(N1-i);
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
            for (size_t m=0; m<M; m++, Y+=2)
            {
                i = cblas_izamin((int)N1,X,(int)K);
                X += 2*i*K; *Y = *X; *(Y+1) = *(X+1);
                X -= 2*((int)((N1-i)*K)-(int)J);
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
