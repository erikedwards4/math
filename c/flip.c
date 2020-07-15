//Flips X along dim (reverses order of elements).
//For d=0, this is flipud. For d=1, this is fliplr.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int flip_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int flip_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int flip_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int flip_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);

int flip_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int flip_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int flip_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int flip_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int flip_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_s: dim must be in [0 3]\n"); return 1; }

    const size_t L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const size_t M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1))
    {
        cblas_scopy((int)(R*C*S*H),X,1,Y,1);
    }
    else if (K==1)
    {
        for (size_t l=0, n1=0; l<L; l++)
        {
            for (size_t m=0, n2=n1+M-1; m<M; m++, n1++, n2--) { Y[n2] = X[n1]; }
        }
    }
    else
    {
        for (size_t l=0, n1=0; l<L; l++)
        {
            for (size_t m=0, n2=n1+M-K; m<M/K; m++, n1+=K, n2-=K)
            {
                cblas_scopy((int)K,&X[n1],1,&Y[n2],1);
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int flip_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_d: dim must be in [0 3]\n"); return 1; }

    const size_t L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const size_t M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1))
    {
        cblas_dcopy((int)(R*C*S*H),X,1,Y,1);
    }
    else if (K==1)
    {
        for (size_t l=0, n1=0; l<L; l++)
        {
            for (size_t m=0, n2=n1+M-1; m<M; m++, n1++, n2--) { Y[n2] = X[n1]; }
        }
    }
    else
    {
        for (size_t l=0, n1=0; l<L; l++)
        {
            for (size_t m=0, n2=n1+M-K; m<M/K; m++, n1+=K, n2-=K)
            {
                cblas_dcopy((int)K,&X[n1],1,&Y[n2],1);
            }
        }
    }

    return 0;
}


int flip_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_c: dim must be in [0 3]\n"); return 1; }

    const size_t L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const size_t M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1))
    {
        cblas_ccopy((int)(R*C*S*H),X,1,Y,1);
    }
    else if (K==1)
    {
        for (size_t l=0, n1=0; l<L; l++)
        {
            for (size_t m=0, n2=n1+2*M-2; m<M; m++, n1+=2, n2-=2) { Y[n2] = X[n1]; Y[n2+1] = X[n1+1]; }
        }
    }
    else
    {
        for (size_t l=0, n1=0; l<L; l++)
        {
            for (size_t m=0, n2=n1+2*(M-K); m<M/K; m++, n1+=2*K, n2-=2*K)
            {
                cblas_ccopy((int)K,&X[n1],1,&Y[n2],1);
            }
        }
    }

    return 0;
}


int flip_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_z: dim must be in [0 3]\n"); return 1; }

    const size_t L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const size_t M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    
    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1))
    {
        cblas_zcopy((int)(R*C*S*H),X,1,Y,1);
    }
    else if (K==1)
    {
        for (size_t l=0, n1=0; l<L; l++)
        {
            for (size_t m=0, n2=n1+2*M-2; m<M; m++, n1+=2, n2-=2) { Y[n2] = X[n1]; Y[n2+1] = X[n1+1]; }
        }
    }
    else
    {
        for (size_t l=0, n1=0; l<L; l++)
        {
            for (size_t m=0, n2=n1+2*(M-K); m<M/K; m++, n1+=2*K, n2-=2*K)
            {
                cblas_zcopy((int)K,&X[n1],1,&Y[n2],1);
            }
        }
    }

    return 0;
}


int flip_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const size_t M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1)) {}
    else if (K==1)
    {
        float x1;
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n1=l*M, n2=n1+M-1; m<M/2; m++, n1++, n2--)
            {
                x1 = X[n1]; X[n1] = X[n2]; X[n2] = x1;
            }
        }
    }
    else
    {
        float *X1;
        if (!(X1=(float *)malloc(K*sizeof(float)))) { fprintf(stderr,"error in flip_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n1=l*M, n2=n1+M-K; m<M/K/2; m++, n1+=K, n2-=K)
            {
                //for (n=0; n<K; n++, n1++, n2--) { x1 = X[n1]; X[n1] = X[n2]; X[n2] = x1; }
                cblas_scopy((int)K,&X[n1],1,X1,1);
                cblas_scopy((int)K,&X[n2],1,&X[n1],1);
                cblas_scopy((int)K,X1,1,&X[n2],1);
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int flip_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const size_t M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1)) {}
    else if (K==1)
    {
        double x1;
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n1=l*M, n2=n1+M-1; m<M/2; m++, n1++, n2--)
            {
                x1 = X[n1]; X[n1] = X[n2]; X[n2] = x1;
            }
        }
    }
    else
    {
        double *X1;
        if (!(X1=(double *)malloc(K*sizeof(double)))) { fprintf(stderr,"error in flip_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n1=l*M, n2=n1+M-K; m<M/K/2; m++, n1+=K, n2-=K)
            {
                cblas_dcopy((int)K,&X[n1],1,X1,1);
                cblas_dcopy((int)K,&X[n2],1,&X[n1],1);
                cblas_dcopy((int)K,X1,1,&X[n2],1);
            }
        }
    }

    return 0;
}


int flip_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_inplace_c: dim must be in [0 3]\n"); return 1; }

    const size_t L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const size_t M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1)) {}
    else if (K==1)
    {
        float x1;
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n1=2*l*M, n2=n1+2*M-2; m<M/2; m++, n1+=2, n2-=2)
            {
                x1 = X[n1]; X[n1] = X[n2]; X[n2] = x1;
                x1 = X[n1+1]; X[n1+1] = X[n2+1]; X[n2+1] = x1;
            }
        }
    }
    else
    {
        float *X1;
        if (!(X1=(float *)malloc(K*sizeof(float)))) { fprintf(stderr,"error in flip_inplace_c: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n1=2*l*M, n2=n1+2*(M-K); m<M/K/2; m++, n1+=2*K, n2-=2*K)
            {
                cblas_ccopy((int)K,&X[n1],1,X1,1);
                cblas_ccopy((int)K,&X[n2],1,&X[n1],1);
                cblas_ccopy((int)K,X1,1,&X[n2],1);
            }
        }
    }

    return 0;
}


int flip_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_inplace_z: dim must be in [0 3]\n"); return 1; }

    const size_t L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const size_t M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1)) {}
    else if (K==1)
    {
        double x1;
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n1=2*l*M, n2=n1+2*M-2; m<M/2; m++, n1+=2, n2-=2)
            {
                x1 = X[n1]; X[n1] = X[n2]; X[n2] = x1;
                x1 = X[n1+1]; X[n1+1] = X[n2+1]; X[n2+1] = x1;
            }
        }
    }
    else
    {
        double *X1;
        if (!(X1=(double *)malloc(K*sizeof(double)))) { fprintf(stderr,"error in flip_inplace_z: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n1=2*l*M, n2=n1+2*(M-K); m<M/K/2; m++, n1+=2*K, n2-=2*K)
            {
                cblas_zcopy((int)K,&X[n1],1,X1,1);
                cblas_zcopy((int)K,&X[n2],1,&X[n1],1);
                cblas_zcopy((int)K,X1,1,&X[n2],1);
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
