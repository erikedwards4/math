//Non-hermitian transpose of matrix X.
//This has in-place and not-in-place versions.

//For the inplace R!=C case, there may be a faster way to do this (with much more code).
//However, it is already ~same speed as Armadillo inplace_strans.

//For the inplace R==C case, using cblas_?swap is slower for small N,
//slightly slower for medium N, and slightly faster for large N (e.g. N=10000).
//But this gain of less than 1% at large N does not seem worth it, so removed.
//However, for complex case, the cblas_?swap is definitely faster for N>3000.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int transpose_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor);
int transpose_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor);
int transpose_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor);
int transpose_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor);

int transpose_inplace_s (float *X, const size_t R, const size_t C, const char iscolmajor);
int transpose_inplace_d (double *X, const size_t R, const size_t C, const char iscolmajor);
int transpose_inplace_c (float *X, const size_t R, const size_t C, const char iscolmajor);
int transpose_inplace_z (double *X, const size_t R, const size_t C, const char iscolmajor);


int transpose_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor)
{
    const size_t N = R*C;

    if (R==1 || C==1)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (iscolmajor)
    {
        for (size_t c=0u; c<C; ++c, Y-=N-1)
        {
            for (size_t r=0u; r<R; ++r, ++X, Y+=C) { *Y = *X; }
        }
    }
    else
    {
        for (size_t r=0u; r<R; ++r, Y-=N-1)
        {
            for (size_t c=0u; c<C; ++c, ++X, Y+=R) { *Y = *X; }
        }
    }

    return 0;
}


int transpose_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor)
{
    const size_t N = R*C;

    if (R==1 || C==1)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (iscolmajor)
    {
        for (size_t c=0u; c<C; ++c, Y-=N-1)
        {
            for (size_t r=0u; r<R; ++r, ++X, Y+=C) { *Y = *X; }
        }
    }
    else
    {
        for (size_t r=0u; r<R; ++r, Y-=N-1)
        {
            for (size_t c=0u; c<C; ++c, ++X, Y+=R) { *Y = *X; }
        }
    }

    return 0;
}


int transpose_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor)
{
    const size_t N = R*C;

    if (R==1 || C==1)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
    }
    else if (iscolmajor)
    {
        for (size_t c=0u; c<C; ++c, Y-=2*N-2)
        {
            for (size_t r=0u; r<R; ++r, ++X, Y+=2*C-1) { *Y = *X; *++Y = *++X; }
        }
    }
    else
    {
        for (size_t r=0u; r<R; ++r, Y-=2*N-2)
        {
            for (size_t c=0u; c<C; ++c, ++X, Y+=2*R-1) { *Y = *X; *++Y = *++X; }
        }
    }

    return 0;
}


int transpose_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor)
{
    const size_t N = R*C;

    if (R==1 || C==1)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
    }
    else if (iscolmajor)
    {
        for (size_t c=0u; c<C; ++c, Y-=2*N-2)
        {
            for (size_t r=0u; r<R; ++r, ++X, Y+=2*C-1) { *Y = *X; *++Y = *++X; }
        }
    }
    else
    {
        for (size_t r=0u; r<R; ++r, Y-=2*N-2)
        {
            for (size_t c=0u; c<C; ++c, ++X, Y+=2*R-1) { *Y = *X; *++Y = *++X; }
        }
    }

    return 0;
}


int transpose_inplace_s (float *X, const size_t R, const size_t C, const char iscolmajor)
{
    if (R==1 || C==1) {}
    else if (R==C)
    {
        size_t k;
        float x1;
        ++X;
        for (size_t c=0u; c<C-1; ++c, X+=c+1)
        {
            for (size_t r=0u; r<R-c-1; ++r, ++X)
            {
                k = (r+1)*R-r-1;
                x1 = *X; *X = *(X+k); *(X+k) = x1;
            }
        }
        // for (size_t c=0u; c<C-1; ++c, X+=R+1) //slightly faster for N>1e5
        // {
        //     cblas_sswap((int)(R-c-1),X,1,&X[R-1],(int)C);
        // }
    }
    else
    {
        const size_t N = R*C;
        float *Xt;
        if (!(Xt=(float *)malloc(N*sizeof(float)))) { fprintf(stderr,"error in transpose_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t n=0u; n<N; ++n, ++X, ++Xt) { *Xt = *X; }
        X -= N; Xt -= N;
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c, X-=N-1)
            {
                for (size_t r=0u; r<R; ++r, ++Xt, X+=C) { *X = *Xt; }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r, X-=N-1)
            {
                for (size_t c=0u; c<C; ++c, ++Xt, X+=R) { *X = *Xt; }
            }
        }
        Xt -= N;
        free(Xt);
    }

    return 0;
}


int transpose_inplace_d (double *X, const size_t R, const size_t C, const char iscolmajor)
{
    if (R==1 || C==1) {}
    else if (R==C)
    {
        size_t k;
        double x1;
        ++X;
        for (size_t c=0u; c<C-1; ++c, X+=c+1)
        {
            for (size_t r=0u; r<R-c-1; ++r, ++X)
            {
                k = (r+1)*R-r-1;
                x1 = *X; *X = *(X+k); *(X+k) = x1;
            }
        }
    }
    else
    {
        const size_t N = R*C;
        double *Xt;
        if (!(Xt=(double *)malloc(N*sizeof(double)))) { fprintf(stderr,"error in transpose_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t n=0u; n<N; ++n, ++X, ++Xt) { *Xt = *X; }
        X -= N; Xt -= N;
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c, X-=N-1)
            {
                for (size_t r=0u; r<R; ++r, ++Xt, X+=C) { *X = *Xt; }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r, X-=N-1)
            {
                for (size_t c=0u; c<C; ++c, ++Xt, X+=R) { *X = *Xt; }
            }
        }
        Xt -= N;
        free(Xt);
    }

    return 0;
}


int transpose_inplace_c (float *X, const size_t R, const size_t C, const char iscolmajor)
{
    if (R==1 || C==1) {}
    else if (R==C)
    {
        if (R<2500)
        {
            size_t k;
            float x1;
            X += 2;
            for (size_t c=0u; c<C-1; ++c, X+=2*c+2)
            {
                for (size_t r=0u; r<R-c-1; ++r, ++X)
                {
                    k = 2*((r+1)*R-r-1);
                    x1 = *X; *X = *(X+k); *(X+k) = x1;
                    ++X;
                    x1 = *X; *X = *(X+k); *(X+k) = x1;
                }
            }
        }
        else
        {
            for (size_t c=0, n=2; c<C-1; ++c, n+=2*R+2)
            {
                cblas_cswap((int)(R-c-1),&X[n],1,&X[2*((c+1)*R+c)],(int)C);
            }
        }
    }
    else
    {
        const size_t N = R*C;
        float *Xt;
        if (!(Xt=(float *)malloc(2*N*sizeof(float)))) { fprintf(stderr,"error in transpose_inplace_c: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t n=0u; n<2*N; ++n, ++X, ++Xt) { *Xt = *X; }
        X -= 2*N; Xt -= 2*N;
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c, X-=2*N-2)
            {
                for (size_t r=0u; r<R; ++r, ++Xt, X+=2*C-1) { *X = *Xt; *++X = *++Xt; }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r, X-=2*N-2)
            {
                for (size_t c=0u; c<C; ++c, ++Xt, X+=2*R-1) { *X = *Xt; *++X = *++Xt; }
            }
        }
        Xt -= 2*N;
        free(Xt);
    }

    return 0;
}


int transpose_inplace_z (double *X, const size_t R, const size_t C, const char iscolmajor)
{
    if (R==1 || C==1) {}
    else if (R==C)
    {
        if (R<2500)
        {
            size_t k;
            double x1;
            X += 2;
            for (size_t c=0u; c<C-1; ++c, X+=2*c+2)
            {
                for (size_t r=0u; r<R-c-1; ++r, ++X)
                {
                    k = 2*((r+1)*R-r-1);
                    x1 = *X; *X = *(X+k); *(X+k) = x1;
                    ++X;
                    x1 = *X; *X = *(X+k); *(X+k) = x1;
                }
            }
        }
        else
        {
            for (size_t c=0, n=2; c<C-1; ++c, n+=2*R+2)
            {
                cblas_zswap((int)(R-c-1),&X[n],1,&X[2*((c+1)*R+c)],(int)C);
            }
        }
    }
    else
    {
        const size_t N = R*C;
        double *Xt;
        if (!(Xt=(double *)malloc(2*N*sizeof(double)))) { fprintf(stderr,"error in transpose_inplace_z: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t n=0u; n<2*N; ++n, ++X, ++Xt) { *Xt = *X; }
        X -= 2*N; Xt -= 2*N;
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c, X-=2*N-2)
            {
                for (size_t r=0u; r<R; ++r, ++Xt, X+=2*C-1) { *X = *Xt; *++X = *++Xt; }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r, X-=2*N-2)
            {
                for (size_t c=0u; c<C; ++c, ++Xt, X+=2*R-1) { *X = *Xt; *++X = *++Xt; }
            }
        }
        Xt -= 2*N;
        free(Xt);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
