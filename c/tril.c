//Zeros all elements above the kth diagonal of matrix X.
//This has in-place and not-in-place versions.

//LAPACKE_?lacpy was surprisingly a bit slower even for square, k==0 case.

#include <stdio.h>
//#include <lapacke.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int tril_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);
int tril_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);
int tril_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);
int tril_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);

int tril_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);
int tril_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);
int tril_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);
int tril_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k);


int tril_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_s: k must be in [1-R C-1]\n"); return 1; }

    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    // if (k==0 && R==C && S*H==1u)
    // {
    //     const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
    //     const int lda = (iscolmajor) ? R : C;
    //     if (LAPACKE_slacpy (Ord,'L',(int)R,(int)C,X,lda,Y,lda))
    //     { fprintf(stderr,"error in tril_s: problem with LAPACKE function\n"); return 1; }
    // }
    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0u;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0u;                           //number of all-X cols
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s)
            {
                for (size_t n=R*CX; n>0u; --n, ++X, ++Y) { *Y = *X; }
                for (size_t c=CX; c<C-C0; ++c)
                {
                    for (int n=(int)c-k; n>0; --n, ++X, ++Y) { *Y = 0.0f; }
                    for (int n=(int)R-(int)c+k; n>0; --n, ++X, ++Y) { *Y = *X; }
                }
                for (size_t n=R*C0; n>0u; --n, ++X, ++Y) { *Y = 0.0f; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0u;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0u;    //number of all-X rows
        for (size_t n=R0*C*SH; n>0u; --n, ++X, ++Y) { *Y = 0.0f; }
        for (size_t r=R0; r<R-RX; ++r)
        {
            for (int n=(int)SH*((int)r+k+1); n>0; --n, ++X, ++Y) { *Y = *X; }
            for (int n=(int)SH*((int)C-(int)r-k-1); n>0; --n, ++X, ++Y) { *Y = 0.0f; }
        }
        for (size_t n=RX*C*SH; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }

    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int tril_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_d: k must be in [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0u;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0u;                           //number of all-X cols
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s)
            {
                for (size_t n=R*CX; n>0u; --n, ++X, ++Y) { *Y = *X; }
                for (size_t c=CX; c<C-C0; ++c)
                {
                    for (int n=(int)c-k; n>0; --n, ++X, ++Y) { *Y = 0.0; }
                    for (int n=(int)R-(int)c+k; n>0; --n, ++X, ++Y) { *Y = *X; }
                }
                for (size_t n=R*C0; n>0u; --n, ++X, ++Y) { *Y = 0.0; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0u;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0u;    //number of all-X rows
        for (size_t n=R0*C*SH; n>0u; --n, ++X, ++Y) { *Y = 0.0; }
        for (size_t r=R0; r<R-RX; ++r)
        {
            for (int n=(int)SH*((int)r+k+1); n>0; --n, ++X, ++Y) { *Y = *X; }
            for (int n=(int)SH*((int)C-(int)r-k-1); n>0; --n, ++X, ++Y) { *Y = 0.0; }
        }
        for (size_t n=RX*C*SH; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }

    return 0;
}


int tril_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_c: k must be in [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0u;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0u;                           //number of all-X cols
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s)
            {
                for (size_t n=R*CX; n>0u; --n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
                for (size_t c=CX; c<C-C0; ++c)
                {
                    for (int n=(int)c-k; n>0; --n, X+=2, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
                    for (int n=(int)R-(int)c+k; n>0; --n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
                }
                for (size_t n=R*C0; n>0u; --n, X+=2, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0u;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0u;    //number of all-X rows
        for (size_t n=R0*C*SH; n>0u; --n, X+=2, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
        for (size_t r=R0; r<R-RX; ++r)
        {
            for (int n=(int)SH*((int)r+k+1); n>0; --n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            for (int n=(int)SH*((int)C-(int)r-k-1); n>0; --n, X+=2, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
        }
        for (size_t n=RX*C*SH; n>0u; --n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
    }

    return 0;
}


int tril_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_z: k must be in [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0u;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0u;                           //number of all-X cols
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s)
            {
                for (size_t n=R*CX; n>0u; --n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
                for (size_t c=CX; c<C-C0; ++c)
                {
                    for (int n=(int)c-k; n>0; --n, X+=2, ++Y) { *Y = 0.0; *++Y = 0.0; }
                    for (int n=(int)R-(int)c+k; n>0; --n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
                }
                for (size_t n=R*C0; n>0u; --n, X+=2, ++Y) { *Y = 0.0; *++Y = 0.0; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0u;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0u;    //number of all-X rows
        for (size_t n=R0*C*SH; n>0u; --n, X+=2, ++Y) { *Y = 0.0; *++Y = 0.0; }
        for (size_t r=R0; r<R-RX; ++r)
        {
            for (int n=(int)SH*((int)r+k+1); n>0; --n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            for (int n=(int)SH*((int)C-(int)r-k-1); n>0; --n, X+=2, ++Y) { *Y = 0.0; *++Y = 0.0; }
        }
        for (size_t n=RX*C*SH; n>0u; --n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
    }

    return 0;
}


int tril_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_inplace_s: k must be in [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0u;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0u;                           //number of all-X cols
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s)
            {
                X += R*CX;
                for (size_t c=CX; c<C-C0; ++c)
                {
                    for (int n=(int)c-k; n>0; --n, ++X) { *X = 0.0f; }
                    X += (int)R - (int)c + k;
                }
                for (size_t n=R*C0; n>0u; --n, ++X) { *X = 0.0f; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0u;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0u;    //number of all-X rows
        for (size_t n=R0*C*SH; n>0u; --n, ++X) { *X = 0.0f; }
        for (size_t r=R0; r<R-RX; ++r)
        {
            X += (int)SH * ((int)r+k+1);
            for (int n=(int)SH*((int)C-(int)r-k-1); n>0; --n, ++X) { *X = 0.0f; }
        }
    }

    return 0;
}


int tril_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_inplace_d: k must be in [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0u;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0u;                           //number of all-X cols
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s)
            {
                X += R*CX;
                for (size_t c=CX; c<C-C0; ++c)
                {
                    for (int n=(int)c-k; n>0; --n, ++X) { *X = 0.0; }
                    X += (int)R - (int)c + k;
                }
                for (size_t n=R*C0; n>0u; --n, ++X) { *X = 0.0; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0u;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0u;    //number of all-X rows
        for (size_t n=R0*C*SH; n>0u; --n, ++X) { *X = 0.0; }
        for (size_t r=R0; r<R-RX; ++r)
        {
            X += (int)SH * ((int)r+k+1);
            for (int n=(int)SH*((int)C-(int)r-k-1); n>0; --n, ++X) { *X = 0.0; }
        }
    }

    return 0;
}


int tril_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_inplace_c: k must be in [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0u;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0u;                           //number of all-X cols
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s)
            {
                X += 2u*R*CX;
                for (size_t c=CX; c<C-C0; ++c)
                {
                    for (int n=(int)c-k; n>0; --n, ++X) { *X = 0.0f; *++X = 0.0f; }
                    X += 2*((int)R-(int)c+k);
                }
                for (size_t n=R*C0; n>0u; --n, ++X) { *X = 0.0f; *++X = 0.0f; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0u;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0u;    //number of all-X rows
        for (size_t n=R0*C*SH; n>0u; --n, ++X) { *X = 0.0f; *++X = 0.0f; }
        for (size_t r=R0; r<R-RX; ++r)
        {
            X += 2*(int)SH*((int)r+k+1);
            for (int n=(int)SH*((int)C-(int)r-k-1); n>0; --n, ++X) { *X = 0.0f; *++X = 0.0f; }
        }
    }

    return 0;
}


int tril_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_inplace_z: k must be in [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0u;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0u;                           //number of all-X cols
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s)
            {
                X += 2u*R*CX;
                for (size_t c=CX; c<C-C0; ++c)
                {
                    for (int n=(int)c-k; n>0; --n, ++X) { *X = 0.0; *++X = 0.0; }
                    X += 2*((int)R-(int)c+k);
                }
                for (size_t n=R*C0; n>0u; --n, ++X) { *X = 0.0; *++X = 0.0; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0u;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0u;    //number of all-X rows
        for (size_t n=R0*C*SH; n>0u; --n, ++X) { *X = 0.0; *++X = 0.0; }
        for (size_t r=R0; r<R-RX; ++r)
        {
            X += 2*(int)SH*((int)r+k+1);
            for (int n=(int)SH*((int)C-(int)r-k-1); n>0; --n, ++X) { *X = 0.0; *++X = 0.0; }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
