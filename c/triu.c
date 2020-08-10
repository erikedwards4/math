//Zeros all elements below the kth diagonal of matrix X.
//This has in-place and not-in-place versions.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int triu_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int triu_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int triu_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int triu_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);

int triu_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int triu_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int triu_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int triu_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);


int triu_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_s: k must be in [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;                                    //number of all-0 cols
        const size_t CX = (k>(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1);    //number of all-X cols
        for (size_t h=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                for (size_t n=0; n<R*C0; ++n, ++X, ++Y) { *Y = 0.0f; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    for (int n=0; n<(int)c-k+1; ++n, ++X, ++Y) { *Y = *X; }
                    for (int n=0; n<(int)R-(int)c+k-1; ++n, ++X, ++Y) { *Y = 0.0f; }
                }
                for (size_t n=0; n<R*CX; ++n, ++X, ++Y) { *Y = *X; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)R-(int)C) ? (size_t)((int)R-(int)C+k) : 0;    //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(1-k);                            //number of all-X rows
        for (size_t n=0; n<RX*C*SH; ++n, ++X, ++Y) { *Y = *X; }
        for (size_t r=RX; r<R-R0; ++r)
        {
            for (int n=0; n<(int)SH*((int)r+k); ++n, ++X, ++Y) { *Y = 0.0f; }
            for (int n=0; n<(int)SH*((int)C-(int)r-k); ++n, ++X, ++Y) { *Y = *X; }
        }
        for (size_t n=0; n<R0*C*SH; ++n, ++Y) { *Y = 0.0f; }
    }

    return 0;
}


int triu_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_d: k must be in [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;                                    //number of all-0 cols
        const size_t CX = (k>(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1);    //number of all-X cols
        for (size_t h=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                for (size_t n=0; n<R*C0; ++n, ++X, ++Y) { *Y = 0.0; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    for (int n=0; n<(int)c-k+1; ++n, ++X, ++Y) { *Y = *X; }
                    for (int n=0; n<(int)R-(int)c+k-1; ++n, ++X, ++Y) { *Y = 0.0; }
                }
                for (size_t n=0; n<R*CX; ++n, ++X, ++Y) { *Y = *X; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)R-(int)C) ? (size_t)((int)R-(int)C+k) : 0;    //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(1-k);                            //number of all-X rows
        for (size_t n=0; n<RX*C*SH; ++n, ++X, ++Y) { *Y = *X; }
        for (size_t r=RX; r<R-R0; ++r)
        {
            for (int n=0; n<(int)SH*((int)r+k); ++n, ++X, ++Y) { *Y = 0.0; }
            for (int n=0; n<(int)SH*((int)C-(int)r-k); ++n, ++X, ++Y) { *Y = *X; }
        }
        for (size_t n=0; n<R0*C*SH; ++n, ++Y) { *Y = 0.0; }
    }

    return 0;
}


int triu_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_c: k must be in [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;                                    //number of all-0 cols
        const size_t CX = (k>(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1);    //number of all-X cols
        for (size_t h=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                for (size_t n=0; n<R*C0; ++n, X+=2, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    for (int n=0; n<(int)c-k+1; ++n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
                    for (int n=0; n<(int)R-(int)c+k-1; ++n, X+=2, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
                }
                for (size_t n=0; n<R*CX; ++n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)R-(int)C) ? (size_t)((int)R-(int)C+k) : 0;    //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(1-k);                            //number of all-X rows
        for (size_t n=0; n<RX*C*SH; ++n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
        for (size_t r=RX; r<R-R0; ++r)
        {
            for (int n=0; n<(int)SH*((int)r+k); ++n, X+=2, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
            for (int n=0; n<(int)SH*((int)C-(int)r-k); ++n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
        }
        for (size_t n=0; n<R0*C*SH; ++n, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
    }

    return 0;
}


int triu_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_z: k must be in [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;                                    //number of all-0 cols
        const size_t CX = (k>(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1);    //number of all-X cols
        for (size_t h=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                for (size_t n=0; n<R*C0; ++n, X+=2, ++Y) { *Y = 0.0; *++Y = 0.0; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    for (int n=0; n<(int)c-k+1; ++n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
                    for (int n=0; n<(int)R-(int)c+k-1; ++n, X+=2, ++Y) { *Y = 0.0; *++Y = 0.0; }
                }
                for (size_t n=0; n<R*CX; ++n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)R-(int)C) ? (size_t)((int)R-(int)C+k) : 0;    //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(1-k);                            //number of all-X rows
        for (size_t n=0; n<RX*C*SH; ++n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
        for (size_t r=RX; r<R-R0; ++r)
        {
            for (int n=0; n<(int)SH*((int)r+k); ++n, X+=2, ++Y) { *Y = 0.0; *++Y = 0.0; }
            for (int n=0; n<(int)SH*((int)C-(int)r-k); ++n, ++X, ++Y) { *Y = *X; *++Y = *++X; }
        }
        for (size_t n=0; n<R0*C*SH; ++n, ++Y) { *Y = 0.0; *++Y = 0.0; }
    }

    return 0;
}


int triu_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_inplace_s: k must be in [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;                                    //number of all-0 cols
        const size_t CX = (k>(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1);    //number of all-X cols
        for (size_t h=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                for (size_t n=0; n<R*C0; ++n, ++X) { *X = 0.0f; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    X += (int)c - k + 1;
                    for (int n=0; n<(int)R-(int)c+k-1; ++n, ++X) { *X = 0.0f; }
                }
                X += R*CX;
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)R-(int)C) ? (size_t)((int)R-(int)C+k) : 0;    //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(1-k);                            //number of all-X rows
        X += RX*C*SH;
        for (size_t r=RX; r<R-R0; ++r)
        {
            for (int n=0; n<(int)SH*((int)r+k); ++n, ++X) { *X = 0.0f; }
            X += (int)SH*((int)C-(int)r-k);
        }
        for (size_t n=0; n<R0*C*SH; ++n, ++X) { *X = 0.0f; }
    }

    return 0;
}


int triu_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_inplace_d: k must be in [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;                                    //number of all-0 cols
        const size_t CX = (k>(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1);    //number of all-X cols
        for (size_t h=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                for (size_t n=0; n<R*C0; ++n, ++X) { *X = 0.0; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    X += (int)c - k + 1;
                    for (int n=0; n<(int)R-(int)c+k-1; ++n, ++X) { *X = 0.0; }
                }
                X += R*CX;
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)R-(int)C) ? (size_t)((int)R-(int)C+k) : 0;    //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(1-k);                            //number of all-X rows
        X += RX*C*SH;
        for (size_t r=RX; r<R-R0; ++r)
        {
            for (int n=0; n<(int)SH*((int)r+k); ++n, ++X) { *X = 0.0; }
            X += (int)SH*((int)C-(int)r-k);
        }
        for (size_t n=0; n<R0*C*SH; ++n, ++X) { *X = 0.0; }
    }

    return 0;
}


int triu_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_inplace_c: k must be in [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;                                    //number of all-0 cols
        const size_t CX = (k>(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1);    //number of all-X cols
        for (size_t h=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                for (size_t n=0; n<R*C0; ++n, ++X) { *X = 0.0f; *++X = 0.0f; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    X += 2*((int)c-k+1);
                    for (int n=0; n<(int)R-(int)c+k-1; ++n, ++X) { *X = 0.0f; *++X = 0.0f; }
                }
                X += 2*R*CX;
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)R-(int)C) ? (size_t)((int)R-(int)C+k) : 0;    //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(1-k);                            //number of all-X rows
        X += 2*RX*C*SH;
        for (size_t r=RX; r<R-R0; ++r)
        {
            for (int n=0; n<(int)SH*((int)r+k); ++n, ++X) { *X = 0.0f; *++X = 0.0f; }
            X += 2*(int)SH*((int)C-(int)r-k);
        }
        for (size_t n=0; n<R0*C*SH; ++n, ++X) { *X = 0.0f; *++X = 0.0f; }
    }

    return 0;
}


int triu_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_inplace_z: k must be in [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;                                    //number of all-0 cols
        const size_t CX = (k>(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1);    //number of all-X cols
        for (size_t h=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                for (size_t n=0; n<R*C0; ++n, ++X) { *X = 0.0; *++X = 0.0; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    X += 2*((int)c-k+1);
                    for (int n=0; n<(int)R-(int)c+k-1; ++n, ++X) { *X = 0.0; *++X = 0.0; }
                }
                X += 2*R*CX;
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)R-(int)C) ? (size_t)((int)R-(int)C+k) : 0;    //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(1-k);                            //number of all-X rows
        X += 2*RX*C*SH;
        for (size_t r=RX; r<R-R0; ++r)
        {
            for (int n=0; n<(int)SH*((int)r+k); ++n, ++X) { *X = 0.0; *++X = 0.0; }
            X += 2*(int)SH*((int)C-(int)r-k);
        }
        for (size_t n=0; n<R0*C*SH; ++n, ++X) { *X = 0.0; *++X = 0.0; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
