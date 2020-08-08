//Gets one column of input X

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int col_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c);
int col_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c);
int col_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c);
int col_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c);


int col_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c)
{
    if (c>C) { fprintf(stderr,"error in col_s: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (H*S==1)
    {
        if (iscolmajor)
        {
            X += c*R;
            for (size_t r=0; r<R; ++r, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            X += c;
            for (size_t r=0; r<R; ++r, X+=C, ++Y) { *Y = *X; }
        }
    }
    else if (iscolmajor)
    {
        X += c*R;
        for (size_t h=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s, X+=R*(C-1))
            {
                for (size_t r=0; r<R; ++r, ++X, ++Y) { *Y = *X; }
            }
        }
    }
    else
    {
        X += c*H*S;
        for (size_t r=0; r<R; ++r, X+=(C-1)*H*S)
        {
            for (size_t s=0; s<H*S; ++s, ++X, ++Y) { *Y = *X; }
        }
    }

    return 0;
}


int col_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c)
{
    if (c>C) { fprintf(stderr,"error in col_d: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (H*S==1)
    {
        if (iscolmajor)
        {
            X += c*R;
            for (size_t r=0; r<R; ++r, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            X += c;
            for (size_t r=0; r<R; ++r, X+=C, ++Y) { *Y = *X; }
        }
    }
    else if (iscolmajor)
    {
        X += c*R;
        for (size_t h=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s, X+=R*(C-1))
            {
                for (size_t r=0; r<R; ++r, ++X, ++Y) { *Y = *X; }
            }
        }
    }
    else
    {
        X += c*H*S;
        for (size_t r=0; r<R; ++r, X+=(C-1)*H*S)
        {
            for (size_t s=0; s<H*S; ++s, ++X, ++Y) { *Y = *X; }
        }
    }

    return 0;
}


int col_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c)
{
    if (c>C) { fprintf(stderr,"error in col_c: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (H*S==1)
    {
        if (iscolmajor)
        {
            X += 2*c*R;
            for (size_t r=0; r<R; ++r, ++X, ++Y) { *Y = *X; *++Y = *++X; }
        }
        else
        {
            X += 2*c;
            for (size_t r=0; r<R; ++r, X+=2*C-1, ++Y) { *Y = *X; *++Y = *++X; }
        }
    }
    else if (iscolmajor)
    {
        X += 2*c*R;
        for (size_t h=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s, X+=2*R*(C-1))
            {
                for (size_t r=0; r<R; ++r, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
    }
    else
    {
        X += 2*c*H*S;
        for (size_t r=0; r<R; ++r, X+=2*(C-1)*H*S)
        {
            for (size_t s=0; s<H*S; ++s, ++X, ++Y) { *Y = *X; *++Y = *++X; }
        }
    }

    return 0;
}


int col_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c)
{
    if (c>C) { fprintf(stderr,"error in col_z: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (H*S==1)
    {
        if (iscolmajor)
        {
            X += 2*c*R;
            for (size_t r=0; r<R; ++r, ++X, ++Y) { *Y = *X; *++Y = *++X; }
        }
        else
        {
            X += 2*c;
            for (size_t r=0; r<R; ++r, X+=2*C-1, ++Y) { *Y = *X; *++Y = *++X; }
        }
    }
    else if (iscolmajor)
    {
        X += 2*c*R;
        for (size_t h=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s, X+=2*R*(C-1))
            {
                for (size_t r=0; r<R; ++r, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
    }
    else
    {
        X += 2*c*H*S;
        for (size_t r=0; r<R; ++r, X+=2*(C-1)*H*S)
        {
            for (size_t s=0; s<H*S; ++s, ++X, ++Y) { *Y = *X; *++Y = *++X; }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
