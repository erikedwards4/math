//Gets one column of input X

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int col_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t c);
int col_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t c);
int col_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t c);
int col_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t c);


int col_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t c)
{
    if (c>C) { fprintf(stderr,"error in col_s: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (H*S==1u)
    {
        if (iscolmajor)
        {
            X += c*R;
            for (size_t r=R; r>0u; --r, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            X += c;
            for (size_t r=R; r>0u; --r, X+=C, ++Y) { *Y = *X; }
        }
    }
    else if (iscolmajor)
    {
        X += c*R;
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s, X+=R*(C-1u))
            {
                for (size_t r=R; r>0u; --r, ++X, ++Y) { *Y = *X; }
            }
        }
    }
    else
    {
        X += c*H*S;
        for (size_t r=R; r>0u; --r, X+=(C-1u)*H*S)
        {
            for (size_t s=H*S; s>0u; --s, ++X, ++Y) { *Y = *X; }
        }
    }

    return 0;
}


int col_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t c)
{
    if (c>C) { fprintf(stderr,"error in col_d: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (H*S==1u)
    {
        if (iscolmajor)
        {
            X += c*R;
            for (size_t r=R; r>0u; --r, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            X += c;
            for (size_t r=R; r>0u; --r, X+=C, ++Y) { *Y = *X; }
        }
    }
    else if (iscolmajor)
    {
        X += c*R;
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s, X+=R*(C-1u))
            {
                for (size_t r=R; r>0u; --r, ++X, ++Y) { *Y = *X; }
            }
        }
    }
    else
    {
        X += c*H*S;
        for (size_t r=R; r>0u; --r, X+=(C-1u)*H*S)
        {
            for (size_t s=H*S; s>0u; --s, ++X, ++Y) { *Y = *X; }
        }
    }

    return 0;
}


int col_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t c)
{
    if (c>C) { fprintf(stderr,"error in col_c: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (H*S==1u)
    {
        if (iscolmajor)
        {
            X += 2u*c*R;
            for (size_t r=R; r>0u; --r, ++X, ++Y) { *Y = *X; *++Y = *++X; }
        }
        else
        {
            X += 2u*c;
            for (size_t r=R; r>0u; --r, X+=2u*C-1u, ++Y) { *Y = *X; *++Y = *++X; }
        }
    }
    else if (iscolmajor)
    {
        X += 2u*c*R;
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s, X+=2u*R*(C-1u))
            {
                for (size_t r=R; r>0u; --r, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
    }
    else
    {
        X += 2u*c*H*S;
        for (size_t r=R; r>0u; --r, X+=2u*(C-1u)*H*S)
        {
            for (size_t s=H*S; s>0u; --s, ++X, ++Y) { *Y = *X; *++Y = *++X; }
        }
    }

    return 0;
}


int col_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t c)
{
    if (c>C) { fprintf(stderr,"error in col_z: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (H*S==1u)
    {
        if (iscolmajor)
        {
            X += 2u*c*R;
            for (size_t r=R; r>0u; --r, ++X, ++Y) { *Y = *X; *++Y = *++X; }
        }
        else
        {
            X += 2u*c;
            for (size_t r=R; r>0u; --r, X+=2u*C-1u, ++Y) { *Y = *X; *++Y = *++X; }
        }
    }
    else if (iscolmajor)
    {
        X += 2u*c*R;
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s, X+=2u*R*(C-1u))
            {
                for (size_t r=R; r>0u; --r, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
    }
    else
    {
        X += 2u*c*H*S;
        for (size_t r=R; r>0u; --r, X+=2u*(C-1u)*H*S)
        {
            for (size_t s=H*S; s>0u; --s, ++X, ++Y) { *Y = *X; *++Y = *++X; }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
