//Gets one row of input X

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int row_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_s: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        X += r;
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s)
            {
                for (size_t c=C; c>0u; --c, X+=R, ++Y) { *Y = *X; }
            }
        }
    }
    else
    {
        X += r*H*S*C;
        for (size_t c=H*S*C; c>0u; --c, ++X, ++Y) { *Y = *X; }
    }

    return 0;
}


int row_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_d: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        X += r;
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s)
            {
                for (size_t c=C; c>0u; --c, X+=R, ++Y) { *Y = *X; }
            }
        }
    }
    else
    {
        X += r*H*S*C;
        for (size_t c=H*S*C; c>0u; --c, ++X, ++Y) { *Y = *X; }
    }

    return 0;
}


int row_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_c: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        X += 2u*r;
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s)
            {
                for (size_t c=C; c>0u; --c, X+=2u*R-1u, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
    }
    else
    {
        X += 2u*r*H*S*C;
        for (size_t c=H*S*C; c>0u; --c, ++X, ++Y) { *Y = *X; *++Y = *++X; }
    }

    return 0;
}


int row_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_z: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        X += 2u*r;
        for (size_t h=H; h>0u; --h)
        {
            for (size_t s=S; s>0u; --s)
            {
                for (size_t c=C; c>0u; --c, X+=2u*R-1u, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
    }
    else
    {
        X += 2u*r*H*S*C;
        for (size_t c=H*S*C; c>0u; --c, ++X, ++Y) { *Y = *X; *++Y = *++X; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
