//Gets one row of input X

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int row_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r);
int row_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r);
int row_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r);
int row_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r);


int row_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_s: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        X += r;
        for (size_t h=0u; h<H; ++h)
        {
            for (size_t s=0u; s<S; ++s)
            {
                for (size_t c=0u; c<C; ++c, X+=R, ++Y) { *Y = *X; }
            }
        }
    }
    else
    {
        X += r*H*S*C;
        for (size_t c=0u; c<H*S*C; ++c, ++X, ++Y) { *Y = *X; }
    }

    return 0;
}


int row_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_d: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        X += r;
        for (size_t h=0u; h<H; ++h)
        {
            for (size_t s=0u; s<S; ++s)
            {
                for (size_t c=0u; c<C; ++c, X+=R, ++Y) { *Y = *X; }
            }
        }
    }
    else
    {
        X += r*H*S*C;
        for (size_t c=0u; c<H*S*C; ++c, ++X, ++Y) { *Y = *X; }
    }

    return 0;
}


int row_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_c: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        X += 2*r;
        for (size_t h=0u; h<H; ++h)
        {
            for (size_t s=0u; s<S; ++s)
            {
                for (size_t c=0u; c<C; ++c, X+=2*R-1, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
    }
    else
    {
        X += 2*r*H*S*C;
        for (size_t c=0u; c<H*S*C; ++c, ++X, ++Y) { *Y = *X; *++Y = *++X; }
    }

    return 0;
}


int row_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_z: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        X += 2*r;
        for (size_t h=0u; h<H; ++h)
        {
            for (size_t s=0u; s<S; ++s)
            {
                for (size_t c=0u; c<C; ++c, X+=2*R-1, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
    }
    else
    {
        X += 2*r*H*S*C;
        for (size_t c=0u; c<H*S*C; ++c, ++X, ++Y) { *Y = *X; *++Y = *++X; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
