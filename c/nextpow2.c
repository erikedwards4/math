//This gets next power-of-2 for each element of X.
//This is for int or uint, so I only make .c function to show.
//Note that Octave returns the power itself, i.e. log2 of what this returns.
//However, I return 0 for X=0, whereas Octave returns log2(1)=0 for X=0. 
//Also, I deal with negative numbers here (-3 -> -4),
//whereas Octave is positive only (-3 -> log2(4)).

#include <stdio.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int nextpow2_i (int *Y, const int *X, const size_t N);
int nextpow2_u (size_t *Y, const size_t *X, const size_t N);


int nextpow2_i (int *Y, const int *X, const size_t N)
{
    int x, y;

    for (size_t n=N; n>0u; --n, ++X, ++Y)
    {
        if (*X==0) { *Y = 0u; }
        else if (*X<0)
        {
            x = -*X; y = 1;
            while (y<x) { y *= 2; }
            *Y = -y;
        }
        else
        {
            x = *X; y = 1;
            while (y<x) { y *= 2; }
            *Y = y;
        }
    }

    return 0;
}


int nextpow2_u (size_t *Y, const size_t *X, const size_t N)
{
    size_t x, y;
    
    for (size_t n=N; n>0u; --n, ++X, ++Y)
    {
        if (*X==0u) { *Y = 0u; }
        else
        {
            x = *X; y = 1u;
            while (y<x) { y *= 2u; }
            *Y = y;
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
