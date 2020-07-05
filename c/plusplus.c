//For each element of X, does x++ (increment by 1).
//This only has in-place version (since that is how used in C code).

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int plusplus_s (float *X, const size_t N);
int plusplus_d (double *X, const size_t N);


int plusplus_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n]++; }

    return 0;
}


int plusplus_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n]++; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
