//For each element of X, does x++ (increment by 1).
//This only has in-place version (since that is how used in C code).

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int plusplus_s (float *X, const int N);
int plusplus_d (double *X, const int N);


int plusplus_s (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in plusplus_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n]++; }

    return 0;
}


int plusplus_d (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in plusplus_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n]++; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
