//For each element of X, does x-- (decrement by 1).
//This only has in-place version (since that is how used in C code).

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int minusminus_s (float *X, const int N);
int minusminus_d (double *X, const int N);


int minusminus_s (float *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in minusminus_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { X[n]--; }

    return 0;
}


int minusminus_d (double *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in minusminus_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { X[n]--; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
