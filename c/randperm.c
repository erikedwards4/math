//Generates integers from 1 to N in random order, and outputs the first M.
//This uses method of M Knuth shuffles.

//to use srand48 and drand48
#ifndef __cplusplus
    #ifndef _XOPEN_SOURCE
        #define _XOPEN_SOURCE
    #endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int randperm_s (float *Y, const size_t M, const size_t N);
int randperm_d (double *Y, const size_t M, const size_t N);


int randperm_s (float *Y, const size_t M, const size_t N)
{
    if (M>N) {  fprintf(stderr, "error in randperm_s: M must be <= N\n"); return 1; }

    //Generate ints 1:N
    size_t *X;
    if (!(X=(size_t *)malloc(N*sizeof(size_t)))) { fprintf(stderr,"error in randperm_s: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t n=0; n<N; n++) { X[n] = n + 1; }

    //Seed rand
    struct timespec ts;
	if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in randperm_s: timespec_get.\n"); perror("timespec_get"); return 1; }
	else { srand48(ts.tv_nsec^ts.tv_sec); }

    //M Knuth shuffles
    for (size_t m=0; m<M; m++)
	{
		size_t n = m + (size_t)((N-m)*drand48());  //random index from m to N-1
        Y[m] = (float)X[n]; X[n] = X[m];           //combine swap and output
	}

    free(X);
    return 0;
}


int randperm_d (double *Y, const size_t M, const size_t N)
{
    if (M>N) {  fprintf(stderr, "error in randperm_d: M must be <= N\n"); return 1; }

    //Generate ints 1:N
    size_t *X;
    if (!(X=(size_t *)malloc(N*sizeof(size_t)))) { fprintf(stderr,"error in randperm_d: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t n=0; n<N; n++) { X[n] = n + 1; }

    //Seed rand
    struct timespec ts;
	if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in randperm_d: timespec_get.\n"); perror("timespec_get"); return 1; }
	else { srand48(ts.tv_nsec^ts.tv_sec); }

    //M Knuth shuffles
    for (size_t m=0; m<M; m++)
	{
		size_t n = m + (size_t)((N-m)*drand48());  //random index from m to N-1
        Y[m] = (double)X[n]; X[n] = X[m];          //combine swap and output
	}

    free(X);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
