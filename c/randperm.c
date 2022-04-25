//Generates integers from 1 to N in random order, and outputs the first M.
//This uses method of M Knuth shuffles.
//For each shuffle, a random index from m to N-1 is needed,
//so I use the PCG-random-based code of randi.c.

//to use srand48 and drand48
// #ifndef __cplusplus
//     #ifndef _XOPEN_SOURCE
//         #define _XOPEN_SOURCE
//     #endif
// #endif

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int randperm_s (float *Y, const size_t M, const size_t N)
{
    if (M>N) { fprintf(stderr, "error in randperm_s: M must be <= N\n"); return 1; }

    uint32_t r, xorshifted, rot, bound, thresh;
    uint64_t state = 0u;
    const uint64_t mul = 6364136223846793005u;
    const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
    struct timespec ts;
    size_t k;

    //Generate ints 1:N
    size_t *X;
    if (!(X=(size_t *)malloc(N*sizeof(size_t)))) { fprintf(stderr,"error in randperm_s: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t n=1u; n<=N; ++n, ++X) { *X = n; }
    X -= N;

    //Init random num generator
    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in randperm_s: timespec_get.\n"); perror("timespec_get"); return 1; }
    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

    //M Knuth shuffles
    for (size_t m=0u; m<M; ++m, ++X, ++Y)
	{
        //Get random index k from m to N-1 (m<=k<N)
		//k = (size_t)((N-m)*drand48());
        state = state*mul + inc;
        xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
        rot = state >> 59u;
        r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
        bound = (uint32_t)(N-m);
        thresh = -bound % bound;
        while (r<thresh)
        {
            state = state*mul + inc;
            xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
            rot = state >> 59u;
            r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
        }
        k = (size_t)(r%bound);

        //Shuffle
        *Y = (float)*(X+k);
        *(X+k) = *X;
	}

    //Free
    X -= M;
    free(X);

    return 0;
}


int randperm_d (double *Y, const size_t M, const size_t N)
{
    if (M>N) { fprintf(stderr, "error in randperm_d: M must be <= N\n"); return 1; }

    uint32_t r, xorshifted, rot, bound, thresh;
    uint64_t state = 0u;
    const uint64_t mul = 6364136223846793005u;
    const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
    struct timespec ts;
    size_t k;

    //Generate ints 1:N
    size_t *X;
    if (!(X=(size_t *)malloc(N*sizeof(size_t)))) { fprintf(stderr,"error in randperm_d: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t n=1u; n<=N; ++n, ++X) { *X = n; }
    X -= N;

    //Init random num generator
    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in randperm_d: timespec_get.\n"); perror("timespec_get"); return 1; }
    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;
    //srand48(ts.tv_nsec^ts.tv_sec);

    //M Knuth shuffles
    for (size_t m=0u; m<M; ++m, ++X, ++Y)
	{
		//Get random index k from m to N-1 (m<=k<N)
		//k = (size_t)((N-m)*drand48());
        state = state*mul + inc;
        xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
        rot = state >> 59u;
        r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
        bound = (uint32_t)(N-m);
        thresh = -bound % bound;
        while (r<thresh)
        {
            state = state*mul + inc;
            xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
            rot = state >> 59u;
            r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
        }
        k = (size_t)(r%bound);

        //Shuffle
        *Y = (double)*(X+k);
        *(X+k) = *X;
	}

    //Free
    X -= M;
    free(X);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
