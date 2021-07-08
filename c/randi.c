//Gets random integers from a uniform distribution on [a,b).
//This uses modified code from PCG randoms minimal C library,
//but remains stand-alone (no install of PCG required).

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int randi_s (float *Y, const int a, const int b, const size_t N);
int randi_d (double *Y, const int a, const int b, const size_t N);
int randi_c (float *Y, const int a, const int b, const size_t N);
int randi_z (double *Y, const int a, const int b, const size_t N);


int randi_s (float *Y, const int a, const int b, const size_t N)
{
    if (a>b) { fprintf(stderr, "error in randi_s: a must be <= b\n"); return 1; }

    if (a==b)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = (float)a; }
    }
    else
    {
        const uint32_t bound = (uint32_t)(b-a+1);
        const uint32_t thresh = -bound % bound;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t mul = 6364136223846793005u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in randi_s: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        for (size_t n=0u; n<N; ++n, ++Y)
        {
            state = state*mul + inc;
            xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
            rot = state >> 59u;
            r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
            while (r<thresh)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
            }
            *Y = (float)((int)(r%bound) + a);
        }
    }

    return 0;
}


int randi_d (double *Y, const int a, const int b, const size_t N)
{
    if (a>b) { fprintf(stderr, "error in randi_d: a must be <= b\n"); return 1; }

    if (a==b)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = (double)a; }
    }
    else
    {
        const uint32_t bound = (uint32_t)(b-a+1);
        const uint32_t thresh = -bound % bound;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t mul = 6364136223846793005u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in randi_d: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        for (size_t n=0u; n<N; ++n, ++Y)
        {
            state = state*mul + inc;
            xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
            rot = state >> 59u;
            r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
            while (r<thresh)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
            }
            *Y = (double)((int)(r%bound) + a);
        }
    }
    
    return 0;
}


int randi_c (float *Y, const int a, const int b, const size_t N)
{
    if (a>b) { fprintf(stderr, "error in randi_c: a must be <= b\n"); return 1; }

    if (a==b)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = (float)a; }
    }
    else
    {
        const uint32_t bound = (uint32_t)(b-a+1);
        const uint32_t thresh = -bound % bound;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t mul = 6364136223846793005u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in randi_c: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        for (size_t n=0u; n<2u*N; ++n, ++Y)
        {
            state = state*mul + inc;
            xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
            rot = state >> 59u;
            r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
            while (r<thresh)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
            }
            *Y = (float)((int)(r%bound) + a);
        }
    }
    
    return 0;
}


int randi_z (double *Y, const int a, const int b, const size_t N)
{
    if (a>b) { fprintf(stderr, "error in randi_z: a must be <= b\n"); return 1; }

    if (a==b)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = (double)a; }
    }
    else
    {
        const uint32_t bound = (uint32_t)(b-a+1);
        const uint32_t thresh = -bound % bound;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t mul = 6364136223846793005u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in randi_z: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        for (size_t n=0u; n<2u*N; ++n, ++Y)
        {
            state = state*mul + inc;
            xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
            rot = state >> 59u;
            r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
            while (r<thresh)
            {
                state = state*mul + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
            }
            *Y = (double)((int)(r%bound) + a);
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
