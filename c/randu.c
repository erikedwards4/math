//Gets random floats from a uniform distribution on [a,b).
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

int randu_s (float *Y, const float a, const float b, const size_t N);
int randu_d (double *Y, const double a, const double b, const size_t N);
int randu_c (float *Y, const float a, const float b, const size_t N);
int randu_z (double *Y, const double a, const double b, const size_t N);


int randu_s (float *Y, const float a, const float b, const size_t N)
{
    if (a>b) { fprintf(stderr, "error in randu_s: a must be <= b\n"); return 1; }

    uint64_t t0, oldstate=0U, state=0U, inc;
    uint32_t xorshifted, rot, yi, threshold;

    threshold = -bound % bound;

    //Init random num generator
    #ifdef __cplusplus
    t0 = (uint64_t)time(nullptr);
    #else
    t0 = (uint64_t)time(NULL);
    #endif
    inc = ((uint64_t)(&inc) << 1u) | 1u;
    state = oldstate*6364136223846793005ULL + inc;
    state += t0;
    oldstate = state;
    state = oldstate*6364136223846793005ULL + inc;

    for (size_t n=0u; n<N; ++n, ++Y)
    {
        oldstate = state;
        state = oldstate*6364136223846793005ULL + inc;
        xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
        rot = oldstate >> 59u;
        yi = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
        *Y = (float)ldexp(yi,-32);
    }

    return 0;
}


int randu_d (double *Y, const double a, const double b, const size_t N)
{
    if (a>b) { fprintf(stderr, "error in randu_d: a must be <= b\n"); return 1; }

    uint64_t t0, oldstate=0U, state=0U, inc;
    uint32_t xorshifted, rot, yi;

    //Init random num generator
    #ifdef __cplusplus
    t0 = (uint64_t)time(nullptr);
    #else
    t0 = (uint64_t)time(NULL);
    #endif
    inc = ((uint64_t)(&inc) << 1u) | 1u;
    state = oldstate*6364136223846793005ULL + inc;
    state += t0;
    oldstate = state;
    state = oldstate*6364136223846793005ULL + inc;

    for (size_t n=0u; n<N; ++n, ++Y)
    {
        oldstate = state;
        state = oldstate*6364136223846793005ULL + inc;
        xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
        rot = oldstate >> 59u;
        yi = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
        *Y = ldexp(yi,-32);
    }
    
    return 0;
}


int randu_c (float *Y, const float a, const float b, const size_t N)
{
    if (a>b) { fprintf(stderr, "error in randu_c: a must be <= b\n"); return 1; }

    uint64_t t0, oldstate=0U, state=0U, inc;
    uint32_t xorshifted, rot, yi;

    //Init random num generator
    #ifdef __cplusplus
    t0 = (uint64_t)time(nullptr);
    #else
    t0 = (uint64_t)time(NULL);
    #endif
    inc = ((uint64_t)(&inc) << 1u) | 1u;
    state = oldstate*6364136223846793005ULL + inc;
    state += t0;
    oldstate = state;
    state = oldstate*6364136223846793005ULL + inc;

    for (size_t n=0u; n<2u*N; ++n, ++Y)
    {
        oldstate = state;
        state = oldstate*6364136223846793005ULL + inc;
        xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
        rot = oldstate >> 59u;
        yi = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
        *Y = (float)ldexp(yi,-32);
    }
    
    return 0;
}


int randu_z (double *Y, const double a, const double b, const size_t N)
{
    if (a>b) { fprintf(stderr, "error in randu_z: a must be <= b\n"); return 1; }

    uint64_t t0, oldstate=0U, state=0U, inc;
    uint32_t xorshifted, rot, yi;

    //Init random num generator
    #ifdef __cplusplus
    t0 = (uint64_t)time(nullptr);
    #else
    t0 = (uint64_t)time(NULL);
    #endif
    inc = ((uint64_t)(&inc) << 1u) | 1u;
    state = oldstate*6364136223846793005ULL + inc;
    state += t0;
    oldstate = state;
    state = oldstate*6364136223846793005ULL + inc;

    for (size_t n=0u; n<2u*N; ++n, ++Y)
    {
        oldstate = state;
        state = oldstate*6364136223846793005ULL + inc;
        xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
        rot = oldstate >> 59u;
        yi = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
        *Y = ldexp(yi,-32);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
