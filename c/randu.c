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

    if (a>=b)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = a; }
    }
    else
    {
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;

        //Init random num generator
        #ifdef __cplusplus
            state = (uint64_t)time(nullptr) + inc;
        #else
            state = (uint64_t)time(NULL) + inc;
        #endif
        state = state*6364136223846793005ull + inc;

        //Generate
        if (a==0.0f && b==1.0f)
        {
            for (size_t n=0u; n<N; ++n, ++Y)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                *Y = ldexp((float)r,-32);
            }
        }
        else
        {
            const float sc = b - a;
            for (size_t n=0u; n<N; ++n, ++Y)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                *Y = ldexp((float)r,-32)*sc + a;
            }
        }
    }

    return 0;
}


int randu_d (double *Y, const double a, const double b, const size_t N)
{
    if (a>b) { fprintf(stderr, "error in randu_d: a must be <= b\n"); return 1; }

    if (a>=b)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = a; }
    }
    else
    {
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;

        //Init random num generator
        #ifdef __cplusplus
            state = (uint64_t)time(nullptr) + inc;
        #else
            state = (uint64_t)time(NULL) + inc;
        #endif
        state = state*6364136223846793005ull + inc;

        //Generate
        if (a==0.0 && b==1.0)
        {
            for (size_t n=0u; n<N; ++n, ++Y)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                *Y = ldexp((double)r,-32);
            }
        }
        else
        {
            const double sc = b - a;
            for (size_t n=0u; n<N; ++n, ++Y)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                *Y = ldexp((double)r,-32)*sc + a;
            }
        }
    }

    return 0;
}


int randu_c (float *Y, const float a, const float b, const size_t N)
{
    if (a>b) { fprintf(stderr, "error in randu_c: a must be <= b\n"); return 1; }

    if (a>=b)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = a; }
    }
    else
    {
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;

        //Init random num generator
        #ifdef __cplusplus
            state = (uint64_t)time(nullptr) + inc;
        #else
            state = (uint64_t)time(NULL) + inc;
        #endif
        state = state*6364136223846793005ull + inc;

        //Generate
        if (a==0.0f && b==1.0f)
        {
            for (size_t n=0u; n<2u*N; ++n, ++Y)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                *Y = ldexp((float)r,-32);
            }
        }
        else
        {
            const float sc = b - a;
            for (size_t n=0u; n<N; ++n, ++Y)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                *Y = ldexp((float)r,-32)*sc + a;
            }
        }
    }
    
    return 0;
}


int randu_z (double *Y, const double a, const double b, const size_t N)
{
    if (a>b) { fprintf(stderr, "error in randu_z: a must be <= b\n"); return 1; }

    if (a>=b)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = a; }
    }
    else
    {
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;

        //Init random num generator
        #ifdef __cplusplus
            state = (uint64_t)time(nullptr) + inc;
        #else
            state = (uint64_t)time(NULL) + inc;
        #endif
        state = state*6364136223846793005ull + inc;

        //Generate
        if (a==0.0 && b==1.0)
        {
            for (size_t n=0u; n<2u*N; ++n, ++Y)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                *Y = ldexp((double)r,-32);
            }
        }
        else
        {
            const double sc = b - a;
            for (size_t n=0u; n<N; ++n, ++Y)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                *Y = ldexp((double)r,-32)*sc + a;
            }
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
