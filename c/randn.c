//Gets random floats from a normal distribution with mu and sig.
//This uses modified code from PCG randoms minimal C library,
//but remains stand-alone (no install of PCG required).
//I use the textbook cos and sin formulae to make Gaussian.

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int randn_s (float *Y, const float mu, const float sig, const size_t N);
int randn_d (double *Y, const double mu, const double sig, const size_t N);
int randn_c (float *Y, const float mu, const float sig, const size_t N);
int randn_z (double *Y, const double mu, const double sig, const size_t N);


int randn_s (float *Y, const float mu, const float sig, const size_t N)
{
    if (sig<0.0f) { fprintf(stderr, "error in randn_s: sig must be nonnegative\n"); return 1; }

    if (sig<=0.0f)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = mu; }
    }
    else
    {
        const float M_2PI = (float)(2.0*M_PI);
        float u1, u2, R;
        size_t n;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in randn_s: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        if (mu==0.0f && sig==1.0f)
        {
            for (n=0u; n<N-1u; n+=2u)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = sqrtf(-2.0f*logf(1.0f-u1));
                *Y++ = R * cosf(M_2PI*u2);
                *Y++ = R * sinf(M_2PI*u2);
            }
            if (n<N)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = sqrtf(-2.0f*logf(1.0f-u1));
                *Y++ = R * cosf(M_2PI*u2);
            }
        }
        else
        {
            for (n=0u; n<N-1u; n+=2u)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = sig * sqrtf(-2.0f*logf(1.0f-u1));
                *Y++ = R*cosf(M_2PI*u2) + mu;
                *Y++ = R*sinf(M_2PI*u2) + mu;
            }
            if (n<N)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = sig * sqrtf(-2.0f*logf(1.0f-u1));
                *Y++ = R*cosf(M_2PI*u2) + mu;
            }
        }
    }

    return 0;
}


int randn_d (double *Y, const double mu, const double sig, const size_t N)
{
    if (sig<0.0) { fprintf(stderr, "error in randn_d: sig must be nonnegative\n"); return 1; }

    if (sig<=0.0)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = mu; }
    }
    else
    {
        const double M_2PI = 2.0*M_PI;
        size_t n;
        double u1, u2, R;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in randn_d: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        if (mu==0.0 && sig==1.0)
        {
            for (n=0u; n<N-1u; n+=2u)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = sqrt(-2.0*log(1.0-u1));
                *Y++ = R * cos(M_2PI*u2);
                *Y++ = R * sin(M_2PI*u2);
            }
            if (n<N)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = sqrt(-2.0*log(1.0-u1));
                *Y++ = R * cos(M_2PI*u2);
            }
        }
        else
        {
            for (n=0u; n<N-1u; n+=2u)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = sig * sqrt(-2.0*log(1.0-u1));
                *Y++ = R*cos(M_2PI*u2) + mu;
                *Y++ = R*sin(M_2PI*u2) + mu;
            }
            if (n<N)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = sig * sqrt(-2.0*log(1.0-u1));
                *Y++ = R*cos(M_2PI*u2) + mu;
            }
        }
    }

    return 0;
}


int randn_c (float *Y, const float mu, const float sig, const size_t N)
{
    if (sig<0.0f) { fprintf(stderr, "error in randn_s: sig must be nonnegative\n"); return 1; }

    if (sig<=0.0f)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = mu; }
    }
    else
    {
        const float M_2PI = (float)(2.0*M_PI);
        float u1, u2, R;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in randn_c: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        if (mu==0.0f && sig==1.0f)
        {
            for (size_t n=0u; n<N; ++n)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = sqrtf(-2.0f*logf(1.0f-u1));
                *Y++ = R * cosf(M_2PI*u2);
                *Y++ = R * sinf(M_2PI*u2);
            }
        }
        else
        {
            for (size_t n=0u; n<N; ++n)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((float)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((float)r,-32);
                R = sig * sqrtf(-2.0f*logf(1.0f-u1));
                *Y++ = R*cosf(M_2PI*u2) + mu;
                *Y++ = R*sinf(M_2PI*u2) + mu;
            }
        }
    }

    return 0;
}


int randn_z (double *Y, const double mu, const double sig, const size_t N)
{
    if (sig<0.0) { fprintf(stderr, "error in randn_z: sig must be nonnegative\n"); return 1; }

    if (sig<=0.0)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = mu; }
    }
    else
    {
        const double M_2PI = 2.0*M_PI;
        double u1, u2, R;
        uint32_t r, xorshifted, rot;
        uint64_t state = 0u;
        const uint64_t inc = ((uint64_t)(&state) << 1u) | 1u;
        struct timespec ts;

        //Init random num generator
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr, "error in randn_z: timespec_get.\n"); perror("timespec_get"); return 1; }
	    state = (uint64_t)(ts.tv_nsec^ts.tv_sec) + inc;

        //Generate
        if (mu==0.0 && sig==1.0)
        {
            for (size_t n=0u; n<N; ++n)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = sqrt(-2.0*log(1.0-u1));
                *Y++ = R * cos(M_2PI*u2);
                *Y++ = R * sin(M_2PI*u2);
            }
        }
        else
        {
            for (size_t n=0u; n<N; ++n)
            {
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u1 = ldexp((double)r,-32);
                state = state*6364136223846793005ull + inc;
                xorshifted = (uint32_t)(((state >> 18u) ^ state) >> 27u);
                rot = state >> 59u;
                r = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
                u2 = ldexp((double)r,-32);
                R = sig * sqrt(-2.0*log(1.0-u1));
                *Y++ = R*cos(M_2PI*u2) + mu;
                *Y++ = R*sin(M_2PI*u2) + mu;
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
