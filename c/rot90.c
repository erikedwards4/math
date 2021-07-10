//Rotates matrix X by 90 degrees K times.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int rot90_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int K);
int rot90_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int K);
int rot90_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int K);
int rot90_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int K);


int rot90_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int K)
{
    const size_t N = R*C;
    const int k = (K<0) ? 4+K%4 : K%4;

    if (N==0u) {}
    else if (k==0 || N==1)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (k==1)
    {
        if (iscolmajor)
        {
            Y += N - C;
            for (size_t r=0u; r<R; ++r, X-=N-1, Y-=2*C)
            {
                for (size_t c=0u; c<C; ++c, X+=R, ++Y) { *Y = *X; }
            }
        }
        else
        {
            Y += R - 1;
            for (size_t r=0u; r<R; ++r, Y-=N+1)
            {
                for (size_t c=0u; c<C; ++c, ++X, Y+=R) { *Y = *X; }
            }
        }
    }
    else if (k==2)
    {
        Y += N - 1;
        for (size_t n=0u; n<N; ++n, ++X, --Y) { *Y = *X; }
    }
    else // (k==3)
    {
        if (iscolmajor)
        {
            Y += C - 1;
            for (size_t c=0u; c<C; ++c, Y-=N+1)
            {
                for (size_t r=0u; r<R; ++r, ++X, Y+=C) { *Y = *X; }
            }
        }
        else
        {
            Y += N - R;
            for (size_t c=0u; c<C; ++c, X-=N-1, Y-=2*R)
            {
                for (size_t r=0u; r<R; ++r, X+=C, ++Y) { *Y = *X; }
            }
        }
    }

    return 0;
}


int rot90_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int K)
{
    const size_t N = R*C;
    const int k = (K<0) ? 4+K%4 : K%4;

    if (N==0u) {}
    else if (k==0 || N==1)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (k==1)
    {
        if (iscolmajor)
        {
            Y += N - C;
            for (size_t r=0u; r<R; ++r, X-=N-1, Y-=2*C)
            {
                for (size_t c=0u; c<C; ++c, X+=R, ++Y) { *Y = *X; }
            }
        }
        else
        {
            Y += R - 1;
            for (size_t r=0u; r<R; ++r, Y-=N+1)
            {
                for (size_t c=0u; c<C; ++c, ++X, Y+=R) { *Y = *X; }
            }
        }
    }
    else if (k==2)
    {
        Y += N - 1;
        for (size_t n=0u; n<N; ++n, ++X, --Y) { *Y = *X; }
    }
    else // (k==3)
    {
        if (iscolmajor)
        {
            Y += C - 1;
            for (size_t c=0u; c<C; ++c, Y-=N+1)
            {
                for (size_t r=0u; r<R; ++r, ++X, Y+=C) { *Y = *X; }
            }
        }
        else
        {
            Y += N - R;
            for (size_t c=0u; c<C; ++c, X-=N-1, Y-=2*R)
            {
                for (size_t r=0u; r<R; ++r, X+=C, ++Y) { *Y = *X; }
            }
        }
    }

    return 0;
}


int rot90_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int K)
{
    const size_t N = R*C;
    const int k = (K<0) ? 4+K%4 : K%4;

    if (N==0u) {}
    else if (k==0 || N==1)
    {
        for (size_t n=0u; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (k==1)
    {
        if (iscolmajor)
        {
            Y += 2*(N-C);
            for (size_t r=0u; r<R; ++r, X-=2*N-2, Y-=4*C)
            {
                for (size_t c=0u; c<C; ++c, X+=2*R-1, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
        else
        {
            Y += 2*(R-1);
            for (size_t r=0u; r<R; ++r, Y-=2*N+2)
            {
                for (size_t c=0u; c<C; ++c, ++X, Y+=2*R-1) { *Y = *X; *++Y = *++X; }
            }
        }
    }
    else if (k==2)
    {
        Y += 2*(N-1);
        for (size_t n=0u; n<N; ++n, ++X, Y-=2) { *Y = *X; *(Y+1) = *++X; }
    }
    else // (k==3)
    {
        if (iscolmajor)
        {
            Y += 2*(C-1);
            for (size_t c=0u; c<C; ++c, Y-=2*N+2)
            {
                for (size_t r=0u; r<R; ++r, ++X, Y+=2*C-1) { *Y = *X; *++Y = *++X; }
            }
        }
        else
        {
            Y += 2*(N-R);
            for (size_t c=0u; c<C; ++c, X-=2*N-2, Y-=4*R)
            {
                for (size_t r=0u; r<R; ++r, X+=2*C-1, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
    }

    return 0;
}


int rot90_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int K)
{
    const size_t N = R*C;
    const int k = (K<0) ? 4+K%4 : K%4;

    if (N==0u) {}
    else if (k==0 || N==1)
    {
        for (size_t n=0u; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (k==1)
    {
        if (iscolmajor)
        {
            Y += 2*(N-C);
            for (size_t r=0u; r<R; ++r, X-=2*N-2, Y-=4*C)
            {
                for (size_t c=0u; c<C; ++c, X+=2*R-1, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
        else
        {
            Y += 2*(R-1);
            for (size_t r=0u; r<R; ++r, Y-=2*N+2)
            {
                for (size_t c=0u; c<C; ++c, ++X, Y+=2*R-1) { *Y = *X; *++Y = *++X; }
            }
        }
    }
    else if (k==2)
    {
        Y += 2*(N-1);
        for (size_t n=0u; n<N; ++n, ++X, Y-=2) { *Y = *X; *(Y+1) = *++X; }
    }
    else // (k==3)
    {
        if (iscolmajor)
        {
            Y += 2*(C-1);
            for (size_t c=0u; c<C; ++c, Y-=2*N+2)
            {
                for (size_t r=0u; r<R; ++r, ++X, Y+=2*C-1) { *Y = *X; *++Y = *++X; }
            }
        }
        else
        {
            Y += 2*(N-R);
            for (size_t c=0u; c<C; ++c, X-=2*N-2, Y-=4*R)
            {
                for (size_t r=0u; r<R; ++r, X+=2*C-1, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
