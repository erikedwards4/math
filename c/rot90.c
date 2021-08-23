//Rotates matrix X by 90 degrees K times.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int rot90_s (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int K);
int rot90_d (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int K);
int rot90_c (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int K);
int rot90_z (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int K);


int rot90_s (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int K)
{
    const size_t N = R*C;
    const int k = (K<0) ? 4+K%4 : K%4;

    if (N==0u) {}
    else if (k==0 || N==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (k==1)
    {
        if (iscolmajor)
        {
            Y += N - C;
            for (size_t r=R; r>0u; --r, X-=N-1u, Y-=2u*C)
            {
                for (size_t c=C; c>0u; --c, X+=R, ++Y) { *Y = *X; }
            }
        }
        else
        {
            Y += R - 1u;
            for (size_t r=R; r>0u; --r, Y-=N+1u)
            {
                for (size_t c=C; c>0u; --c, ++X, Y+=R) { *Y = *X; }
            }
        }
    }
    else if (k==2)
    {
        Y += N - 1u;
        for (size_t n=N; n>0u; --n, ++X, --Y) { *Y = *X; }
    }
    else // (k==3)
    {
        if (iscolmajor)
        {
            Y += C - 1u;
            for (size_t c=C; c>0u; --c, Y-=N+1u)
            {
                for (size_t r=R; r>0u; --r, ++X, Y+=C) { *Y = *X; }
            }
        }
        else
        {
            Y += N - R;
            for (size_t c=C; c>0u; --c, X-=N-1u, Y-=2u*R)
            {
                for (size_t r=R; r>0u; --r, X+=C, ++Y) { *Y = *X; }
            }
        }
    }

    return 0;
}


int rot90_d (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int K)
{
    const size_t N = R*C;
    const int k = (K<0) ? 4+K%4 : K%4;

    if (N==0u) {}
    else if (k==0 || N==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (k==1)
    {
        if (iscolmajor)
        {
            Y += N - C;
            for (size_t r=R; r>0u; --r, X-=N-1u, Y-=2u*C)
            {
                for (size_t c=C; c>0u; --c, X+=R, ++Y) { *Y = *X; }
            }
        }
        else
        {
            Y += R - 1u;
            for (size_t r=R; r>0u; --r, Y-=N+1u)
            {
                for (size_t c=C; c>0u; --c, ++X, Y+=R) { *Y = *X; }
            }
        }
    }
    else if (k==2)
    {
        Y += N - 1u;
        for (size_t n=N; n>0u; --n, ++X, --Y) { *Y = *X; }
    }
    else // (k==3)
    {
        if (iscolmajor)
        {
            Y += C - 1u;
            for (size_t c=C; c>0u; --c, Y-=N+1u)
            {
                for (size_t r=R; r>0u; --r, ++X, Y+=C) { *Y = *X; }
            }
        }
        else
        {
            Y += N - R;
            for (size_t c=C; c>0u; --c, X-=N-1u, Y-=2u*R)
            {
                for (size_t r=R; r>0u; --r, X+=C, ++Y) { *Y = *X; }
            }
        }
    }

    return 0;
}


int rot90_c (float *Y, const float *X, const size_t R, const size_t C, const int iscolmajor, const int K)
{
    const size_t N = R*C;
    const int k = (K<0) ? 4+K%4 : K%4;

    if (N==0u) {}
    else if (k==0 || N==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (k==1)
    {
        if (iscolmajor)
        {
            Y += 2u*(N-C);
            for (size_t r=R; r>0u; --r, X-=2u*N-2u, Y-=4u*C)
            {
                for (size_t c=C; c>0u; --c, X+=2u*R-1u, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
        else
        {
            Y += 2u*(R-1u);
            for (size_t r=R; r>0u; --r, Y-=2u*N+2)
            {
                for (size_t c=C; c>0u; --c, ++X, Y+=2u*R-1u) { *Y = *X; *++Y = *++X; }
            }
        }
    }
    else if (k==2)
    {
        Y += 2u*(N-1u);
        for (size_t n=N; n>0u; --n, ++X, Y-=2) { *Y = *X; *(Y+1) = *++X; }
    }
    else // (k==3)
    {
        if (iscolmajor)
        {
            Y += 2u*(C-1u);
            for (size_t c=C; c>0u; --c, Y-=2u*N+2)
            {
                for (size_t r=R; r>0u; --r, ++X, Y+=2u*C-1u) { *Y = *X; *++Y = *++X; }
            }
        }
        else
        {
            Y += 2u*(N-R);
            for (size_t c=C; c>0u; --c, X-=2u*N-2u, Y-=4u*R)
            {
                for (size_t r=R; r>0u; --r, X+=2u*C-1u, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
    }

    return 0;
}


int rot90_z (double *Y, const double *X, const size_t R, const size_t C, const int iscolmajor, const int K)
{
    const size_t N = R*C;
    const int k = (K<0) ? 4+K%4 : K%4;

    if (N==0u) {}
    else if (k==0 || N==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (k==1)
    {
        if (iscolmajor)
        {
            Y += 2u*(N-C);
            for (size_t r=R; r>0u; --r, X-=2u*N-2u, Y-=4u*C)
            {
                for (size_t c=C; c>0u; --c, X+=2u*R-1u, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
        else
        {
            Y += 2u*(R-1u);
            for (size_t r=R; r>0u; --r, Y-=2u*N+2)
            {
                for (size_t c=C; c>0u; --c, ++X, Y+=2u*R-1u) { *Y = *X; *++Y = *++X; }
            }
        }
    }
    else if (k==2)
    {
        Y += 2u*(N-1u);
        for (size_t n=N; n>0u; --n, ++X, Y-=2) { *Y = *X; *(Y+1) = *++X; }
    }
    else // (k==3)
    {
        if (iscolmajor)
        {
            Y += 2u*(C-1u);
            for (size_t c=C; c>0u; --c, Y-=2u*N+2)
            {
                for (size_t r=R; r>0u; --r, ++X, Y+=2u*C-1u) { *Y = *X; *++Y = *++X; }
            }
        }
        else
        {
            Y += 2u*(N-R);
            for (size_t c=C; c>0u; --c, X-=2u*N-2u, Y-=4u*R)
            {
                for (size_t r=R; r>0u; --r, X+=2u*C-1u, ++Y) { *Y = *X; *++Y = *++X; }
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
