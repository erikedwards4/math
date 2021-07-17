//Linear algebra function.
//Gets Kronecker matrix product for input matrices X1 and X2.
//For 3D and 4D tensors, this is a generalized Kronecker tensor product.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int kronecker_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int kronecker_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int kronecker_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int kronecker_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);


int kronecker_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = R1*R2, C = C1*C2, S = S1*S2, H = H1*H2;
    const size_t N = R*C*S*H, N2 = R2*C2*S2*H2;

    if (N==0u) {}
    else if (S*H==1u)
    {
        if (iscolmajor)
        {
            for (size_t c1=0u; c1<C1; ++c1)
            {
                for (size_t c2=0u; c2<C2; ++c2)
                {
                    for (size_t r1=0u; r1<R1; ++r1, ++X1)
                    {
                        for (size_t r2=0u; r2<R2; ++r2, ++X2, ++Y)
                        {
                            *Y = *X1 * *X2;
                        }
                        if (r1<R1-1u) { X2 -= R2; }
                    }
                    if (c2<C2-1u) { X1 -= R1; }
                }
                if (c1<C1-1u) { X2 -= N2; }
            }
        }
        else
        {
            for (size_t r1=0u; r1<R1; ++r1)
            {
                for (size_t r2=0u; r2<R2; ++r2)
                {
                    for (size_t c1=0u; c1<C1; ++c1, ++X1)
                    {
                        for (size_t c2=0u; c2<C2; ++c2, ++X2, ++Y)
                        {
                            *Y = *X1 * *X2;
                        }
                        if (c1<C1-1u) { X2 -= C2; }
                    }
                    if (r2<R2-1u) { X1 -= C1; }
                }
                if (r1<R1-1u) { X2 -= N2; }
            }
        }
    }
    else if (iscolmajor)
    {
        for (size_t h1=0; h1<H1; ++h1, X1+=R1*C1*S1, X2-=N2)
        {
            for (size_t h2=0; h2<H2; ++h2, X1-=R1*C1*S1, X2+=R2*C2*S2)
            {
                for (size_t s1=0; s1<S1; ++s1, X1+=R1*C1, X2-=R2*C2*S2)
                {
                    for (size_t s2=0; s2<S2; ++s2, X1-=R1*C1, X2+=R2*C2)
                    {
                        for (size_t c1=0u; c1<C1; ++c1, X1+=R1, X2-=R2*C2)
                        {
                            for (size_t c2=0u; c2<C2; ++c2, X1-=R1, X2+=R2)
                            {
                                for (size_t r1=0u; r1<R1; ++r1, ++X1, X2-=R2)
                                {
                                    for (size_t r2=0u; r2<R2; ++r2, ++X2, ++Y)
                                    {
                                        *Y = *X1 * *X2;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        for (size_t r1=0u; r1<R1; ++r1, X1+=C1*S1*H1, X2-=N2)
        {
            for (size_t r2=0u; r2<R2; ++r2, X1-=C1*S1*H1, X2+=C2*S2*H2)
            {
                for (size_t c1=0u; c1<C1; ++c1, X1+=S1*H1, X2-=C2*S2*H2)
                {
                    for (size_t c2=0u; c2<C2; ++c2, X1-=S1*H1, X2+=S2*H2)
                    {
                        for (size_t s1=0; s1<S1; ++s1, X1+=H1, X2-=S2*H2)
                        {
                            for (size_t s2=0; s2<S2; ++s2, X1-=H1, X2+=H2)
                            {
                                for (size_t h1=0; h1<H1; ++h1, ++X1, X2-=H2)
                                {
                                    for (size_t h2=0; h2<H2; ++h2, ++X2, ++Y)
                                    {
                                        *Y = *X1 * *X2;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    return 0;
}


int kronecker_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = R1*R2, C = C1*C2, S = S1*S2, H = H1*H2;
    const size_t N = R*C*S*H, N2 = R2*C2*S2*H2;

    if (N==0u) {}
    else if (S*H==1u)
    {
        if (iscolmajor)
        {
            for (size_t c1=0u; c1<C1; ++c1)
            {
                for (size_t c2=0u; c2<C2; ++c2)
                {
                    for (size_t r1=0u; r1<R1; ++r1, ++X1)
                    {
                        for (size_t r2=0u; r2<R2; ++r2, ++X2, ++Y)
                        {
                            *Y = *X1 * *X2;
                        }
                        if (r1<R1-1u) { X2 -= R2; }
                    }
                    if (c2<C2-1u) { X1 -= R1; }
                }
                if (c1<C1-1u) { X2 -= N2; }
            }
        }
        else
        {
            for (size_t r1=0u; r1<R1; ++r1)
            {
                for (size_t r2=0u; r2<R2; ++r2)
                {
                    for (size_t c1=0u; c1<C1; ++c1, ++X1)
                    {
                        for (size_t c2=0u; c2<C2; ++c2, ++X2, ++Y)
                        {
                            *Y = *X1 * *X2;
                        }
                        if (c1<C1-1u) { X2 -= C2; }
                    }
                    if (r2<R2-1u) { X1 -= C1; }
                }
                if (r1<R1-1u) { X2 -= N2; }
            }
        }
    }
    else if (iscolmajor)
    {
        for (size_t h1=0; h1<H1; ++h1, X1+=R1*C1*S1, X2-=N2)
        {
            for (size_t h2=0; h2<H2; ++h2, X1-=R1*C1*S1, X2+=R2*C2*S2)
            {
                for (size_t s1=0; s1<S1; ++s1, X1+=R1*C1, X2-=R2*C2*S2)
                {
                    for (size_t s2=0; s2<S2; ++s2, X1-=R1*C1, X2+=R2*C2)
                    {
                        for (size_t c1=0u; c1<C1; ++c1, X1+=R1, X2-=R2*C2)
                        {
                            for (size_t c2=0u; c2<C2; ++c2, X1-=R1, X2+=R2)
                            {
                                for (size_t r1=0u; r1<R1; ++r1, ++X1, X2-=R2)
                                {
                                    for (size_t r2=0u; r2<R2; ++r2, ++X2, ++Y)
                                    {
                                        *Y = *X1 * *X2;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        for (size_t r1=0u; r1<R1; ++r1, X1+=C1*S1*H1, X2-=N2)
        {
            for (size_t r2=0u; r2<R2; ++r2, X1-=C1*S1*H1, X2+=C2*S2*H2)
            {
                for (size_t c1=0u; c1<C1; ++c1, X1+=S1*H1, X2-=C2*S2*H2)
                {
                    for (size_t c2=0u; c2<C2; ++c2, X1-=S1*H1, X2+=S2*H2)
                    {
                        for (size_t s1=0; s1<S1; ++s1, X1+=H1, X2-=S2*H2)
                        {
                            for (size_t s2=0; s2<S2; ++s2, X1-=H1, X2+=H2)
                            {
                                for (size_t h1=0; h1<H1; ++h1, ++X1, X2-=H2)
                                {
                                    for (size_t h2=0; h2<H2; ++h2, ++X2, ++Y)
                                    {
                                        *Y = *X1 * *X2;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}


int kronecker_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = R1*R2, C = C1*C2, S = S1*S2, H = H1*H2;
    const size_t N = R*C*S*H, N2 = R2*C2*S2*H2;
    float x1r, x1i, x2r, x2i;

    if (N==0u) {}
    else if (S*H==1u)
    {
        if (iscolmajor)
        {
            for (size_t c1=0u; c1<C1; ++c1)
            {
                for (size_t c2=0u; c2<C2; ++c2)
                {
                    for (size_t r1=0u; r1<R1; ++r1, ++X1)
                    {
                        x1r = *X1; x1i = *++X1;
                        for (size_t r2=0u; r2<R2; ++r2, ++X2, ++Y)
                        {
                            x2r = *X2; x2i = *++X2;
                            *Y = x1r*x2r - x1i*x2i;
                            *++Y = x1r*x2i + x1i*x2r;
                        }
                        if (r1<R1-1u) { X2 -= 2u*R2; }
                    }
                    if (c2<C2-1u) { X1 -= 2u*R1; }
                }
                if (c1<C1-1u) { X2 -= 2u*N2; }
            }
        }
        else
        {
            for (size_t r1=0u; r1<R1; ++r1)
            {
                for (size_t r2=0u; r2<R2; ++r2)
                {
                    for (size_t c1=0u; c1<C1; ++c1, ++X1)
                    {
                        x1r = *X1; x1i = *++X1;
                        for (size_t c2=0u; c2<C2; ++c2, ++X2, ++Y)
                        {
                            x2r = *X2; x2i = *++X2;
                            *Y = x1r*x2r - x1i*x2i;
                            *++Y = x1r*x2i + x1i*x2r;
                        }
                        if (c1<C1-1u) { X2 -= 2u*C2; }
                    }
                    if (r2<R2-1u) { X1 -= 2u*C1; }
                }
                if (r1<R1-1u) { X2 -= 2u*N2; }
            }
        }
    }
    else if (iscolmajor)
    {
        for (size_t h1=0; h1<H1; ++h1, X1+=2u*R1*C1*S1, X2-=2u*N2)
        {
            for (size_t h2=0; h2<H2; ++h2, X1-=2u*R1*C1*S1, X2+=2u*R2*C2*S2)
            {
                for (size_t s1=0; s1<S1; ++s1, X1+=2u*R1*C1, X2-=2u*R2*C2*S2)
                {
                    for (size_t s2=0; s2<S2; ++s2, X1-=2u*R1*C1, X2+=2u*R2*C2)
                    {
                        for (size_t c1=0u; c1<C1; ++c1, X1+=2u*R1, X2-=2u*R2*C2)
                        {
                            for (size_t c2=0u; c2<C2; ++c2, X1-=2u*R1, X2+=2u*R2)
                            {
                                for (size_t r1=0u; r1<R1; ++r1, ++X1, X2-=2u*R2)
                                {
                                    x1r = *X1; x1i = *++X1;
                                    for (size_t r2=0u; r2<R2; ++r2, ++X2, ++Y)
                                    {
                                        x2r = *X2; x2i = *++X2;
                                        *Y = x1r*x2r - x1i*x2i;
                                        *++Y = x1r*x2i + x1i*x2r;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        for (size_t r1=0u; r1<R1; ++r1, X1+=2u*C1*S1*H1, X2-=2u*N2)
        {
            for (size_t r2=0u; r2<R2; ++r2, X1-=2u*C1*S1*H1, X2+=2u*C2*S2*H2)
            {
                for (size_t c1=0u; c1<C1; ++c1, X1+=2u*S1*H1, X2-=2u*C2*S2*H2)
                {
                    for (size_t c2=0u; c2<C2; ++c2, X1-=2u*S1*H1, X2+=2u*S2*H2)
                    {
                        for (size_t s1=0; s1<S1; ++s1, X1+=2u*H1, X2-=2u*S2*H2)
                        {
                            for (size_t s2=0; s2<S2; ++s2, X1-=2u*H1, X2+=2u*H2)
                            {
                                for (size_t h1=0; h1<H1; ++h1, ++X1, X2-=2u*H2)
                                {
                                    x1r = *X1; x1i = *++X1;
                                    for (size_t h2=0; h2<H2; ++h2, ++X2, ++Y)
                                    {
                                        x2r = *X2; x2i = *++X2;
                                        *Y = x1r*x2r - x1i*x2i;
                                        *++Y = x1r*x2i + x1i*x2r;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    return 0;
}


int kronecker_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = R1*R2, C = C1*C2, S = S1*S2, H = H1*H2;
    const size_t N = R*C*S*H, N2 = R2*C2*S2*H2;
    double x1r, x1i, x2r, x2i;

    if (N==0u) {}
    else if (S*H==1u)
    {
        if (iscolmajor)
        {
            for (size_t c1=0u; c1<C1; ++c1)
            {
                for (size_t c2=0u; c2<C2; ++c2)
                {
                    for (size_t r1=0u; r1<R1; ++r1, ++X1)
                    {
                        x1r = *X1; x1i = *++X1;
                        for (size_t r2=0u; r2<R2; ++r2, ++X2, ++Y)
                        {
                            x2r = *X2; x2i = *++X2;
                            *Y = x1r*x2r - x1i*x2i;
                            *++Y = x1r*x2i + x1i*x2r;
                        }
                        if (r1<R1-1u) { X2 -= 2u*R2; }
                    }
                    if (c2<C2-1u) { X1 -= 2u*R1; }
                }
                if (c1<C1-1u) { X2 -= 2u*N2; }
            }
        }
        else
        {
            for (size_t r1=0u; r1<R1; ++r1)
            {
                for (size_t r2=0u; r2<R2; ++r2)
                {
                    for (size_t c1=0u; c1<C1; ++c1, ++X1)
                    {
                        x1r = *X1; x1i = *++X1;
                        for (size_t c2=0u; c2<C2; ++c2, ++X2, ++Y)
                        {
                            x2r = *X2; x2i = *++X2;
                            *Y = x1r*x2r - x1i*x2i;
                            *++Y = x1r*x2i + x1i*x2r;
                        }
                        if (c1<C1-1u) { X2 -= 2u*C2; }
                    }
                    if (r2<R2-1u) { X1 -= 2u*C1; }
                }
                if (r1<R1-1u) { X2 -= 2u*N2; }
            }
        }
    }
    else if (iscolmajor)
    {
        for (size_t h1=0; h1<H1; ++h1, X1+=2u*R1*C1*S1, X2-=2u*N2)
        {
            for (size_t h2=0; h2<H2; ++h2, X1-=2u*R1*C1*S1, X2+=2u*R2*C2*S2)
            {
                for (size_t s1=0; s1<S1; ++s1, X1+=2u*R1*C1, X2-=2u*R2*C2*S2)
                {
                    for (size_t s2=0; s2<S2; ++s2, X1-=2u*R1*C1, X2+=2u*R2*C2)
                    {
                        for (size_t c1=0u; c1<C1; ++c1, X1+=2u*R1, X2-=2u*R2*C2)
                        {
                            for (size_t c2=0u; c2<C2; ++c2, X1-=2u*R1, X2+=2u*R2)
                            {
                                for (size_t r1=0u; r1<R1; ++r1, ++X1, X2-=2u*R2)
                                {
                                    x1r = *X1; x1i = *++X1;
                                    for (size_t r2=0u; r2<R2; ++r2, ++X2, ++Y)
                                    {
                                        x2r = *X2; x2i = *++X2;
                                        *Y = x1r*x2r - x1i*x2i;
                                        *++Y = x1r*x2i + x1i*x2r;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        for (size_t r1=0u; r1<R1; ++r1, X1+=2u*C1*S1*H1, X2-=2u*N2)
        {
            for (size_t r2=0u; r2<R2; ++r2, X1-=2u*C1*S1*H1, X2+=2u*C2*S2*H2)
            {
                for (size_t c1=0u; c1<C1; ++c1, X1+=2u*S1*H1, X2-=2u*C2*S2*H2)
                {
                    for (size_t c2=0u; c2<C2; ++c2, X1-=2u*S1*H1, X2+=2u*S2*H2)
                    {
                        for (size_t s1=0; s1<S1; ++s1, X1+=2u*H1, X2-=2u*S2*H2)
                        {
                            for (size_t s2=0; s2<S2; ++s2, X1-=2u*H1, X2+=2u*H2)
                            {
                                for (size_t h1=0; h1<H1; ++h1, ++X1, X2-=2u*H2)
                                {
                                    x1r = *X1; x1i = *++X1;
                                    for (size_t h2=0; h2<H2; ++h2, ++X2, ++Y)
                                    {
                                        x2r = *X2; x2i = *++X2;
                                        *Y = x1r*x2r - x1i*x2i;
                                        *++Y = x1r*x2i + x1i*x2r;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
