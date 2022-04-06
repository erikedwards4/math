//2-input elementwise function.
//Makes complex output Y from real (X1) and imaginary (X2) parts.
//This has no in-place version.

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int complex_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;

    if (N1==1u)
    {
        for (size_t n=N; n>0u; --n, ++X2, ++Y) { *Y = *X1; *++Y = *X2; }
    }    
    else if (N2==1u)
    {
        for (size_t n=N; n>0u; --n, ++X1, ++Y) { *Y = *X1; *++Y = *X2; }
    }
    else if (N1==N2)
    {
        for (size_t n=N; n>0u; --n, ++X1, ++X2, ++Y) { *Y = *X1; *++Y = *X2; }
    }
    else if (iscolmajor)
    {
        const int r1i = (R1>1u), r2i = (R2>1u);
        const int c1i = (int)R1*((C1>1u)-(R1>1u)), c2i = (int)R2*((C2>1u)-(R2>1u));
        const int s1i = (int)(R1*C1)*((S1>1u)-(C1>1u)), s2i = (int)(R2*C2)*((S2>1u)-(C2>1u));
        const int h1i = (int)(R1*C1*S1)*((H1>1u)-(S1>1u)), h2i = (int)(R2*C2*S2)*((H2>1u)-(S2>1u));
        for (size_t h=H; h>0u; --h, X1+=h1i, X2+=h2i)
        {
            for (size_t s=S; s>0u; --s, X1+=s1i, X2+=s2i)
            {
                for (size_t c=C; c>0u; --c, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=R; r>0u; --r, X1+=r1i, X2+=r2i, ++Y)
                    {
                        *Y = *X1; *++Y = *X2;
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = (H1>1u), h2i = (H2>1u);
        const int s1i = (int)H1*((S1>1u)-(H1>1u)), s2i = (int)H2*((S2>1u)-(H2>1u));
        const int c1i = (int)(H1*S1)*((C1>1u)-(S1>1u)), c2i = (int)(H2*S2)*((C2>1u)-(S2>1u));
        const int r1i = (int)(H1*S1*C1)*((R1>1u)-(C1>1u)), r2i = (int)(H2*S2*C2)*((R2>1u)-(C2>1u));
        for (size_t r=R; r>0u; --r, X1+=r1i, X2+=r2i)
        {
            for (size_t c=C; c>0u; --c, X1+=c1i, X2+=c2i)
            {
                for (size_t s=S; s>0u; --s, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=H; h>0u; --h, X1+=h1i, X2+=h2i, ++Y)
                    {
                        *Y = *X1; *++Y = *X2;
                    }
                }
            }
        }
    }

    return 0;
}


int complex_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;   

    if (N1==1u)
    {
        for (size_t n=N; n>0u; --n, ++X2, ++Y) { *Y = *X1; *++Y = *X2; }
    }    
    else if (N2==1u)
    {
        for (size_t n=N; n>0u; --n, ++X1, ++Y) { *Y = *X1; *++Y = *X2; }
    }
    else if (N1==N2)
    {
        for (size_t n=N; n>0u; --n, ++X1, ++X2, ++Y) { *Y = *X1; *++Y = *X2; }
    }
    else if (iscolmajor)
    {
        const int r1i = (R1>1u), r2i = (R2>1u);
        const int c1i = (int)R1*((C1>1u)-(R1>1u)), c2i = (int)R2*((C2>1u)-(R2>1u));
        const int s1i = (int)(R1*C1)*((S1>1u)-(C1>1u)), s2i = (int)(R2*C2)*((S2>1u)-(C2>1u));
        const int h1i = (int)(R1*C1*S1)*((H1>1u)-(S1>1u)), h2i = (int)(R2*C2*S2)*((H2>1u)-(S2>1u));
        for (size_t h=H; h>0u; --h, X1+=h1i, X2+=h2i)
        {
            for (size_t s=S; s>0u; --s, X1+=s1i, X2+=s2i)
            {
                for (size_t c=C; c>0u; --c, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=R; r>0u; --r, X1+=r1i, X2+=r2i, ++Y)
                    {
                        *Y = *X1; *++Y = *X2;
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = (H1>1u), h2i = (H2>1u);
        const int s1i = (int)H1*((S1>1u)-(H1>1u)), s2i = (int)H2*((S2>1u)-(H2>1u));
        const int c1i = (int)(H1*S1)*((C1>1u)-(S1>1u)), c2i = (int)(H2*S2)*((C2>1u)-(S2>1u));
        const int r1i = (int)(H1*S1*C1)*((R1>1u)-(C1>1u)), r2i = (int)(H2*S2*C2)*((R2>1u)-(C2>1u));
        for (size_t r=R; r>0u; --r, X1+=r1i, X2+=r2i)
        {
            for (size_t c=C; c>0u; --c, X1+=c1i, X2+=c2i)
            {
                for (size_t s=S; s>0u; --s, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=H; h>0u; --h, X1+=h1i, X2+=h2i, ++Y)
                    {
                        *Y = *X1; *++Y = *X2;
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
