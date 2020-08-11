//2-input elementwise function.
//Makes complex output Y from real (X1) and imaginary (X2) parts.
//This has no in-place version.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int complex_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int complex_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);


int complex_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;

    if (N1==1)
    {
        for (size_t n=0; n<N; ++n, ++X2, ++Y) { *Y = *X1; *++Y = *X2; }
    }    
    else if (N2==1)
    {
        for (size_t n=0; n<N; ++n, ++X1, ++Y) { *Y = *X1; *++Y = *X2; }
    }
    else if (N1==N2)
    {
        for (size_t n=0; n<N; ++n, ++X1, ++X2, ++Y) { *Y = *X1; *++Y = *X2; }
    }
    else if (iscolmajor)
    {
        const int r1i = (int)(R1>1), r2i = (int)(R2>1);
        const int c1i = (int)R1*((int)(C1>1)-(int)(R1>1)), c2i = (int)R2*((int)(C2>1)-(int)(R2>1));
        const int s1i = (int)(R1*C1)*((int)(S1>1)-(int)(C1>1)), s2i = (int)(R2*C2)*((int)(S2>1)-(int)(C2>1));
        const int h1i = (int)(R1*C1*S1)*((int)(H1>1)-(int)(S1>1)), h2i = (int)(R2*C2*S2)*((int)(H2>1)-(int)(S2>1));
        for (size_t h=0; h<H; ++h, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0; s<S; ++s, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0; c<C; ++c, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0; r<R; ++r, X1+=r1i, X2+=r2i, ++Y)
                    {
                        *Y = *X1; *++Y = *X2;
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = (int)(H1>1), h2i = (int)(H2>1);
        const int s1i = (int)H1*((int)(S1>1)-(int)(H1>1)), s2i = (int)H2*((int)(S2>1)-(int)(H2>1));
        const int c1i = (int)(H1*S1)*((int)(C1>1)-(int)(S1>1)), c2i = (int)(H2*S2)*((int)(C2>1)-(int)(S2>1));
        const int r1i = (int)(H1*S1*C1)*((int)(R1>1)-(int)(C1>1)), r2i = (int)(H2*S2*C2)*((int)(R2>1)-(int)(C2>1));
        for (size_t r=0; r<R; ++r, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0; c<C; ++c, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0; s<S; ++s, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0; h<H; ++h, X1+=h1i, X2+=h2i, ++Y)
                    {
                        *Y = *X1; *++Y = *X2;
                    }
                }
            }
        }
    }

    return 0;
}


int complex_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;   

    if (N1==1)
    {
        for (size_t n=0; n<N; ++n, ++X2, ++Y) { *Y = *X1; *++Y = *X2; }
    }    
    else if (N2==1)
    {
        for (size_t n=0; n<N; ++n, ++X1, ++Y) { *Y = *X1; *++Y = *X2; }
    }
    else if (N1==N2)
    {
        for (size_t n=0; n<N; ++n, ++X1, ++X2, ++Y) { *Y = *X1; *++Y = *X2; }
    }
    else if (iscolmajor)
    {
        const int r1i = (int)(R1>1), r2i = (int)(R2>1);
        const int c1i = (int)R1*((int)(C1>1)-(int)(R1>1)), c2i = (int)R2*((int)(C2>1)-(int)(R2>1));
        const int s1i = (int)(R1*C1)*((int)(S1>1)-(int)(C1>1)), s2i = (int)(R2*C2)*((int)(S2>1)-(int)(C2>1));
        const int h1i = (int)(R1*C1*S1)*((int)(H1>1)-(int)(S1>1)), h2i = (int)(R2*C2*S2)*((int)(H2>1)-(int)(S2>1));
        for (size_t h=0; h<H; ++h, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0; s<S; ++s, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0; c<C; ++c, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0; r<R; ++r, X1+=r1i, X2+=r2i, ++Y)
                    {
                        *Y = *X1; *++Y = *X2;
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = (int)(H1>1), h2i = (int)(H2>1);
        const int s1i = (int)H1*((int)(S1>1)-(int)(H1>1)), s2i = (int)H2*((int)(S2>1)-(int)(H2>1));
        const int c1i = (int)(H1*S1)*((int)(C1>1)-(int)(S1>1)), c2i = (int)(H2*S2)*((int)(C2>1)-(int)(S2>1));
        const int r1i = (int)(H1*S1*C1)*((int)(R1>1)-(int)(C1>1)), r2i = (int)(H2*S2*C2)*((int)(R2>1)-(int)(C2>1));
        for (size_t r=0; r<R; ++r, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0; c<C; ++c, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0; s<S; ++s, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0; h<H; ++h, X1+=h1i, X2+=h2i, ++Y)
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
