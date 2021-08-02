//2-input elementwise function.
//Makes complex output Y from real-valued magnitude (X1) and phase (X2) parts.
//Magnitude is the modulus (abs) and phase is the argument (arg), both real-valued.
//This has no in-place version.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int polar_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);
int polar_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor);


int polar_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;

    if (N1==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X2, ++Y) { *Y = *X1 * cosf(*X2); *++Y = *X1 * sinf(*X2); }
        //for (size_t n=0, n2=0; n<N; ++n, n2+=2) { Y[n2] = *X1 * cosf(*X2); Y[n2+1] = *X1 * sinf(*X2); }
    }
    else if (N2==1u)
    {
        const float cosx = cosf(*X2), sinx = sinf(*X2);
        for (size_t n=0u; n<N; ++n, ++X1, ++Y) { *Y = *X1 * cosx; *++Y = *X1 * sinx; }
        //for (size_t n=0, n2=0; n<N; ++n, n2+=2) { Y[n2] = *X1 * cosx; Y[n2+1] = *X1 * sinx; }
        //cblas_scopy((int)N,X1,1,Y,2); cblas_scopy((int)N,X1,1,&Y[1],2);
        //cblas_sscal((int)N,cosx,Y,2); cblas_sscal((int)N,sinx,&Y[1],2);
    }
    else if (N1==N2)
    {
        for (size_t n=0u; n<N; ++n, ++X1, ++X2, ++Y) { *Y = *X1 * cosf(*X2); *++Y = *X1 * sinf(*X2); }
        //for (size_t n=0, n2=0; n<N; ++n, n2+=2) { Y[n2] = *X1 * cosf(*X2); Y[n2+1] = *X1 * sinf(*X2); }
    }
    else if (iscolmajor)
    {
        const int r1i = (int)(R1>1u), r2i = (int)(R2>1u);
        const int c1i = (int)R1*((int)(C1>1u)-(int)(R1>1u)), c2i = (int)R2*((int)(C2>1u)-(int)(R2>1u));
        const int s1i = (int)(R1*C1)*((int)(S1>1u)-(int)(C1>1u)), s2i = (int)(R2*C2)*((int)(S2>1u)-(int)(C2>1u));
        const int h1i = (int)(R1*C1*S1)*((int)(H1>1u)-(int)(S1>1u)), h2i = (int)(R2*C2*S2)*((int)(H2>1u)-(int)(S2>1u));
        for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i, ++Y)
                    {
                        *Y = *X1 * cosf(*X2);
                        *++Y = *X1 * sinf(*X2);
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = (int)(H1>1u), h2i = (int)(H2>1u);
        const int s1i = (int)H1*((int)(S1>1u)-(int)(H1>1u)), s2i = (int)H2*((int)(S2>1u)-(int)(H2>1u));
        const int c1i = (int)(H1*S1)*((int)(C1>1u)-(int)(S1>1u)), c2i = (int)(H2*S2)*((int)(C2>1u)-(int)(S2>1u));
        const int r1i = (int)(H1*S1*C1)*((int)(R1>1u)-(int)(C1>1u)), r2i = (int)(H2*S2*C2)*((int)(R2>1u)-(int)(C2>1u));
        for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i, ++Y)
                    {
                        *Y = *X1 * cosf(*X2);
                        *++Y = *X1 * sinf(*X2);
                    }
                }
            }
        }
    }

    return 0;
}


int polar_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;

    if (N1==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X2, ++Y) { *Y = *X1 * cos(*X2); *++Y = *X1 * sin(*X2); }
    }    
    else if (N2==1u)
    {
        const double cosx = cos(*X2), sinx = sin(*X2);
        for (size_t n=0u; n<N; ++n, ++X1, ++Y) { *Y = *X1 * cosx; *++Y = *X1 * sinx; }
    }
    else if (N1==N2)
    {
        for (size_t n=0u; n<N; ++n, ++X1, ++X2, ++Y) { *Y = *X1 * cos(*X2); *++Y = *X1 * sin(*X2); }
    }
    else if (iscolmajor)
    {
        const int r1i = (int)(R1>1u), r2i = (int)(R2>1u);
        const int c1i = (int)R1*((int)(C1>1u)-(int)(R1>1u)), c2i = (int)R2*((int)(C2>1u)-(int)(R2>1u));
        const int s1i = (int)(R1*C1)*((int)(S1>1u)-(int)(C1>1u)), s2i = (int)(R2*C2)*((int)(S2>1u)-(int)(C2>1u));
        const int h1i = (int)(R1*C1*S1)*((int)(H1>1u)-(int)(S1>1u)), h2i = (int)(R2*C2*S2)*((int)(H2>1u)-(int)(S2>1u));
        for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i, ++Y)
                    {
                        *Y = *X1 * cos(*X2);
                        *++Y = *X1 * sin(*X2);
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = (int)(H1>1u), h2i = (int)(H2>1u);
        const int s1i = (int)H1*((int)(S1>1u)-(int)(H1>1u)), s2i = (int)H2*((int)(S2>1u)-(int)(H2>1u));
        const int c1i = (int)(H1*S1)*((int)(C1>1u)-(int)(S1>1u)), c2i = (int)(H2*S2)*((int)(C2>1u)-(int)(S2>1u));
        const int r1i = (int)(H1*S1*C1)*((int)(R1>1u)-(int)(C1>1u)), r2i = (int)(H2*S2*C2)*((int)(R2>1u)-(int)(C2>1u));
        for (size_t r=0u; r<R; ++r, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0u; c<C; ++c, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0u; s<S; ++s, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0u; h<H; ++h, X1+=h1i, X2+=h2i, ++Y)
                    {
                        *Y = *X1 * cos(*X2);
                        *++Y = *X1 * sin(*X2);
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
