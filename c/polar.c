//2-input elementwise function.
//Makes complex output Y from real-valued magnitude (X1) and phase (X2) parts.
//Magnitude is the modulus (abs) and phase is the argument (arg), both real-valued.
//This has no in-place version.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int polar_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int polar_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);


int polar_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;

    if (N1==1)
    {
        for (size_t n=0, n2=0; n<N; n++, n2+=2) { Y[n2] = X1[0]*cosf(X2[n]); Y[n2+1] = X1[0]*sinf(X2[n]); }
    }    
    else if (N2==1)
    {
        const float cosx = cosf(X2[0]), sinx = sinf(X2[0]);
        for (size_t n=0, n2=0; n<N; n++, n2+=2) { Y[n2] = X1[n]*cosx; Y[n2+1] = X1[n]*sinx; }
        //cblas_scopy((int)N,X1,1,Y,2); cblas_scopy((int)N,X1,1,&Y[1],2);
        //cblas_sscal((int)N,cosx,Y,2); cblas_sscal((int)N,sinx,&Y[1],2);
    }
    else if (N==N1 && N==N2)
    {
        for (size_t n=0, n2=0; n<N; n++, n2+=2) { Y[n2] = X1[n]*cosf(X2[n]); Y[n2+1] = X1[n]*sinf(X2[n]); }
    }
    else if (iscolmajor)
    {
        const size_t RCSH1 = R1*C1*S1*(H1>1), RCSH2 = R2*C2*S2*(H2>1);
        const size_t RCS1 = R1*C1*(S1>1), RCS2 = R2*C2*(S2>1);
        const size_t RC1 = R1*(C1>1), RC2 = R2*(C2>1);
        const size_t r1i = (R1>1), r2i = (R2>1);
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                for (size_t c=0; c<C; c++)
                {
                    size_t n1 = h*RCSH1 + s*RCS1 + c*RC1;
                    size_t n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (size_t r=0; r<R; r++, n+=2, n1+=r1i, n2+=r2i)
                    {
                        Y[n] = X1[n1]*cosf(X2[n2]); Y[n+1] = X1[n1]*sinf(X2[n2]);
                    }
                }
            }
        }
    }
    else
    {
        const size_t HSCR1 = H1*S1*C1*(R1>1), HSCR2 = H2*S2*C2*(R2>1);
        const size_t HSC1 = H1*S1*(C1>1), HSC2 = H2*S2*(C2>1);
        const size_t HS1 = H1*(S1>1), HS2 = H2*(S2>1);
        const size_t h1i = (H1>1), h2i = (H2>1);
        for (size_t r=0, n=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++)
            {
                for (size_t s=0; s<S; s++)
                {
                    size_t n1 = r*HSCR1 + c*HSC1 + s*HS1;
                    size_t n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (size_t h=0; h<H; h++, n+=2, n1+=h1i, n2+=h2i)
                    {
                        Y[n] = X1[n1]*cosf(X2[n2]); Y[n+1] = X1[n1]*sinf(X2[n2]);
                    }
                }
            }
        }
    }

    return 0;
}


int polar_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;

    if (N1==1)
    {
        for (size_t n=0, n2=0; n<N; n++, n2+=2) { Y[n2] = X1[0]*cos(X2[n]); Y[n2+1] = X1[0]*sin(X2[n]); }
    }    
    else if (N2==1)
    {
        const double cosx = cos(X2[0]), sinx = sin(X2[0]);
        for (size_t n=0, n2=0; n<N; n++, n2+=2) { Y[n2] = X1[n]*cosx; Y[n2+1] = X1[n]*sinx; }
    }
    else if (N==N1 && N==N2)
    {
        for (size_t n=0, n2=0; n<N; n++, n2+=2) { Y[n2] = X1[n]*cos(X2[n]); Y[n2+1] = X1[n]*sin(X2[n]); }
    }
    else if (iscolmajor)
    {
        const size_t RCSH1 = R1*C1*S1*(H1>1), RCSH2 = R2*C2*S2*(H2>1);
        const size_t RCS1 = R1*C1*(S1>1), RCS2 = R2*C2*(S2>1);
        const size_t RC1 = R1*(C1>1), RC2 = R2*(C2>1);
        const size_t r1i = (R1>1), r2i = (R2>1);
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                for (size_t c=0; c<C; c++)
                {
                    size_t n1 = h*RCSH1 + s*RCS1 + c*RC1;
                    size_t n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (size_t r=0; r<R; r++, n+=2, n1+=r1i, n2+=r2i)
                    {
                        Y[n] = X1[n1]*cos(X2[n2]); Y[n+1] = X1[n1]*sin(X2[n2]);
                    }
                }
            }
        }
    }
    else
    {
        const size_t HSCR1 = H1*S1*C1*(R1>1), HSCR2 = H2*S2*C2*(R2>1);
        const size_t HSC1 = H1*S1*(C1>1), HSC2 = H2*S2*(C2>1);
        const size_t HS1 = H1*(S1>1), HS2 = H2*(S2>1);
        const size_t h1i = (H1>1), h2i = (H2>1);
        for (size_t r=0, n=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++)
            {
                for (size_t s=0; s<S; s++)
                {
                    size_t n1 = r*HSCR1 + c*HSC1 + s*HS1;
                    size_t n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (size_t h=0; h<H; h++, n+=2, n1+=h1i, n2+=h2i)
                    {
                        Y[n] = X1[n1]*cos(X2[n2]); Y[n+1] = X1[n1]*sin(X2[n2]);
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
