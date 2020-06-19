//2-input elementwise function.
//Makes complex output Y from real (X1) and imaginary (X2) parts.
//This has no in-place version.

#include <stdio.h>
#include <cblas.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int complex_s (float *Y, const float *X1, const float *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);
int complex_d (double *Y, const double *X1, const double *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);


int complex_s (float *Y, const float *X1, const float *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    int r, c, s, h, n1, n2, n = 0;

    //Checks
    if (R1<0) { fprintf(stderr,"error in complex_s: R1 (num rows X1) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in complex_s: C1 (num cols X1) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in complex_s: S1 (num slices X1) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in complex_s: H1 (num hyperslices X1) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in complex_s: R2 (num rows X2) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in complex_s: C2 (num cols X2) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in complex_s: S2 (num slices X2) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in complex_s: H2 (num hyperslices X2) must be nonnegative\n"); return 1; }

    if (N1==1)
    {
        cblas_scopy(N,X1,0,Y,2);
        cblas_scopy(N,X2,1,&Y[1],2);
    }    
    else if (N2==1)
    {
        cblas_scopy(N,X1,1,Y,2);
        cblas_scopy(N,X2,0,&Y[1],2);
    }
    else if (N==N1 && N==N2)
    {
        cblas_scopy(N,X1,1,Y,2);
        cblas_scopy(N,X2,1,&Y[1],2);
    }
    else if (iscolmajor)
    {
        const int RCSH1 = R1*C1*S1*(H1>1), RCSH2 = R2*C2*S2*(H2>1);
        const int RCS1 = R1*C1*(S1>1), RCS2 = R2*C2*(S2>1);
        const int RC1 = R1*(C1>1), RC2 = R2*(C2>1);
        const int r1i = (R1>1), r2i = (R2>1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                for (c=0; c<C; c++)
                {
                    n1 = h*RCSH1 + s*RCS1 + c*RC1;
                    n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (r=0; r<R; r++)
                    {
                        Y[n] = X1[n1]; Y[n+1] = X2[n2];
                        n+=2; n1 += r1i; n2 += r2i;
                    }
                }
            }
        }
    }
    else
    {
        const int HSCR1 = H1*S1*C1*(R1>1), HSCR2 = H2*S2*C2*(R2>1);
        const int HSC1 = H1*S1*(C1>1), HSC2 = H2*S2*(C2>1);
        const int HS1 = H1*(S1>1), HS2 = H2*(S2>1);
        const int h1i = (H1>1), h2i = (H2>1);
        for (r=0; r<R; r++)
        {
            for (c=0; c<C; c++)
            {
                for (s=0; s<S; s++)
                {
                    n1 = r*HSCR1 + c*HSC1 + s*HS1;
                    n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (h=0; h<H; h++)
                    {
                        Y[n] = X1[n1]; Y[n+1] = X2[n2];
                        n+=2; n1 += h1i; n2 += h2i;
                    }
                }
            }
        }
    }

    return 0;
}


int complex_d (double *Y, const double *X1, const double *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    int r, c, s, h, n1, n2, n = 0;

    //Checks
    if (R1<0) { fprintf(stderr,"error in complex_d: R1 (num rows X1) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in complex_d: C1 (num cols X1) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in complex_d: S1 (num slices X1) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in complex_d: H1 (num hyperslices X1) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in complex_d: R2 (num rows X2) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in complex_d: C2 (num cols X2) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in complex_d: S2 (num slices X2) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in complex_d: H2 (num hyperslices X2) must be nonnegative\n"); return 1; }

    if (N1==1)
    {
        cblas_dcopy(N,X1,0,Y,2);
        cblas_dcopy(N,X2,1,&Y[1],2);
    }    
    else if (N2==1)
    {
        cblas_dcopy(N,X1,1,Y,2);
        cblas_dcopy(N,X2,0,&Y[1],2);
    }
    else if (N==N1 && N==N2)
    {
        cblas_dcopy(N,X1,1,Y,2);
        cblas_dcopy(N,X2,1,&Y[1],2);
    }
    else if (iscolmajor)
    {
        const int RCSH1 = R1*C1*S1*(H1>1), RCSH2 = R2*C2*S2*(H2>1);
        const int RCS1 = R1*C1*(S1>1), RCS2 = R2*C2*(S2>1);
        const int RC1 = R1*(C1>1), RC2 = R2*(C2>1);
        const int r1i = (R1>1), r2i = (R2>1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                for (c=0; c<C; c++)
                {
                    n1 = h*RCSH1 + s*RCS1 + c*RC1;
                    n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (r=0; r<R; r++)
                    {
                        Y[n] = X1[n1]; Y[n+1] = X2[n2];
                        n+=2; n1 += r1i; n2 += r2i;
                    }
                }
            }
        }
    }
    else
    {
        const int HSCR1 = H1*S1*C1*(R1>1), HSCR2 = H2*S2*C2*(R2>1);
        const int HSC1 = H1*S1*(C1>1), HSC2 = H2*S2*(C2>1);
        const int HS1 = H1*(S1>1), HS2 = H2*(S2>1);
        const int h1i = (H1>1), h2i = (H2>1);
        for (r=0; r<R; r++)
        {
            for (c=0; c<C; c++)
            {
                for (s=0; s<S; s++)
                {
                    n1 = r*HSCR1 + c*HSC1 + s*HS1;
                    n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (h=0; h<H; h++)
                    {
                        Y[n] = X1[n1]; Y[n+1] = X2[n2];
                        n+=2; n1 += h1i; n2 += h2i;
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
