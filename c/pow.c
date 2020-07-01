//2-input elementwise function.
//Raises each element of X to the power of P.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int pow_s (float *Y, const float *X, const float *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);
int pow_d (double *Y, const double *X, const double *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);
int pow_c (float *Y, const float *X, const float *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);
int pow_z (double *Y, const double *X, const double *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);

int pow_inplace_s (float *X, const float *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);
int pow_inplace_d (double *X, const double *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);
int pow_inplace_c (float *X, const float *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);
int pow_inplace_z (double *X, const double *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);


int pow_s (float *Y, const float *X, const float *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    if (R1<0) { fprintf(stderr,"error in pow_s: R1 (num rows X) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in pow_s: C1 (num cols X) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in pow_s: S1 (num slices X) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in pow_s: H1 (num hyperslices X) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in pow_s: R2 (num rows P) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in pow_s: C2 (num cols P) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in pow_s: S2 (num slices P) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in pow_s: H2 (num hyperslices P) must be nonnegative\n"); return 1; }

    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    int n1, n2, n = 0;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    if (N1==1)
    {
        while (n<N) { Y[n] = powf(X[0],P[n]); n++; }
    }    
    else if (N2==1)
    {
        while (n<N) { Y[n] = powf(X[n],P[0]); n++; }
    }
    else if (N==N1 && N==N2)
    {
        while (n<N) { Y[n] = powf(X[n],P[n]); n++; }
    }
    else if (iscolmajor)
    {
        const int RCSH1 = R1*C1*S1*(H1>1), RCSH2 = R2*C2*S2*(H2>1);
        const int RCS1 = R1*C1*(S1>1), RCS2 = R2*C2*(S2>1);
        const int RC1 = R1*(C1>1), RC2 = R2*(C2>1);
        const int r1i = (R1>1), r2i = (R2>1);
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                for (int c=0; c<C; c++)
                {
                    n1 = h*RCSH1 + s*RCS1 + c*RC1;
                    n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (int r=0; r<R; r++)
                    {
                        Y[n] = powf(X[n1],P[n2]);
                        n++; n1 += r1i; n2 += r2i;
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
        for (int r=0; r<R; r++)
        {
            for (int c=0; c<C; c++)
            {
                for (int s=0; s<S; s++)
                {
                    n1 = r*HSCR1 + c*HSC1 + s*HS1;
                    n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (int h=0; h<H; h++)
                    {
                        Y[n] = powf(X[n1],P[n2]);
                        n++; n1 += h1i; n2 += h2i;
                    }
                }
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int pow_d (double *Y, const double *X, const double *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    if (R1<0) { fprintf(stderr,"error in pow_d: R1 (num rows X) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in pow_d: C1 (num cols X) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in pow_d: S1 (num slices X) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in pow_d: H1 (num hyperslices X) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in pow_d: R2 (num rows P) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in pow_d: C2 (num cols P) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in pow_d: S2 (num slices P) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in pow_d: H2 (num hyperslices P) must be nonnegative\n"); return 1; }

    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    int n1, n2, n = 0;

    if (N1==1)
    {
        while (n<N) { Y[n] = pow(X[0],P[n]); n++; }
    }    
    else if (N2==1)
    {
        while (n<N) { Y[n] = pow(X[n],P[0]); n++; }
    }
    else if (N==N1 && N==N2)
    {
        while (n<N) { Y[n] = pow(X[n],P[n]); n++; }
    }
    else if (iscolmajor)
    {
        const int RCSH1 = R1*C1*S1*(H1>1), RCSH2 = R2*C2*S2*(H2>1);
        const int RCS1 = R1*C1*(S1>1), RCS2 = R2*C2*(S2>1);
        const int RC1 = R1*(C1>1), RC2 = R2*(C2>1);
        const int r1i = (R1>1), r2i = (R2>1);
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                for (int c=0; c<C; c++)
                {
                    n1 = h*RCSH1 + s*RCS1 + c*RC1;
                    n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (int r=0; r<R; r++)
                    {
                        Y[n] = pow(X[n1],P[n2]);
                        n++; n1 += r1i; n2 += r2i;
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
        for (int r=0; r<R; r++)
        {
            for (int c=0; c<C; c++)
            {
                for (int s=0; s<S; s++)
                {
                    n1 = r*HSCR1 + c*HSC1 + s*HS1;
                    n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (int h=0; h<H; h++)
                    {
                        Y[n] = pow(X[n1],P[n2]);
                        n++; n1 += h1i; n2 += h2i;
                    }
                }
            }
        }
    }

    return 0;
}


int pow_c (float *Y, const float *X, const float *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    if (R1<0) { fprintf(stderr,"error in pow_c: R1 (num rows X) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in pow_c: C1 (num cols X) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in pow_c: S1 (num slices X) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in pow_c: H1 (num hyperslices X) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in pow_c: R2 (num rows P) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in pow_c: C2 (num cols P) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in pow_c: S2 (num slices P) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in pow_c: H2 (num hyperslices P) must be nonnegative\n"); return 1; }

    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    int n1, n2, n = 0;
    _Complex float y;

    if (N1==1)
    {
        const _Complex float x = X[0] + 1.0if*X[1];
        while (n<N) { y = cpowf(x,P[n]+1.0if*P[n+1]); memcpy(&Y[n],(float *)&y,2*sizeof(float)); n+=2; }
    }    
    else if (N2==1)
    {
        const _Complex float p = P[0] + 1.0if*P[1];
        while (n<N) { y = cpowf(X[n]+1.0if*X[n+1],p); memcpy(&Y[n],(float *)&y,2*sizeof(float)); n+=2; }
    }
    else if (N==N1 && N==N2)
    {
        while (n<N) { y = cpowf(X[n]+1.0if*X[n+1],P[n]+1.0if*P[n+1]); memcpy(&Y[n],(float *)&y,2*sizeof(float)); n+=2; }
    }
    else if (iscolmajor)
    {
        const int RCSH1 = 2*R1*C1*S1*(H1>1), RCSH2 = 2*R2*C2*S2*(H2>1);
        const int RCS1 = 2*R1*C1*(S1>1), RCS2 = 2*R2*C2*(S2>1);
        const int RC1 = 2*R1*(C1>1), RC2 = 2*R2*(C2>1);
        const int r1i = 2*(R1>1), r2i = 2*(R2>1);
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                for (int c=0; c<C; c++)
                {
                    n1 = h*RCSH1 + s*RCS1 + c*RC1;
                    n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (int r=0; r<R; r++)
                    {
                        y = cpowf(X[n1]+1.0if*X[n1+1],P[n2]+1.0if*P[n2+1]);
                        memcpy(&Y[n],(float *)&y,2*sizeof(float));
                        n+=2; n1 += r1i; n2 += r2i;
                    }
                }
            }
        }
    }
    else
    {
        const int HSCR1 = 2*H1*S1*C1*(R1>1), HSCR2 = 2*H2*S2*C2*(R2>1);
        const int HSC1 = 2*H1*S1*(C1>1), HSC2 = 2*H2*S2*(C2>1);
        const int HS1 = 2*H1*(S1>1), HS2 = 2*H2*(S2>1);
        const int h1i = 2*(H1>1), h2i = 2*(H2>1);
        for (int r=0; r<R; r++)
        {
            for (int c=0; c<C; c++)
            {
                for (int s=0; s<S; s++)
                {
                    n1 = r*HSCR1 + c*HSC1 + s*HS1;
                    n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (int h=0; h<H; h++)
                    {
                        y = cpowf(X[n1]+1.0if*X[n1+1],P[n2]+1.0if*P[n2+1]);
                        memcpy(&Y[n],(float *)&y,2*sizeof(float));
                        n+=2; n1 += h1i; n2 += h2i;
                    }
                }
            }
        }
    }

    return 0;
}


int pow_z (double *Y, const double *X, const double *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    if (R1<0) { fprintf(stderr,"error in pow_z: R1 (num rows X) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in pow_z: C1 (num cols X) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in pow_z: S1 (num slices X) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in pow_z: H1 (num hyperslices X) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in pow_z: R2 (num rows P) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in pow_z: C2 (num cols P) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in pow_z: S2 (num slices P) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in pow_z: H2 (num hyperslices P) must be nonnegative\n"); return 1; }

    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    int n1, n2, n = 0;
    _Complex double y;

    if (N1==1)
    {
        const _Complex double x = X[0] + 1.0i*X[1];
        while (n<N) { y = cpow(x,P[n]+1.0i*P[n+1]); memcpy(&Y[n],(double *)&y,2*sizeof(double)); n+=2; }
    }    
    else if (N2==1)
    {
        const _Complex double p = P[0] + 1.0i*P[1];
        while (n<N) { y = cpow(X[n]+1.0i*X[n+1],p); memcpy(&Y[n],(double *)&y,2*sizeof(double)); n+=2; }
    }
    else if (N==N1 && N==N2)
    {
        while (n<N) { y = cpow(X[n]+1.0i*X[n+1],P[n]+1.0i*P[n+1]); memcpy(&Y[n],(double *)&y,2*sizeof(double)); n+=2; }
    }
    else if (iscolmajor)
    {
        const int RCSH1 = 2*R1*C1*S1*(H1>1), RCSH2 = 2*R2*C2*S2*(H2>1);
        const int RCS1 = 2*R1*C1*(S1>1), RCS2 = 2*R2*C2*(S2>1);
        const int RC1 = 2*R1*(C1>1), RC2 = 2*R2*(C2>1);
        const int r1i = 2*(R1>1), r2i = 2*(R2>1);
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                for (int c=0; c<C; c++)
                {
                    n1 = h*RCSH1 + s*RCS1 + c*RC1;
                    n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (int r=0; r<R; r++)
                    {
                        y = cpow(X[n1]+1.0i*X[n1+1],P[n2]+1.0i*P[n2+1]);
                        memcpy(&Y[n],(double *)&y,2*sizeof(double));
                        n+=2; n1 += r1i; n2 += r2i;
                    }
                }
            }
        }
    }
    else
    {
        const int HSCR1 = 2*H1*S1*C1*(R1>1), HSCR2 = 2*H2*S2*C2*(R2>1);
        const int HSC1 = 2*H1*S1*(C1>1), HSC2 = 2*H2*S2*(C2>1);
        const int HS1 = 2*H1*(S1>1), HS2 = 2*H2*(S2>1);
        const int h1i = 2*(H1>1), h2i = 2*(H2>1);
        for (int r=0; r<R; r++)
        {
            for (int c=0; c<C; c++)
            {
                for (int s=0; s<S; s++)
                {
                    n1 = r*HSCR1 + c*HSC1 + s*HS1;
                    n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (int h=0; h<H; h++)
                    {
                        y = cpow(X[n1]+1.0i*X[n1+1],P[n2]+1.0i*P[n2+1]);
                        memcpy(&Y[n],(double *)&y,2*sizeof(double));
                        n+=2; n1 += h1i; n2 += h2i;
                    }
                }
            }
        }
    }

    return 0;
}


int pow_inplace_s (float *X, const float *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    if (R1<0) { fprintf(stderr,"error in pow_inplace_s: R1 (num rows X) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in pow_inplace_s: C1 (num cols X) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in pow_inplace_s: S1 (num slices X) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in pow_inplace_s: H1 (num hyperslices X) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in pow_inplace_s: R2 (num rows P) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in pow_inplace_s: C2 (num cols P) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in pow_inplace_s: S2 (num slices P) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in pow_inplace_s: H2 (num hyperslices P) must be nonnegative\n"); return 1; }

    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    if (N1!=N) { fprintf(stderr,"error in pow_inplace_s: first input (X) cannot be broadcast for inplace version\n"); return 1; }
    int n2, n = 0;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    if (N2==1)
    {
        while (n<N) { X[n] = powf(X[n],P[0]); n++; }
    }
    else if (N==N2)
    {
        while (n<N) { X[n] = powf(X[n],P[n]); n++; }
    }
    else if (iscolmajor)
    {
        const int RCSH2 = R2*C2*S2*(H2>1), RCS2 = R2*C2*(S2>1), RC2 = R2*(C2>1), r2i = (R2>1);
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                for (int c=0; c<C; c++)
                {
                    n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (int r=0; r<R; r++)
                    {
                        X[n] = powf(X[n],P[n2]);
                        n++; n2 += r2i;
                    }
                }
            }
        }
    }
    else
    {
        const int HSCR2 = H2*S2*C2*(R2>1), HSC2 = H2*S2*(C2>1), HS2 = H2*(S2>1), h2i = (H2>1);
        for (int r=0; r<R; r++)
        {
            for (int c=0; c<C; c++)
            {
                for (int s=0; s<S; s++)
                {
                    n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (int h=0; h<H; h++)
                    {
                        X[n] = powf(X[n],P[n2]);
                        n++; n2 += h2i;
                    }
                }
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int pow_inplace_d (double *X, const double *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    if (R1<0) { fprintf(stderr,"error in pow_inplace_d: R1 (num rows X) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in pow_inplace_d: C1 (num cols X) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in pow_inplace_d: S1 (num slices X) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in pow_inplace_d: H1 (num hyperslices X) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in pow_inplace_d: R2 (num rows P) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in pow_inplace_d: C2 (num cols P) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in pow_inplace_d: S2 (num slices P) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in pow_inplace_d: H2 (num hyperslices P) must be nonnegative\n"); return 1; }

    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    if (N1!=N) { fprintf(stderr,"error in pow_inplace_d: first input (X) cannot be broadcast for inplace version\n"); return 1; }
    int n2, n = 0;

    if (N2==1)
    {
        while (n<N) { X[n] = pow(X[n],P[0]); n++; }
    }
    else if (N==N2)
    {
        while (n<N) { X[n] = pow(X[n],P[n]); n++; }
    }
    else if (iscolmajor)
    {
        const int RCSH2 = R2*C2*S2*(H2>1), RCS2 = R2*C2*(S2>1), RC2 = R2*(C2>1), r2i = (R2>1);
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                for (int c=0; c<C; c++)
                {
                    n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (int r=0; r<R; r++)
                    {
                        X[n] = pow(X[n],P[n2]);
                        n++; n2 += r2i;
                    }
                }
            }
        }
    }
    else
    {
        const int HSCR2 = H2*S2*C2*(R2>1), HSC2 = H2*S2*(C2>1), HS2 = H2*(S2>1), h2i = (H2>1);
        for (int r=0; r<R; r++)
        {
            for (int c=0; c<C; c++)
            {
                for (int s=0; s<S; s++)
                {
                    n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (int h=0; h<H; h++)
                    {
                        X[n] = pow(X[n],P[n2]);
                        n++; n2 += h2i;
                    }
                }
            }
        }
    }

    return 0;
}


int pow_inplace_c (float *X, const float *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    if (R1<0) { fprintf(stderr,"error in pow_inplace_c: R1 (num rows X) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in pow_inplace_c: C1 (num cols X) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in pow_inplace_c: S1 (num slices X) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in pow_inplace_c: H1 (num hyperslices X) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in pow_inplace_c: R2 (num rows P) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in pow_inplace_c: C2 (num cols P) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in pow_inplace_c: S2 (num slices P) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in pow_inplace_c: H2 (num hyperslices P) must be nonnegative\n"); return 1; }

    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    if (N1!=N) { fprintf(stderr,"error in pow_inplace_c: first input (X) cannot be broadcast for inplace version\n"); return 1; }
    int n2, n = 0;
    _Complex float x;

    if (N2==1)
    {
        const _Complex float p = P[0] + 1.0if*P[1];
        while (n<N) { x = cpowf(X[n]+1.0if*X[n+1],p); memcpy(&X[n],(float *)&x,2*sizeof(float)); n+=2; }
    }
    else if (N==N2)
    {
        while (n<N) { x = cpowf(X[n]+1.0if*X[n+1],P[n]+1.0if*P[n+1]); memcpy(&X[n],(float *)&x,2*sizeof(float)); n+=2; }
    }
    else if (iscolmajor)
    {
        const int RCSH2 = 2*R2*C2*S2*(H2>1), RCS2 = 2*R2*C2*(S2>1), RC2 = 2*R2*(C2>1), r2i = 2*(R2>1);
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                for (int c=0; c<C; c++)
                {
                    n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (int r=0; r<R; r++)
                    {
                        x = cpowf(X[n]+1.0if*X[n+1],P[n2]+1.0if*P[n2+1]);
                        memcpy(&X[n],(float *)&x,2*sizeof(float));
                        n+=2; n2 += r2i;
                    }
                }
            }
        }
    }
    else
    {
        const int HSCR2 = 2*H2*S2*C2*(R2>1), HSC2 = 2*H2*S2*(C2>1), HS2 = 2*H2*(S2>1), h2i = 2*(H2>1);
        for (int r=0; r<R; r++)
        {
            for (int c=0; c<C; c++)
            {
                for (int s=0; s<S; s++)
                {
                    n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (int h=0; h<H; h++)
                    {
                        x = cpowf(X[n]+1.0if*X[n+1],P[n2]+1.0if*P[n2+1]);
                        memcpy(&X[n],(float *)&x,2*sizeof(float));
                        n+=2; n2 += h2i;
                    }
                }
            }
        }
    }

    return 0;
}


int pow_inplace_z (double *X, const double *P, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    if (R1<0) { fprintf(stderr,"error in pow_inplace_z: R1 (num rows X) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in pow_inplace_z: C1 (num cols X) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in pow_inplace_z: S1 (num slices X) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in pow_inplace_z: H1 (num hyperslices X) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in pow_inplace_z: R2 (num rows P) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in pow_inplace_z: C2 (num cols P) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in pow_inplace_z: S2 (num slices P) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in pow_inplace_z: H2 (num hyperslices P) must be nonnegative\n"); return 1; }
    
    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    if (N1!=N) { fprintf(stderr,"error in pow_inplace_z: first input (X) cannot be broadcast for inplace version\n"); return 1; }
    int n2, n = 0;
    _Complex double x;

    if (N2==1)
    {
        const _Complex double p = P[0] + 1.0i*P[1];
        while (n<N) { x = cpow(X[n]+1.0i*X[n+1],p); memcpy(&X[n],(double *)&x,2*sizeof(double)); n+=2; }
    }
    else if (N==N2)
    {
        while (n<N) { x = cpow(X[n]+1.0i*X[n+1],P[n]+1.0i*P[n+1]); memcpy(&X[n],(double *)&x,2*sizeof(double)); n+=2; }
    }
    else if (iscolmajor)
    {
        const int RCSH2 = 2*R2*C2*S2*(H2>1), RCS2 = 2*R2*C2*(S2>1), RC2 = 2*R2*(C2>1), r2i = 2*(R2>1);
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                for (int c=0; c<C; c++)
                {
                    n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (int r=0; r<R; r++)
                    {
                        x = cpow(X[n]+1.0i*X[n+1],P[n2]+1.0i*P[n2+1]);
                        memcpy(&X[n],(double *)&x,2*sizeof(double));
                        n+=2; n2 += r2i;
                    }
                }
            }
        }
    }
    else
    {
        const int HSCR2 = 2*H2*S2*C2*(R2>1), HSC2 = 2*H2*S2*(C2>1), HS2 = 2*H2*(S2>1), h2i = 2*(H2>1);
        for (int r=0; r<R; r++)
        {
            for (int c=0; c<C; c++)
            {
                for (int s=0; s<S; s++)
                {
                    n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (int h=0; h<H; h++)
                    {
                        x = cpow(X[n]+1.0i*X[n+1],P[n2]+1.0i*P[n2+1]);
                        memcpy(&X[n],(double *)&x,2*sizeof(double));
                        n+=2; n2 += h2i;
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
