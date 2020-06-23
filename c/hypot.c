//2-input elementwise function.
//Gets hypotenuse for corresponding elements in X11, X12: y = sqrt(x1^2 + x2^2)
//For complex input, output is real-valued: y = sqrt(|x1|^2 + |x2|^2)
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int hypot_s (float *Y, const float *X1, const float *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);
int hypot_d (double *Y, const double *X1, const double *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);
int hypot_c (float *Y, const float *X1, const float *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);
int hypot_z (double *Y, const double *X1, const double *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);

int hypot_inplace_s (float *X1, const float *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);
int hypot_inplace_d (double *X1, const double *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);
int hypot_inplace_c (float *X1, const float *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);
int hypot_inplace_z (double *X1, const double *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor);


int hypot_s (float *Y, const float *X1, const float *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    int r, c, s, h, n1, n2, n = 0;
    //struct timespec tic, toc;

    //Checks
    if (R1<0) { fprintf(stderr,"error in hypot_s: R1 (num rows X1) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in hypot_s: C1 (num cols X1) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in hypot_s: S1 (num slices X1) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in hypot_s: H1 (num hyperslices X1) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in hypot_s: R2 (num rows X2) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in hypot_s: C2 (num cols X2) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in hypot_s: S2 (num slices X2) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in hypot_s: H2 (num hyperslices X2) must be nonnegative\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&tic);
    if (N1==1)
    {
        while (n<N) { Y[n] = hypotf(X1[0],X2[n]); n++; }
    }    
    else if (N2==1)
    {
        while (n<N) { Y[n] = hypotf(X1[n],X2[0]); n++; }
    }
    else if (N==N1 && N==N2)
    {
        while (n<N) { Y[n] = hypotf(X1[n],X2[n]); n++; }
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
                        Y[n] = hypotf(X1[n1],X2[n2]);
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
                        Y[n] = hypotf(X1[n1],X2[n2]);
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


int hypot_d (double *Y, const double *X1, const double *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    int r, c, s, h, n1, n2, n = 0;

    //Checks
    if (R1<0) { fprintf(stderr,"error in hypot_d: R1 (num rows X1) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in hypot_d: C1 (num cols X1) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in hypot_d: S1 (num slices X1) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in hypot_d: H1 (num hyperslices X1) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in hypot_d: R2 (num rows X2) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in hypot_d: C2 (num cols X2) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in hypot_d: S2 (num slices X2) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in hypot_d: H2 (num hyperslices X2) must be nonnegative\n"); return 1; }

    if (N1==1)
    {
        while (n<N) { Y[n] = hypot(X1[0],X2[n]); n++; }
    }    
    else if (N2==1)
    {
        while (n<N) { Y[n] = hypot(X1[n],X2[0]); n++; }
    }
    else if (N==N1 && N==N2)
    {
        while (n<N) { Y[n] = hypot(X1[n],X2[n]); n++; }
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
                        Y[n] = hypot(X1[n1],X2[n2]);
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
                        Y[n] = hypot(X1[n1],X2[n2]);
                        n++; n1 += h1i; n2 += h2i;
                    }
                }
            }
        }
    }

    return 0;
}


int hypot_c (float *Y, const float *X1, const float *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    int r, c, s, h, n1, n2 = 0, n = 0;

    //Checks
    if (R1<0) { fprintf(stderr,"error in hypot_c: R1 (num rows X1) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in hypot_c: C1 (num cols X1) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in hypot_c: S1 (num slices X1) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in hypot_c: H1 (num hyperslices X1) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in hypot_c: R2 (num rows X2) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in hypot_c: C2 (num cols X2) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in hypot_c: S2 (num slices X2) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in hypot_c: H2 (num hyperslices X2) must be nonnegative\n"); return 1; }

    if (N1==1)
    {
        const float x1 = X1[0]*X1[0] + X1[1]*X1[1];
        while (n<N) { Y[n] = sqrtf(x1 + X2[n2]*X2[n2]+X2[n2+1]*X2[n2+1]); n++; n2+=2; }
    }    
    else if (N2==1)
    {
        const float x2 = X2[0]*X2[0] + X2[1]*X2[1];
        while (n<N) { Y[n] = sqrtf(X1[n2]*X1[n2]+X1[n2+1]*X1[n2+1] + x2); n++; n2+=2; }
    }
    else if (N==N1 && N==N2)
    {
        while (n<N) { Y[n] = sqrtf(X1[n2]*X1[n2]+X1[n2+1]*X1[n2+1] + X2[n2]*X2[n2]+X2[n2+1]*X2[n2+1]); n++; n2+=2; }
    }
    else if (iscolmajor)
    {
        const int RCSH1 = 2*R1*C1*S1*(H1>1), RCSH2 = 2*R2*C2*S2*(H2>1);
        const int RCS1 = 2*R1*C1*(S1>1), RCS2 = 2*R2*C2*(S2>1);
        const int RC1 = 2*R1*(C1>1), RC2 = 2*R2*(C2>1);
        const int r1i = 2*(R1>1), r2i = 2*(R2>1);
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
                        Y[n] = sqrtf(X1[n1]*X1[n1]+X1[n1+1]*X1[n1+1] + X2[n2]*X2[n2]+X2[n2+1]*X2[n2+1]);
                        n++; n1 += r1i; n2 += r2i;
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
                        Y[n] = sqrtf(X1[n1]*X1[n1]+X1[n1+1]*X1[n1+1] + X2[n2]*X2[n2]+X2[n2+1]*X2[n2+1]);
                        n++; n1 += h1i; n2 += h2i;
                    }
                }
            }
        }
    }

    return 0;
}


int hypot_z (double *Y, const double *X1, const double *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    int r, c, s, h, n1, n2 = 0, n = 0;

    //Checks
    if (R1<0) { fprintf(stderr,"error in hypot_z: R1 (num rows X1) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in hypot_z: C1 (num cols X1) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in hypot_z: S1 (num slices X1) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in hypot_z: H1 (num hyperslices X1) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in hypot_z: R2 (num rows X2) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in hypot_z: C2 (num cols X2) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in hypot_z: S2 (num slices X2) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in hypot_z: H2 (num hyperslices X2) must be nonnegative\n"); return 1; }

    if (N1==1)
    {
        const double x1 = X1[0]*X1[0] + X1[1]*X1[1];
        while (n<N) { Y[n] = sqrt(x1 + X2[n2]*X2[n2]+X2[n2+1]*X2[n2+1]); n++; n2+=2; }
    }    
    else if (N2==1)
    {
        const double x2 = X2[0]*X2[0] + X2[1]*X2[1];
        while (n<N) { Y[n] = sqrt(X1[n2]*X1[n2]+X1[n2+1]*X1[n2+1] + x2); n++; n2+=2; }
    }
    else if (N==N1 && N==N2)
    {
        while (n<N) { Y[n] = sqrt(X1[n2]*X1[n2]+X1[n2+1]*X1[n2+1] + X2[n2]*X2[n2]+X2[n2+1]*X2[n2+1]); n++; n2+=2; }
    }
    else if (iscolmajor)
    {
        const int RCSH1 = 2*R1*C1*S1*(H1>1), RCSH2 = 2*R2*C2*S2*(H2>1);
        const int RCS1 = 2*R1*C1*(S1>1), RCS2 = 2*R2*C2*(S2>1);
        const int RC1 = 2*R1*(C1>1), RC2 = 2*R2*(C2>1);
        const int r1i = 2*(R1>1), r2i = 2*(R2>1);
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
                        Y[n] = sqrt(X1[n1]*X1[n1]+X1[n1+1]*X1[n1+1] + X2[n2]*X2[n2]+X2[n2+1]*X2[n2+1]);
                        n++; n1 += r1i; n2 += r2i;
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
                        Y[n] = sqrt(X1[n1]*X1[n1]+X1[n1+1]*X1[n1+1] + X2[n2]*X2[n2]+X2[n2+1]*X2[n2+1]);
                        n++; n1 += h1i; n2 += h2i;
                    }
                }
            }
        }
    }

    return 0;
}


int hypot_inplace_s (float *X1, const float *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    int r, c, s, h, n2, n = 0;
    //struct timespec tic, toc;

    //Checks
    if (R1<0) { fprintf(stderr,"error in hypot_inplace_s: R1 (num rows X1) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in hypot_inplace_s: C1 (num cols X1) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in hypot_inplace_s: S1 (num slices X1) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in hypot_inplace_s: H1 (num hyperslices X1) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in hypot_inplace_s: R2 (num rows X2) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in hypot_inplace_s: C2 (num cols X2) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in hypot_inplace_s: S2 (num slices X2) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in hypot_inplace_s: H2 (num hyperslices X2) must be nonnegative\n"); return 1; }
    if (N1!=N) { fprintf(stderr,"error in hypot_inplace_s: first input (X1) cannot be broadcast for inplace version\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&tic);
    if (N2==1)
    {
        while (n<N) { X1[n] = hypotf(X1[n],X2[0]); n++; }
    }
    else if (N==N2)
    {
        while (n<N) { X1[n] = hypotf(X1[n],X2[n]); n++; }
    }
    else if (iscolmajor)
    {
        const int RCSH2 = R2*C2*S2*(H2>1), RCS2 = R2*C2*(S2>1), RC2 = R2*(C2>1), r2i = (R2>1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                for (c=0; c<C; c++)
                {
                    n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (r=0; r<R; r++)
                    {
                        X1[n] = hypotf(X1[n],X2[n2]);
                        n++; n2 += r2i;
                    }
                }
            }
        }
    }
    else
    {
        const int HSCR2 = H2*S2*C2*(R2>1), HSC2 = H2*S2*(C2>1), HS2 = H2*(S2>1), h2i = (H2>1);
        for (r=0; r<R; r++)
        {
            for (c=0; c<C; c++)
            {
                for (s=0; s<S; s++)
                {
                    n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (h=0; h<H; h++)
                    {
                        X1[n] = hypotf(X1[n],X2[n2]);
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


int hypot_inplace_d (double *X1, const double *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    int r, c, s, h, n2, n = 0;

    //Checks
    if (R1<0) { fprintf(stderr,"error in hypot_inplace_d: R1 (num rows X1) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in hypot_inplace_d: C1 (num cols X1) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in hypot_inplace_d: S1 (num slices X1) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in hypot_inplace_d: H1 (num hyperslices X1) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in hypot_inplace_d: R2 (num rows X2) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in hypot_inplace_d: C2 (num cols X2) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in hypot_inplace_d: S2 (num slices X2) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in hypot_inplace_d: H2 (num hyperslices X2) must be nonnegative\n"); return 1; }
    if (N1!=N) { fprintf(stderr,"error in hypot_inplace_d: first input (X1) cannot be broadcast for inplace version\n"); return 1; }

    if (N2==1)
    {
        while (n<N) { X1[n] = hypot(X1[n],X2[0]); n++; }
    }
    else if (N==N2)
    {
        while (n<N) { X1[n] = hypot(X1[n],X2[n]); n++; }
    }
    else if (iscolmajor)
    {
        const int RCSH2 = R2*C2*S2*(H2>1), RCS2 = R2*C2*(S2>1), RC2 = R2*(C2>1), r2i = (R2>1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                for (c=0; c<C; c++)
                {
                    n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (r=0; r<R; r++)
                    {
                        X1[n] = hypot(X1[n],X2[n2]);
                        n++; n2 += r2i;
                    }
                }
            }
        }
    }
    else
    {
        const int HSCR2 = H2*S2*C2*(R2>1), HSC2 = H2*S2*(C2>1), HS2 = H2*(S2>1), h2i = (H2>1);
        for (r=0; r<R; r++)
        {
            for (c=0; c<C; c++)
            {
                for (s=0; s<S; s++)
                {
                    n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (h=0; h<H; h++)
                    {
                        X1[n] = hypot(X1[n],X2[n2]);
                        n++; n2 += h2i;
                    }
                }
            }
        }
    }

    return 0;
}


int hypot_inplace_c (float *X1, const float *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    int r, c, s, h, n1 = 0, n2 = 0, n = 0;

    //Checks
    if (R1<0) { fprintf(stderr,"error in hypot_inplace_c: R1 (num rows X1) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in hypot_inplace_c: C1 (num cols X1) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in hypot_inplace_c: S1 (num slices X1) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in hypot_inplace_c: H1 (num hyperslices X1) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in hypot_inplace_c: R2 (num rows X2) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in hypot_inplace_c: C2 (num cols X2) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in hypot_inplace_c: S2 (num slices X2) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in hypot_inplace_c: H2 (num hyperslices X2) must be nonnegative\n"); return 1; }
    if (N1!=N) { fprintf(stderr,"error in hypot_inplace_c: first input (X1) cannot be broadcast for inplace version\n"); return 1; }

    if (N2==1)
    {
        const float x2 = X2[0]*X2[0] + X2[1]*X2[1];
        while (n<N) { X1[n] = sqrtf(X1[n2]*X1[n2]+X1[n2+1]*X1[n2+1] + x2); n++; n2+=2; }
    }
    else if (N==N1 && N==N2)
    {
        while (n<N) { X1[n] = sqrtf(X1[n2]*X1[n2]+X1[n2+1]*X1[n2+1] + X2[n2]*X2[n2]+X2[n2+1]*X2[n2+1]); n++; n2+=2; }
    }
    else if (iscolmajor)
    {
        const int RCSH2 = 2*R2*C2*S2*(H2>1), RCS2 = 2*R2*C2*(S2>1), RC2 = 2*R2*(C2>1), r2i = 2*(R2>1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                for (c=0; c<C; c++)
                {
                    n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (r=0; r<R; r++)
                    {
                        X1[n] = sqrtf(X1[n1]*X1[n1]+X1[n1+1]*X1[n1+1] + X2[n2]*X2[n2]+X2[n2+1]*X2[n2+1]);
                        n++; n1 += 2; n2 += r2i;
                    }
                }
            }
        }
    }
    else
    {
        const int HSCR2 = 2*H2*S2*C2*(R2>1), HSC2 = 2*H2*S2*(C2>1), HS2 = 2*H2*(S2>1), h2i = 2*(H2>1);
        for (r=0; r<R; r++)
        {
            for (c=0; c<C; c++)
            {
                for (s=0; s<S; s++)
                {
                    n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (h=0; h<H; h++)
                    {
                        X1[n] = sqrtf(X1[n1]*X1[n1]+X1[n1+1]*X1[n1+1] + X2[n2]*X2[n2]+X2[n2+1]*X2[n2+1]);
                        n++; n1 += 2; n2 += h2i;
                    }
                }
            }
        }
    }

    return 0;
}


int hypot_inplace_z (double *X1, const double *X2, const int R1, const int C1, const int S1, const int H1, const int R2, const int C2, const int S2, const int H2, const char iscolmajor)
{
    const int R = (R1>R2) ? R1 : R2;
    const int C = (C1>C2) ? C1 : C2;
    const int S = (S1>S2) ? S1 : S2;
    const int H = (H1>H2) ? H1 : H2;
    const int N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    int r, c, s, h, n1 = 0, n2 = 0, n = 0;

    //Checks
    if (R1<0) { fprintf(stderr,"error in hypot_inplace_z: R1 (num rows X1) must be nonnegative\n"); return 1; }
    if (C1<0) { fprintf(stderr,"error in hypot_inplace_z: C1 (num cols X1) must be nonnegative\n"); return 1; }
    if (S1<0) { fprintf(stderr,"error in hypot_inplace_z: S1 (num slices X1) must be nonnegative\n"); return 1; }
    if (H1<0) { fprintf(stderr,"error in hypot_inplace_z: H1 (num hyperslices X1) must be nonnegative\n"); return 1; }
    if (R2<0) { fprintf(stderr,"error in hypot_inplace_z: R2 (num rows X2) must be nonnegative\n"); return 1; }
    if (C2<0) { fprintf(stderr,"error in hypot_inplace_z: C2 (num cols X2) must be nonnegative\n"); return 1; }
    if (S2<0) { fprintf(stderr,"error in hypot_inplace_z: S2 (num slices X2) must be nonnegative\n"); return 1; }
    if (H2<0) { fprintf(stderr,"error in hypot_inplace_z: H2 (num hyperslices X2) must be nonnegative\n"); return 1; }
    if (N1!=N) { fprintf(stderr,"error in hypot_inplace_z: first input (X1) cannot be broadcast for inplace version\n"); return 1; }

    if (N2==1)
    {
        const double x2 = X2[0]*X2[0] + X2[1]*X2[1];
        while (n<N) { X1[n] = sqrt(X1[n2]*X1[n2]+X1[n2+1]*X1[n2+1] + x2); n++; n2+=2; }
    }
    else if (N==N1 && N==N2)
    {
        while (n<N) { X1[n] = sqrt(X1[n2]*X1[n2]+X1[n2+1]*X1[n2+1] + X2[n2]*X2[n2]+X2[n2+1]*X2[n2+1]); n++; n2+=2; }
    }
    else if (iscolmajor)
    {
        const int RCSH2 = 2*R2*C2*S2*(H2>1), RCS2 = 2*R2*C2*(S2>1), RC2 = 2*R2*(C2>1), r2i = 2*(R2>1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                for (c=0; c<C; c++)
                {
                    n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (r=0; r<R; r++)
                    {
                        X1[n] = sqrt(X1[n1]*X1[n1]+X1[n1+1]*X1[n1+1] + X2[n2]*X2[n2]+X2[n2+1]*X2[n2+1]);
                        n++; n1 += 2; n2 += r2i;
                    }
                }
            }
        }
    }
    else
    {
        const int HSCR2 = 2*H2*S2*C2*(R2>1), HSC2 = 2*H2*S2*(C2>1), HS2 = 2*H2*(S2>1), h2i = 2*(H2>1);
        for (r=0; r<R; r++)
        {
            for (c=0; c<C; c++)
            {
                for (s=0; s<S; s++)
                {
                    n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (h=0; h<H; h++)
                    {
                        X1[n] = sqrt(X1[n1]*X1[n1]+X1[n1+1]*X1[n1+1] + X2[n2]*X2[n2]+X2[n2+1]*X2[n2+1]);
                        n++; n1 += 2; n2 += h2i;
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
