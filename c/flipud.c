//Flips X along cols (reverses each col).
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int flipud_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor);
int flipud_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor);
int flipud_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor);
int flipud_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor);

int flipud_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor);
int flipud_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor);
int flipud_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor);
int flipud_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor);


int flipud_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1 = 0, n2;
    //struct timespec tic, toc;

    //Checks
    if (R<0) { fprintf(stderr,"error in flipud_s: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flipud_s: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flipud_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flipud_s: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&tic);
    if (iscolmajor)
    {
        const int RC = R*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                for (c=0; c<C; c++)
                {
                    n2 = h*RCS + s*RC + c*R + R - 1;
                    for (r=0; r<R; r++)
                    {
                        Y[n1] = X[n2];
                        n1++; n2--;
                    }
                }
            }
        }
    }
    else
    {
        const int CSH = C*S*H;
        n2 = (R-1)*CSH;
        for (r=0; r<R; r++)
        {
            cblas_scopy(CSH,&X[n2],1,&Y[n1],1);
            n1 += CSH;  n2 -= CSH;
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int flipud_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1 = 0, n2;

    //Checks
    if (R<0) { fprintf(stderr,"error in flipud_d: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flipud_d: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flipud_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flipud_d: H (num hyperslices X) must be nonnegative\n"); return 1; }

    if (iscolmajor)
    {
        const int RC = R*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                for (c=0; c<C; c++)
                {
                    n2 = h*RCS + s*RC + c*R + R - 1;
                    for (r=0; r<R; r++)
                    {
                        Y[n1] = X[n2];
                        n1++; n2--;
                    }
                }
            }
        }
    }
    else
    {
        const int CSH = C*S*H;
        n2 = (R-1)*CSH;
        for (r=0; r<R; r++)
        {
            cblas_dcopy(CSH,&X[n2],1,&Y[n1],1);
            n1 += CSH;  n2 -= CSH;
        }
    }
    
    return 0;
}


int flipud_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1 = 0, n2;

    //Checks
    if (R<0) { fprintf(stderr,"error in flipud_c: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flipud_c: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flipud_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flipud_c: H (num hyperslices X) must be nonnegative\n"); return 1; }

    if (iscolmajor)
    {
        const int R2 = 2*R, RC = R2*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                for (c=0; c<C; c++)
                {
                    n2 = h*RCS + s*RC + c*R2 + R2 - 1;
                    for (r=0; r<R; r++)
                    {
                        Y[n1] = X[n2];
                        n1++; n2--;
                        Y[n1] = X[n2];
                        n1++; n2--;
                    }
                }
            }
        }
    }
    else
    {
        const int CSH = 2*C*S*H;
        n2 = (R-1)*CSH;
        for (r=0; r<R; r++)
        {
            cblas_scopy(CSH,&X[n2],1,&Y[n1],1);
            n1 += CSH;  n2 -= CSH;
        }
    }
    
    return 0;
}


int flipud_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1 = 0, n2;

    //Checks
    if (R<0) { fprintf(stderr,"error in flipud_z: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flipud_z: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flipud_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flipud_z: H (num hyperslices X) must be nonnegative\n"); return 1; }

    if (iscolmajor)
    {
        const int R2 = 2*R, RC = R2*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                for (c=0; c<C; c++)
                {
                    n2 = h*RCS + s*RC + c*R2 + R2 - 1;
                    for (r=0; r<R; r++)
                    {
                        Y[n1] = X[n2];
                        n1++; n2--;
                        Y[n1] = X[n2];
                        n1++; n2--;
                    }
                }
            }
        }
    }
    else
    {
        const int CSH = 2*C*S*H;
        n2 = (R-1)*CSH;
        for (r=0; r<R; r++)
        {
            cblas_dcopy(CSH,&X[n2],1,&Y[n1],1);
            n1 += CSH;  n2 -= CSH;
        }
    }
    
    return 0;
}


int flipud_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1, n2;
    float x;
    //struct timespec tic, toc;

    //Checks
    if (R<0) { fprintf(stderr,"error in flipud_s: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flipud_s: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flipud_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flipud_s: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&tic);
    if (iscolmajor)
    {
        const int RC = R*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                for (c=0; c<C; c++)
                {
                    n1 = h*RCS + s*RC + c*R;
                    n2 = n1 + R - 1;
                    for (r=0; r<R/2; r++)
                    {
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2--;
                    }
                }
            }
        }
    }
    else
    {
        const int CSH = C*S*H;
        for (r=0; r<R/2; r++)
        {
            n1 = r*CSH;
            n2 = (R-r-1)*CSH;
            for (c=0; c<C; c++)
            {
                for (s=0; s<S; s++)
                {
                    for (h=0; h<H; h++)
                    {
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2++;
                    }
                }
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int flipud_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1, n2;
    double x;

    //Checks
    if (R<0) { fprintf(stderr,"error in flipud_d: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flipud_d: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flipud_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flipud_d: H (num hyperslices X) must be nonnegative\n"); return 1; }

    if (iscolmajor)
    {
        const int RC = R*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                for (c=0; c<C; c++)
                {
                    n1 = h*RCS + s*RC + c*R;
                    n2 = n1 + R - 1;
                    for (r=0; r<R/2; r++)
                    {
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2--;
                    }
                }
            }
        }
    }
    else
    {
        const int CSH = C*S*H;
        for (r=0; r<R/2; r++)
        {
            n1 = r*CSH;
            n2 = (R-r-1)*CSH;
            for (c=0; c<C; c++)
            {
                for (s=0; s<S; s++)
                {
                    for (h=0; h<H; h++)
                    {
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2++;
                    }
                }
            }
        }
    }
    
    return 0;
}


int flipud_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1, n2;
    float x;

    //Checks
    if (R<0) { fprintf(stderr,"error in flipud_c: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flipud_c: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flipud_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flipud_c: H (num hyperslices X) must be nonnegative\n"); return 1; }

    if (iscolmajor)
    {
        const int R2 = 2*R, RC = R2*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                for (c=0; c<C; c++)
                {
                    n1 = h*RCS + s*RC + c*R2;
                    n2 = n1 + R2 - 2;
                    for (r=0; r<R/2; r++)
                    {
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        x = X[n1+1]; X[n1+1] = X[n2+1]; X[n2+1] = x;
                        n1+=2; n2-=2;
                    }
                }
            }
        }
    }
    else
    {
        const int CSH = 2*C*S*H;
        for (r=0; r<R/2; r++)
        {
            n1 = r*CSH;
            n2 = (R-r-1)*CSH;
            for (c=0; c<C; c++)
            {
                for (s=0; s<S; s++)
                {
                    for (h=0; h<H; h++)
                    {
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2++;
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2++;
                    }
                }
            }
        }
    }

    return 0;
}


int flipud_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1, n2;
    double x;

    //Checks
    if (R<0) { fprintf(stderr,"error in flipud_z: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flipud_z: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flipud_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flipud_z: H (num hyperslices X) must be nonnegative\n"); return 1; }

    if (iscolmajor)
    {
        const int R2 = 2*R, RC = R2*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                for (c=0; c<C; c++)
                {
                    n1 = h*RCS + s*RC + c*R2;
                    n2 = n1 + R2 - 2;
                    for (r=0; r<R/2; r++)
                    {
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        x = X[n1+1]; X[n1+1] = X[n2+1]; X[n2+1] = x;
                        n1+=2; n2-=2;
                    }
                }
            }
        }
    }
    else
    {
        const int CSH = 2*C*S*H;
        for (r=0; r<R/2; r++)
        {
            n1 = r*CSH;
            n2 = (R-r-1)*CSH;
            for (c=0; c<C; c++)
            {
                for (s=0; s<S; s++)
                {
                    for (h=0; h<H; h++)
                    {
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2++;
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2++;
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
