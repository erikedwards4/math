//Flips X along rows (reverses each row).
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int fliplr_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor);
int fliplr_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor);
int fliplr_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor);
int fliplr_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor);

int fliplr_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor);
int fliplr_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor);
int fliplr_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor);
int fliplr_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor);


int fliplr_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1 = 0, n2;
    //struct timespec tic, toc;

    //Checks
    if (R<0) { fprintf(stderr,"error in fliplr_s: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in fliplr_s: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in fliplr_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in fliplr_s: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&tic);
    if (iscolmajor)
    {
        const int RC = R*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                n2 = h*RCS + s*RC + (C-1)*R;
                for (c=0; c<C; c++)
                {
                    cblas_scopy(R,&X[n2],1,&Y[n1],1);
                    n1 += R; n2 -= R;
                }
            }
        }
    }
    else
    {
        const int SH = S*H, CSH = C*SH;
        for (r=0; r<R; r++)
        {
            n2 = r*CSH + (C-1)*SH;
            for (c=0; c<C; c++)
            {
                cblas_scopy(SH,&X[n2],1,&Y[n1],1);
                n1 += SH; n2 -= SH;
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int fliplr_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1 = 0, n2;

    //Checks
    if (R<0) { fprintf(stderr,"error in fliplr_d: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in fliplr_d: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in fliplr_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in fliplr_d: H (num hyperslices X) must be nonnegative\n"); return 1; }

    if (iscolmajor)
    {
        const int RC = R*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                n2 = h*RCS + s*RC + (C-1)*R;
                for (c=0; c<C; c++)
                {
                    cblas_dcopy(R,&X[n2],1,&Y[n1],1);
                    n1 += R; n2 -= R;
                }
            }
        }
    }
    else
    {
        const int SH = S*H, CSH = C*SH;
        for (r=0; r<R; r++)
        {
            n2 = r*CSH + (C-1)*SH;
            for (c=0; c<C; c++)
            {
                cblas_dcopy(SH,&X[n2],1,&Y[n1],1);
                n1 += SH; n2 -= SH;
            }
        }
    }
    
    return 0;
}


int fliplr_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1 = 0, n2;

    //Checks
    if (R<0) { fprintf(stderr,"error in fliplr_c: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in fliplr_c: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in fliplr_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in fliplr_c: H (num hyperslices X) must be nonnegative\n"); return 1; }

    if (iscolmajor)
    {
        const int R2 = 2*R, RC = R2*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                n2 = h*RCS + s*RC + (C-1)*R2;
                for (c=0; c<C; c++)
                {
                    cblas_scopy(R2,&X[n2],1,&Y[n1],1);
                    n1 += R2; n2 -= R2;
                }
            }
        }
    }
    else
    {
        const int SH = 2*S*H, CSH = C*SH;
        for (r=0; r<R; r++)
        {
            n2 = r*CSH + (C-1)*SH;
            for (c=0; c<C; c++)
            {
                cblas_scopy(SH,&X[n2],1,&Y[n1],1);
                n1 += SH; n2 -= SH;
            }
        }
    }
    
    return 0;
}


int fliplr_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1 = 0, n2;

    //Checks
    if (R<0) { fprintf(stderr,"error in fliplr_z: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in fliplr_z: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in fliplr_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in fliplr_z: H (num hyperslices X) must be nonnegative\n"); return 1; }

    if (iscolmajor)
    {
        const int R2 = 2*R, RC = R2*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                n2 = h*RCS + s*RC + (C-1)*R2;
                for (c=0; c<C; c++)
                {
                    cblas_dcopy(R2,&X[n2],1,&Y[n1],1);
                    n1 += R2; n2 -= R2;
                }
            }
        }
    }
    else
    {
        const int SH = 2*S*H, CSH = C*SH;
        for (r=0; r<R; r++)
        {
            n2 = r*CSH + (C-1)*SH;
            for (c=0; c<C; c++)
            {
                cblas_dcopy(SH,&X[n2],1,&Y[n1],1);
                n1 += SH; n2 -= SH;
            }
        }
    }
    
    return 0;
}


int fliplr_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1, n2;
    float x;
    //struct timespec tic, toc;

    //Checks
    if (R<0) { fprintf(stderr,"error in fliplr_inplace_s: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in fliplr_inplace_s: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in fliplr_inplace_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in fliplr_inplace_s: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&tic);
    if (iscolmajor)
    {
        const int RC = R*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                n1 = h*RCS + s*RC; n2 = n1 + (C-1)*R;
                for (c=0; c<C/2; c++)
                {
                    for (r=0; r<R; r++)
                    {
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2++;
                    }
                    n2 -= 2*R;
                }
            }
        }
    }
    else
    {
        const int SH = S*H, CSH = C*SH;
        for (r=0; r<R; r++)
        {
            n1 = r*CSH; n2 = n1 + (C-1)*SH;
            for (c=0; c<C/2; c++)
            {
                for (s=0; s<S; s++)
                {
                    for (h=0; h<H; h++)
                    {
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2++;
                    }
                }
                n2 -= 2*SH;
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int fliplr_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1, n2;
    double x;

    //Checks
    if (R<0) { fprintf(stderr,"error in fliplr_inplace_d: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in fliplr_inplace_d: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in fliplr_inplace_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in fliplr_inplace_d: H (num hyperslices X) must be nonnegative\n"); return 1; }

    if (iscolmajor)
    {
        const int RC = R*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                n1 = h*RCS + s*RC; n2 = n1 + (C-1)*R;
                for (c=0; c<C/2; c++)
                {
                    for (r=0; r<R; r++)
                    {
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2++;
                    }
                    n2 -= 2*R;
                }
            }
        }
    }
    else
    {
        const int SH = S*H, CSH = C*SH;
        for (r=0; r<R; r++)
        {
            n1 = r*CSH; n2 = n1 + (C-1)*SH;
            for (c=0; c<C/2; c++)
            {
                for (s=0; s<S; s++)
                {
                    for (h=0; h<H; h++)
                    {
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2++;
                    }
                }
                n2 -= 2*SH;
            }
        }
    }
    
    return 0;
}


int fliplr_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1, n2;
    float x;

    //Checks
    if (R<0) { fprintf(stderr,"error in fliplr_inplace_c: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in fliplr_inplace_c: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in fliplr_inplace_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in fliplr_inplace_c: H (num hyperslices X) must be nonnegative\n"); return 1; }

    if (iscolmajor)
    {
        const int R2 = 2*R, RC = R2*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                n1 = h*RCS + s*RC; n2 = n1 + (C-1)*R2;
                for (c=0; c<C/2; c++)
                {
                    for (r=0; r<R; r++)
                    {
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2++;
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2++;
                    }
                    n2 -= 2*R2;
                }
            }
        }
    }
    else
    {
        const int SH = 2*S*H, CSH = C*SH;
        for (r=0; r<R; r++)
        {
            n1 = r*CSH; n2 = n1 + (C-1)*SH;
            for (c=0; c<C/2; c++)
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
                n2 -= 2*SH;
            }
        }
    }

    return 0;
}


int fliplr_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor)
{
    int r, c, s, h, n1, n2;
    double x;

    //Checks
    if (R<0) { fprintf(stderr,"error in fliplr_inplace_z: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in fliplr_inplace_z: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in fliplr_inplace_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in fliplr_inplace_z: H (num hyperslices X) must be nonnegative\n"); return 1; }

    if (iscolmajor)
    {
        const int R2 = 2*R, RC = R2*C, RCS = RC*S;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                n1 = h*RCS + s*RC; n2 = n1 + (C-1)*R2;
                for (c=0; c<C/2; c++)
                {
                    for (r=0; r<R; r++)
                    {
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2++;
                        x = X[n1]; X[n1] = X[n2]; X[n2] = x;
                        n1++; n2++;
                    }
                    n2 -= 2*R2;
                }
            }
        }
    }
    else
    {
        const int SH = 2*S*H, CSH = C*SH;
        for (r=0; r<R; r++)
        {
            n1 = r*CSH; n2 = n1 + (C-1)*SH;
            for (c=0; c<C/2; c++)
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
                n2 -= 2*SH;
            }
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
