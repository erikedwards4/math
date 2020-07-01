//Zeros all elements above the kth diagonal of matrix X.
//This has in-place and not-in-place versions.

//LAPACKE_?lacpy was surprisingly a bit slower even for square, k==0 case.

#include <stdio.h>
#include <cblas.h>
#include <lapacke.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int tril_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);
int tril_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);
int tril_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);
int tril_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);

int tril_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);
int tril_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);
int tril_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);
int tril_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);


int tril_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in tril_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in tril_s: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in tril_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in tril_s: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in tril_s: k must be in [1-R C-1]\n"); return 1; }

    const float z = 0.0f;
    int n = 0;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    // if (k==0 && R==C && S*H==1)
    // {
    //     const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
    //     const int lda = (iscolmajor) ? R : C;
    //     if (LAPACKE_slacpy (Ord,'L',R,C,X,lda,Y,lda))
    //     { fprintf(stderr,"error in tril_s: problem with LAPACKE function\n"); return 1; }
    // }
    if (iscolmajor)
    {
        const int C0 = (k<C-R) ? C-R-k : 0; //number of all-0 cols
        const int CX = (k>-1) ? k+1 : 0;    //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (CX>0) { cblas_scopy(R*CX,&X[n],1,&Y[n],1); n += R*CX; }
                for (int c=CX; c<C-C0; c++)
                {
                    cblas_scopy(c-k,&z,0,&Y[n],1); n += c-k;
                    cblas_scopy(R-c+k,&X[n],1,&Y[n],1); n += R-c+k;
                }
                if (C0>0) { cblas_scopy(R*C0,&z,0,&Y[n],1); n += R*C0; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (k<0) ? -k : 0;          //number of all-0 rows
        const int RX = (k>C-R-1) ? R-C+k+1 : 0; //number of all-X rows
        if (R0>0) { cblas_scopy(R0*C*SH,&z,0,&Y[n],1); n += R0*C*SH; }
        for (int r=R0; r<R-RX; r++)
        {
            cblas_scopy((r+k+1)*SH,&X[n],1,&Y[n],1); n += (r+k+1)*SH;
            cblas_scopy((C-r-k-1)*SH,&z,0,&Y[n],1); n += (C-r-k-1)*SH;
        }
        if (RX>0) { cblas_scopy(RX*C*SH,&X[n],1,&Y[n],1); n += RX*C*SH; }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int tril_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in tril_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in tril_d: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in tril_d: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in tril_d: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in tril_d: k must be in [1-R C-1]\n"); return 1; }

    const double z = 0.0;
    int n = 0;

    if (iscolmajor)
    {
        const int C0 = (C-R-k>0) ? C-R-k : 0; //number of all-0 cols
        const int CX = (k+1>0) ? k+1 : 0;     //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (CX>0) { cblas_dcopy(R*CX,&X[n],1,&Y[n],1); n += R*CX; }
                for (int c=CX; c<C-C0; c++)
                {
                    cblas_dcopy(c-k,&z,0,&Y[n],1); n += c-k;
                    cblas_dcopy(R-c+k,&X[n],1,&Y[n],1); n += R-c+k;
                }
                if (C0>0) { cblas_dcopy(R*C0,&z,0,&Y[n],1); n += R*C0; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (-k>0) ? -k : 0;           //number of all-0 rows
        const int RX = (R-C+k+1>0) ? R-C+k+1 : 0; //number of all-X rows
        if (R0>0) { cblas_dcopy(R0*C*SH,&z,0,&Y[n],1); n += R0*C*SH; }
        for (int r=R0; r<R-RX; r++)
        {
            cblas_dcopy((r+k+1)*SH,&X[n],1,&Y[n],1); n += (r+k+1)*SH;
            cblas_dcopy((C-r-k-1)*SH,&z,0,&Y[n],1); n += (C-r-k-1)*SH;
        }
        if (RX>0) { cblas_dcopy(RX*C*SH,&X[n],1,&Y[n],1); n += RX*C*SH; }
    }

    return 0;
}


int tril_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in tril_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in tril_c: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in tril_c: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in tril_c: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in tril_c: k must be in [1-R C-1]\n"); return 1; }

    const float z[2] = {0.0f,0.0f};
    int n = 0;

    if (iscolmajor)
    {
        const int C0 = (C-R-k>0) ? C-R-k : 0; //number of all-0 cols
        const int CX = (k+1>0) ? k+1 : 0;     //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (CX>0) { cblas_ccopy(R*CX,&X[n],1,&Y[n],1); n += 2*R*CX; }
                for (int c=CX; c<C-C0; c++)
                {
                    cblas_ccopy(c-k,z,0,&Y[n],1); n += 2*(c-k);
                    cblas_ccopy(R-c+k,&X[n],1,&Y[n],1); n += 2*(R-c+k);
                }
                if (C0>0) { cblas_ccopy(R*C0,z,0,&Y[n],1); n += 2*R*C0; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (-k>0) ? -k : 0;           //number of all-0 rows
        const int RX = (R-C+k+1>0) ? R-C+k+1 : 0; //number of all-X rows
        if (R0>0) { cblas_ccopy(R0*C*SH,z,0,&Y[n],1); n += 2*R0*C*SH; }
        for (int r=R0; r<R-RX; r++)
        {
            cblas_ccopy((r+k+1)*SH,&X[n],1,&Y[n],1); n += 2*(r+k+1)*SH;
            cblas_ccopy((C-r-k-1)*SH,z,0,&Y[n],1); n += 2*(C-r-k-1)*SH;
        }
        if (RX>0) { cblas_ccopy(RX*C*SH,&X[n],1,&Y[n],1); n += 2*RX*C*SH; }
    }

    return 0;
}


int tril_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in tril_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in tril_z: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in tril_z: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in tril_z: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in tril_z: k must be in [1-R C-1]\n"); return 1; }

    const double z[2] = {0.0,0.0};
    int n = 0;

    if (iscolmajor)
    {
        const int C0 = (C-R-k>0) ? C-R-k : 0; //number of all-0 cols
        const int CX = (k+1>0) ? k+1 : 0;     //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (CX>0) { cblas_zcopy(R*CX,&X[n],1,&Y[n],1); n += 2*R*CX; }
                for (int c=CX; c<C-C0; c++)
                {
                    cblas_zcopy(c-k,z,0,&Y[n],1); n += 2*(c-k);
                    cblas_zcopy(R-c+k,&X[n],1,&Y[n],1); n += 2*(R-c+k);
                }
                if (C0>0) { cblas_zcopy(R*C0,z,0,&Y[n],1); n += 2*R*C0; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (-k>0) ? -k : 0;           //number of all-0 rows
        const int RX = (R-C+k+1>0) ? R-C+k+1 : 0; //number of all-X rows
        if (R0>0) { cblas_zcopy(R0*C*SH,z,0,&Y[n],1); n += 2*R0*C*SH; }
        for (int r=R0; r<R-RX; r++)
        {
            cblas_zcopy((r+k+1)*SH,&X[n],1,&Y[n],1); n += 2*(r+k+1)*SH;
            cblas_zcopy((C-r-k-1)*SH,z,0,&Y[n],1); n += 2*(C-r-k-1)*SH;
        }
        if (RX>0) { cblas_zcopy(RX*C*SH,&X[n],1,&Y[n],1); n += 2*RX*C*SH; }
    }

    return 0;
}


int tril_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in tril_inplace_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in tril_inplace_s: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in tril_inplace_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in tril_inplace_s: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in tril_inplace_s: k must be in [1-R C-1]\n"); return 1; }

    const float z = 0.0f;
    int n = 0;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    if (iscolmajor)
    {
        const int C0 = (C-R-k>0) ? C-R-k : 0; //number of all-0 cols
        const int CX = (k+1>0) ? k+1 : 0;     //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (CX>0) { n += R*CX; }
                for (int c=CX; c<C-C0; c++)
                {
                    cblas_scopy(c-k,&z,0,&X[n],1);
                    n += R;
                }
                if (C0>0) { cblas_scopy(R*C0,&z,0,&X[n],1); n += R*C0; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (-k>0) ? -k : 0;           //number of all-0 rows
        const int RX = (R-C+k+1>0) ? R-C+k+1 : 0; //number of all-X rows
        if (R0>0) { cblas_scopy(R0*C*SH,&z,0,&X[n],1); n += R0*C*SH; }
        for (int r=R0; r<R-RX; r++)
        {
            n += (r+k+1)*SH;
            cblas_scopy((C-r-k-1)*SH,&z,0,&X[n],1);
            n += SH;
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int tril_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in tril_inplace_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in tril_inplace_d: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in tril_inplace_d: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in tril_inplace_d: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in tril_inplace_d: k must be in [1-R C-1]\n"); return 1; }

    const double z = 0.0;
    int n = 0;

    if (iscolmajor)
    {
        const int C0 = (C-R-k>0) ? C-R-k : 0; //number of all-0 cols
        const int CX = (k+1>0) ? k+1 : 0;     //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (CX>0) { n += R*CX; }
                for (int c=CX; c<C-C0; c++)
                {
                    cblas_dcopy(c-k,&z,0,&X[n],1);
                    n += R;
                }
                if (C0>0) { cblas_dcopy(R*C0,&z,0,&X[n],1); n += R*C0; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (-k>0) ? -k : 0;           //number of all-0 rows
        const int RX = (R-C+k+1>0) ? R-C+k+1 : 0; //number of all-X rows
        if (R0>0) { cblas_dcopy(R0*C*SH,&z,0,&X[n],1); n += R0*C*SH; }
        for (int r=R0; r<R-RX; r++)
        {
            n += (r+k+1)*SH;
            cblas_dcopy((C-r-k-1)*SH,&z,0,&X[n],1);
            n += (C-r-k-1)*SH;
        }
    }

    return 0;
}


int tril_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in tril_inplace_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in tril_inplace_c: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in tril_inplace_c: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in tril_inplace_c: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in tril_inplace_c: k must be in [1-R C-1]\n"); return 1; }

    const float z[2] = {0.0f,0.0f};
    int n = 0;

    if (iscolmajor)
    {
        const int C0 = (C-R-k>0) ? C-R-k : 0; //number of all-0 cols
        const int CX = (k+1>0) ? k+1 : 0;     //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (CX>0) { n += 2*R*CX; }
                for (int c=CX; c<C-C0; c++)
                {
                    cblas_ccopy(c-k,z,0,&X[n],1);
                    n += 2*R;
                }
                if (C0>0) { cblas_ccopy(R*C0,z,0,&X[n],1); n += 2*R*C0; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (-k>0) ? -k : 0;           //number of all-0 rows
        const int RX = (R-C+k+1>0) ? R-C+k+1 : 0; //number of all-X rows
        if (R0>0) { cblas_ccopy(R0*C*SH,z,0,&X[n],1); n += 2*R0*C*SH; }
        for (int r=R0; r<R-RX; r++)
        {
            n += 2*(r+k+1)*SH;
            cblas_ccopy((C-r-k-1)*SH,z,0,&X[n],1);
            n += 2*(C-r-k-1)*SH;
        }
    }

    return 0;
}


int tril_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in tril_inplace_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in tril_inplace_z: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in tril_inplace_z: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in tril_inplace_z: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in tril_inplace_z: k must be in [1-R C-1]\n"); return 1; }

    const double z[2] = {0.0,0.0};
    int n = 0;

    if (iscolmajor)
    {
        const int C0 = (C-R-k>0) ? C-R-k : 0; //number of all-0 cols
        const int CX = (k+1>0) ? k+1 : 0;     //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (CX>0) { n += 2*R*CX; }
                for (int c=CX; c<C-C0; c++)
                {
                    cblas_zcopy(c-k,z,0,&X[n],1);
                    n += 2*R;
                }
                if (C0>0) { cblas_zcopy(R*C0,z,0,&X[n],1); n += 2*R*C0; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (-k>0) ? -k : 0;           //number of all-0 rows
        const int RX = (R-C+k+1>0) ? R-C+k+1 : 0; //number of all-X rows
        if (R0>0) { cblas_zcopy(R0*C*SH,z,0,&X[n],1); n += 2*R0*C*SH; }
        for (int r=R0; r<R-RX; r++)
        {
            n += 2*(r+k+1)*SH;
            cblas_zcopy((C-r-k-1)*SH,z,0,&X[n],1);
            n += 2*(C-r-k-1)*SH;
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
