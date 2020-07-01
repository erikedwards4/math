//Zeros all elements below the kth diagonal of matrix X.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int triu_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);
int triu_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);
int triu_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);
int triu_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);

int triu_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);
int triu_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);
int triu_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);
int triu_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k);


int triu_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in triu_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in triu_s: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in triu_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in triu_s: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in triu_s: k must be in [1-R C-1]\n"); return 1; }

    const float z = 0.0f;
    int n = 0;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    if (iscolmajor)
    {
        const int C0 = (k>0) ? k : 0;         //number of all-0 cols
        const int CX = (k>C-R) ? 0 : C-R+1-k; //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (C0>0) { cblas_scopy(R*C0,&z,0,&Y[n],1); n += R*C0; }
                for (int c=C0; c<C-CX; c++)
                {
                    cblas_scopy(c-k+1,&X[n],1,&Y[n],1); n += c-k+1;
                    cblas_scopy(R-c+k-1,&z,0,&Y[n],1); n += R-c+k-1;
                }
                if (CX>0) { cblas_scopy(R*CX,&X[n],1,&Y[n],1); n += R*CX; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (k>C-R) ? R-C+k : 0;  //number of all-0 rows
        const int RX = (k>0) ? 0 : k+1;      //number of all-X rows
        if (RX>0) { cblas_scopy(RX*C*SH,&X[n],1,&Y[n],1); n += RX*C*SH; }
        for (int r=RX; r<R-R0; r++)
        {
            cblas_scopy((r+k)*SH,&z,0,&Y[n],1); n += (r+k)*SH;
            cblas_scopy((C-r-k)*SH,&X[n],0,&Y[n],1); n += (C-r-k)*SH;
        }
        if (R0>0) { cblas_scopy(R0*C*SH,&z,0,&Y[n],1); n += R0*C*SH; }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int triu_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in triu_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in triu_d: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in triu_d: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in triu_d: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in triu_d: k must be in [1-R C-1]\n"); return 1; }

    const double z = 0.0;
    int n = 0;

    if (iscolmajor)
    {
        const int C0 = (k>0) ? k : 0;         //number of all-0 cols
        const int CX = (k>C-R) ? 0 : C-R+1-k; //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (C0>0) { cblas_dcopy(R*C0,&z,0,&Y[n],1); n += R*C0; }
                for (int c=C0; c<C-CX; c++)
                {
                    cblas_dcopy(c-k+1,&X[n],1,&Y[n],1); n += c-k+1;
                    cblas_dcopy(R-c+k-1,&z,0,&Y[n],1); n += R-c+k-1;
                }
                if (CX>0) { cblas_dcopy(R*CX,&X[n],1,&Y[n],1); n += R*CX; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (k>C-R) ? R-C+k : 0;  //number of all-0 rows
        const int RX = (k>0) ? 0 : k+1;      //number of all-X rows
        if (RX>0) { cblas_dcopy(RX*C*SH,&X[n],1,&Y[n],1); n += RX*C*SH; }
        for (int r=RX; r<R-R0; r++)
        {
            cblas_dcopy((r+k)*SH,&z,0,&Y[n],1); n += (r+k)*SH;
            cblas_dcopy((C-r-k)*SH,&X[n],0,&Y[n],1); n += (C-r-k)*SH;
        }
        if (R0>0) { cblas_dcopy(R0*C*SH,&z,0,&Y[n],1); n += R0*C*SH; }
    }

    return 0;
}


int triu_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in triu_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in triu_c: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in triu_c: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in triu_c: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in triu_c: k must be in [1-R C-1]\n"); return 1; }

    const float z[2] = {0.0f,0.0f};
    int n = 0;

    if (iscolmajor)
    {
        const int C0 = (k>0) ? k : 0;         //number of all-0 cols
        const int CX = (k>C-R) ? 0 : C-R+1-k; //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (C0>0) { cblas_ccopy(R*C0,z,0,&Y[n],1); n += 2*R*C0; }
                for (int c=C0; c<C-CX; c++)
                {
                    cblas_ccopy(c-k+1,&X[n],1,&Y[n],1); n += 2*(c-k+1);
                    cblas_ccopy(R-c+k-1,z,0,&Y[n],1); n += 2*(R-c+k-1);
                }
                if (CX>0) { cblas_ccopy(R*CX,&X[n],1,&Y[n],1); n += 2*R*CX; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (k>C-R) ? R-C+k : 0;  //number of all-0 rows
        const int RX = (k>0) ? 0 : k+1;      //number of all-X rows
        if (RX>0) { cblas_ccopy(RX*C*SH,&X[n],1,&Y[n],1); n += 2*RX*C*SH; }
        for (int r=RX; r<R-R0; r++)
        {
            cblas_ccopy((r+k)*SH,z,0,&Y[n],1); n += 2*(r+k)*SH;
            cblas_ccopy((C-r-k)*SH,&X[n],0,&Y[n],1); n += 2*(C-r-k)*SH;
        }
        if (R0>0) { cblas_ccopy(R0*C*SH,z,0,&Y[n],1); n += 2*R0*C*SH; }
    }

    return 0;
}


int triu_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in triu_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in triu_z: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in triu_z: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in triu_z: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in triu_z: k must be in [1-R C-1]\n"); return 1; }

    const double z[2] = {0.0,0.0};
    int n = 0;

    if (iscolmajor)
    {
        const int C0 = (k>0) ? k : 0;         //number of all-0 cols
        const int CX = (k>C-R) ? 0 : C-R+1-k; //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (C0>0) { cblas_zcopy(R*C0,z,0,&Y[n],1); n += 2*R*C0; }
                for (int c=C0; c<C-CX; c++)
                {
                    cblas_zcopy(c-k+1,&X[n],1,&Y[n],1); n += 2*(c-k+1);
                    cblas_zcopy(R-c+k-1,z,0,&Y[n],1); n += 2*(R-c+k-1);
                }
                if (CX>0) { cblas_zcopy(R*CX,&X[n],1,&Y[n],1); n += 2*R*CX; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (k>C-R) ? R-C+k : 0;  //number of all-0 rows
        const int RX = (k>0) ? 0 : k+1;      //number of all-X rows
        if (RX>0) { cblas_zcopy(RX*C*SH,&X[n],1,&Y[n],1); n += 2*RX*C*SH; }
        for (int r=RX; r<R-R0; r++)
        {
            cblas_zcopy((r+k)*SH,z,0,&Y[n],1); n += 2*(r+k)*SH;
            cblas_zcopy((C-r-k)*SH,&X[n],0,&Y[n],1); n += 2*(C-r-k)*SH;
        }
        if (R0>0) { cblas_zcopy(R0*C*SH,z,0,&Y[n],1); n += 2*R0*C*SH; }
    }

    return 0;
}


int triu_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in triu_inplace_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in triu_inplace_s: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in triu_inplace_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in triu_inplace_s: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in triu_inplace_s: k must be in [1-R C-1]\n"); return 1; }

    const float z = 0.0f;
    int n = 0;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    if (iscolmajor)
    {
        const int C0 = (k>0) ? k : 0;         //number of all-0 cols
        const int CX = (k>C-R) ? 0 : C-R+1-k; //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (C0>0) { cblas_scopy(R*C0,&z,0,&X[n],1); n += R*C0; }
                for (int c=C0; c<C-CX; c++)
                {
                    n += c-k+1;
                    cblas_scopy(R-c+k-1,&z,0,&X[n],1);
                    n += R-c+k-1;
                }
                if (CX>0) { n += R*CX; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (k>C-R) ? R-C+k : 0;  //number of all-0 rows
        const int RX = (k>0) ? 0 : k+1;      //number of all-X rows
        if (RX>0) { n += RX*C*SH; }
        for (int r=RX; r<R-R0; r++)
        {
            cblas_scopy((r+k)*SH,&z,0,&X[n],1);
            n += C*SH;
        }
        if (R0>0) { cblas_scopy(R0*C*SH,&z,0,&X[n],1); n += R0*C*SH; }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int triu_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in triu_inplace_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in triu_inplace_d: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in triu_inplace_d: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in triu_inplace_d: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in triu_inplace_d: k must be in [1-R C-1]\n"); return 1; }

    const double z = 0.0;
    int n = 0;

    if (iscolmajor)
    {
        const int C0 = (k>0) ? k : 0;         //number of all-0 cols
        const int CX = (k>C-R) ? 0 : C-R+1-k; //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (C0>0) { cblas_dcopy(R*C0,&z,0,&X[n],1); n += R*C0; }
                for (int c=C0; c<C-CX; c++)
                {
                    n += c-k+1;
                    cblas_dcopy(R-c+k-1,&z,0,&X[n],1);
                    n += R-c+k-1;
                }
                if (CX>0) { n += R*CX; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (k>C-R) ? R-C+k : 0;  //number of all-0 rows
        const int RX = (k>0) ? 0 : k+1;      //number of all-X rows
        if (RX>0) { n += RX*C*SH; }
        for (int r=RX; r<R-R0; r++)
        {
            cblas_dcopy((r+k)*SH,&z,0,&X[n],1);
            n += C*SH;
        }
        if (R0>0) { cblas_dcopy(R0*C*SH,&z,0,&X[n],1); n += R0*C*SH; }
    }

    return 0;
}


int triu_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in triu_inplace_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in triu_inplace_c: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in triu_inplace_c: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in triu_inplace_c: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in triu_inplace_c: k must be in [1-R C-1]\n"); return 1; }

    const float z[2] = {0.0f,0.0f};
    int n = 0;

    if (iscolmajor)
    {
        const int C0 = (k>0) ? k : 0;         //number of all-0 cols
        const int CX = (k>C-R) ? 0 : C-R+1-k; //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (C0>0) { cblas_ccopy(R*C0,z,0,&X[n],1); n += 2*R*C0; }
                for (int c=C0; c<C-CX; c++)
                {
                    n += 2*(c-k+1);
                    cblas_ccopy(R-c+k-1,z,0,&X[n],1);
                    n += 2*(R-c+k-1);
                }
                if (CX>0) { n += 2*R*CX; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (k>C-R) ? R-C+k : 0;  //number of all-0 rows
        const int RX = (k>0) ? 0 : k+1;      //number of all-X rows
        if (RX>0) { n += 2*RX*C*SH; }
        for (int r=RX; r<R-R0; r++)
        {
            cblas_ccopy((r+k)*SH,z,0,&X[n],1);
            n += 2*C*SH;
        }
        if (R0>0) { cblas_ccopy(R0*C*SH,z,0,&X[n],1); n += 2*R0*C*SH; }
    }

    return 0;
}


int triu_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int k)
{
    //Checks
    if (R<1) { fprintf(stderr,"error in triu_inplace_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in triu_inplace_z: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in triu_inplace_z: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in triu_inplace_z: H (num hyperslices X) must be positive\n"); return 1; }
    if (k<=-R || k>=C) { fprintf(stderr,"error in triu_inplace_z: k must be in [1-R C-1]\n"); return 1; }

    const double z[2] = {0.0,0.0};
    int n = 0;

    if (iscolmajor)
    {
        const int C0 = (k>0) ? k : 0;         //number of all-0 cols
        const int CX = (k>C-R) ? 0 : C-R+1-k; //number of all-X cols
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                if (C0>0) { cblas_zcopy(R*C0,z,0,&X[n],1); n += 2*R*C0; }
                for (int c=C0; c<C-CX; c++)
                {
                    n += 2*(c-k+1);
                    cblas_zcopy(R-c+k-1,z,0,&X[n],1);
                    n += 2*(R-c+k-1);
                }
                if (CX>0) { n += 2*R*CX; }
            }
        }
    }
    else
    {
        const int SH = S*H;
        const int R0 = (k>C-R) ? R-C+k : 0;  //number of all-0 rows
        const int RX = (k>0) ? 0 : k+1;      //number of all-X rows
        if (RX>0) { n += 2*RX*C*SH; }
        for (int r=RX; r<R-R0; r++)
        {
            cblas_zcopy((r+k)*SH,z,0,&X[n],1);
            n += 2*C*SH;
        }
        if (R0>0) { cblas_zcopy(R0*C*SH,z,0,&X[n],1); n += 2*R0*C*SH; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
