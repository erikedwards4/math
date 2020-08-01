//Zeros all elements below the kth diagonal of matrix X.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int triu_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int triu_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int triu_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int triu_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);

int triu_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int triu_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int triu_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int triu_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);


int triu_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_s: k must be in [1-R C-1]\n"); return 1; }

    const float z = 0.0f;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;         //number of all-0 cols
        const size_t CX = (k<(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1); //number of all-X cols
        for (size_t h=0, n=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                if (C0>0) { cblas_scopy((int)(R*C0),&z,0,&Y[n],1); n += R*C0; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    cblas_scopy((int)c-k+1,&X[n],1,&Y[n],1); n += (size_t)((int)c-k+1);
                    cblas_scopy((int)R-(int)c+k-1,&z,0,&Y[n],1); n += (size_t)((int)R-(int)c+k-1);
                }
                if (CX>0) { cblas_scopy((int)(R*CX),&X[n],1,&Y[n],1); n += R*CX; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)C-(int)R) ? (size_t)((int)R-(int)C+k) : 0;  //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(k+1);      //number of all-X rows
        if (RX>0) { cblas_scopy((int)(RX*C*SH),X,1,Y,1); }
        for (size_t r=RX, n=RX*C*SH; r<R-R0; ++r)
        {
            cblas_scopy((int)SH*((int)r+k),&z,0,&Y[n],1); n += (size_t)((int)SH*((int)r+k));
            cblas_scopy((int)SH*((int)C-(int)r-k),&X[n],0,&Y[n],1); n += (size_t)((int)SH*((int)C-(int)r-k));
        }
        if (R0>0) { cblas_scopy((int)(R0*C*SH),&z,0,&Y[C*SH*(R-R0)],1); }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int triu_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_d: k must be in [1-R C-1]\n"); return 1; }

    const double z = 0.0;

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;         //number of all-0 cols
        const size_t CX = (k<(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1); //number of all-X cols
        for (size_t h=0, n=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                if (C0>0) { cblas_dcopy((int)(R*C0),&z,0,&Y[n],1); n += R*C0; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    cblas_dcopy((int)c-k+1,&X[n],1,&Y[n],1); n += (size_t)((int)c-k+1);
                    cblas_dcopy((int)R-(int)c+k-1,&z,0,&Y[n],1); n += (size_t)((int)R-(int)c+k-1);
                }
                if (CX>0) { cblas_dcopy((int)(R*CX),&X[n],1,&Y[n],1); n += R*CX; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)C-(int)R) ? (size_t)((int)R-(int)C+k) : 0;  //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(k+1);      //number of all-X rows
        if (RX>0) { cblas_dcopy((int)(RX*C*SH),X,1,Y,1); }
        for (size_t r=RX, n=RX*C*SH; r<R-R0; ++r)
        {
            cblas_dcopy((int)SH*((int)r+k),&z,0,&Y[n],1); n += (size_t)((int)SH*((int)r+k));
            cblas_dcopy((int)SH*((int)C-(int)r-k),&X[n],0,&Y[n],1); n += (size_t)((int)SH*((int)C-(int)r-k));
        }
        if (R0>0) { cblas_dcopy((int)(R0*C*SH),&z,0,&Y[C*SH*(R-R0)],1); }
    }

    return 0;
}


int triu_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_c: k must be in [1-R C-1]\n"); return 1; }

    const float z[2] = {0.0f,0.0f};

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;         //number of all-0 cols
        const size_t CX = (k<(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1); //number of all-X cols
        for (size_t h=0, n=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                if (C0>0) { cblas_ccopy((int)(R*C0),z,0,&Y[n],1); n += 2*R*C0; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    cblas_ccopy((int)c-k+1,&X[n],1,&Y[n],1); n += 2*(size_t)((int)c-k+1);
                    cblas_ccopy((int)R-(int)c+k-1,z,0,&Y[n],1); n += 2*(size_t)((int)R-(int)c+k-1);
                }
                if (CX>0) { cblas_ccopy((int)(R*CX),&X[n],1,&Y[n],1); n += 2*R*CX; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)C-(int)R) ? (size_t)((int)R-(int)C+k) : 0;  //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(k+1);      //number of all-X rows
        if (RX>0) { cblas_ccopy((int)(RX*C*SH),X,1,Y,1); }
        for (size_t r=RX, n=2*RX*C*SH; r<R-R0; ++r)
        {
            cblas_ccopy((int)SH*((int)r+k),z,0,&Y[n],1); n += 2*(size_t)((int)SH*((int)r+k));
            cblas_ccopy((int)SH*((int)C-(int)r-k),&X[n],0,&Y[n],1); n += 2*(size_t)((int)SH*((int)C-(int)r-k));
        }
        if (R0>0) { cblas_ccopy((int)(R0*C*SH),z,0,&Y[2*C*SH*(R-R0)],1); }
    }

    return 0;
}


int triu_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_z: k must be in [1-R C-1]\n"); return 1; }

    const double z[2] = {0.0,0.0};

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;         //number of all-0 cols
        const size_t CX = (k<(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1); //number of all-X cols
        for (size_t h=0, n=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                if (C0>0) { cblas_zcopy((int)(R*C0),z,0,&Y[n],1); n += 2*R*C0; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    cblas_zcopy((int)c-k+1,&X[n],1,&Y[n],1); n += 2*(size_t)((int)c-k+1);
                    cblas_zcopy((int)R-(int)c+k-1,z,0,&Y[n],1); n += 2*(size_t)((int)R-(int)c+k-1);
                }
                if (CX>0) { cblas_zcopy((int)(R*CX),&X[n],1,&Y[n],1); n += 2*R*CX; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)C-(int)R) ? (size_t)((int)R-(int)C+k) : 0;  //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(k+1);      //number of all-X rows
        if (RX>0) { cblas_zcopy((int)(RX*C*SH),X,1,Y,1); }
        for (size_t r=RX, n=2*RX*C*SH; r<R-R0; ++r)
        {
            cblas_zcopy((int)SH*((int)r+k),z,0,&Y[n],1); n += 2*(size_t)((int)SH*((int)r+k));
            cblas_zcopy((int)SH*((int)C-(int)r-k),&X[n],0,&Y[n],1); n += 2*(size_t)((int)SH*((int)C-(int)r-k));
        }
        if (R0>0) { cblas_zcopy((int)(R0*C*SH),z,0,&Y[2*C*SH*(R-R0)],1); }
    }

    return 0;
}


int triu_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_inplace_s: k must be in [1-R C-1]\n"); return 1; }

    const float z = 0.0f;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;         //number of all-0 cols
        const size_t CX = (k<(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1); //number of all-X cols
        for (size_t h=0, n=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                if (C0>0) { cblas_scopy((int)(R*C0),&z,0,&X[n],1); n += R*C0; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    n += (size_t)((int)c-k+1);
                    cblas_scopy((int)R-(int)c+k-1,&z,0,&X[n],1);
                    n += (size_t)((int)R-(int)c+k-1);
                }
                if (CX>0) { n += R*CX; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)C-(int)R) ? (size_t)((int)R-(int)C+k) : 0;  //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(k+1);      //number of all-X rows
        for (size_t r=RX, n=RX*C*SH; r<R-R0; ++r)
        {
            cblas_scopy((int)SH*((int)r+k),&z,0,&X[n],1);
            n += C*SH;
        }
        if (R0>0) { cblas_scopy((int)(R0*C*SH),&z,0,&X[C*SH*(R-R0)],1); }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int triu_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_inplace_d: k must be in [1-R C-1]\n"); return 1; }

    const double z = 0.0;

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;         //number of all-0 cols
        const size_t CX = (k<(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1); //number of all-X cols
        for (size_t h=0, n=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                if (C0>0) { cblas_dcopy((int)(R*C0),&z,0,&X[n],1); n += R*C0; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    n += (size_t)((int)c-k+1);
                    cblas_dcopy((int)R-(int)c+k-1,&z,0,&X[n],1);
                    n += (size_t)((int)R-(int)c+k-1);
                }
                if (CX>0) { n += R*CX; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)C-(int)R) ? (size_t)((int)R-(int)C+k) : 0;  //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(k+1);      //number of all-X rows
        for (size_t r=RX, n=RX*C*SH; r<R-R0; ++r)
        {
            cblas_dcopy((int)SH*((int)r+k),&z,0,&X[n],1);
            n += C*SH;
        }
        if (R0>0) { cblas_dcopy((int)(R0*C*SH),&z,0,&X[C*SH*(R-R0)],1); }
    }

    return 0;
}


int triu_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_inplace_c: k must be in [1-R C-1]\n"); return 1; }

    const float z[2] = {0.0f,0.0f};

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;         //number of all-0 cols
        const size_t CX = (k<(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1); //number of all-X cols
        for (size_t h=0, n=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                if (C0>0) { cblas_ccopy((int)(R*C0),z,0,&X[n],1); n += 2*R*C0; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    n += 2*(size_t)((int)c-k+1);
                    cblas_ccopy((int)R-(int)c+k-1,z,0,&X[n],1);
                    n += 2*(size_t)((int)R-(int)c+k-1);
                }
                if (CX>0) { n += 2*R*CX; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)C-(int)R) ? (size_t)((int)R-(int)C+k) : 0;  //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(k+1);      //number of all-X rows
        for (size_t r=RX, n=2*RX*C*SH; r<R-R0; ++r)
        {
            cblas_ccopy((int)SH*((int)r+k),z,0,&X[n],1);
            n += 2*C*SH;
        }
        if (R0>0) { cblas_ccopy((int)(R0*C*SH),z,0,&X[2*C*SH*(R-R0)],1); }
    }

    return 0;
}


int triu_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in triu_inplace_z: k must be in [1-R C-1]\n"); return 1; }

    const double z[2] = {0.0,0.0};

    if (iscolmajor)
    {
        const size_t C0 = (k>0) ? (size_t)k : 0;         //number of all-0 cols
        const size_t CX = (k<(int)C-(int)R+1) ? 0 : (size_t)((int)C-(int)R-k+1); //number of all-X cols
        for (size_t h=0, n=0; h<H; ++h)
        {
            for (size_t s=0; s<S; ++s)
            {
                if (C0>0) { cblas_zcopy((int)(R*C0),z,0,&X[n],1); n += 2*R*C0; }
                for (size_t c=C0; c<C-CX; ++c)
                {
                    n += 2*(size_t)((int)c-k+1);
                    cblas_zcopy((int)R-(int)c+k-1,z,0,&X[n],1);
                    n += 2*(size_t)((int)R-(int)c+k-1);
                }
                if (CX>0) { n += 2*R*CX; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<(int)C-(int)R) ? (size_t)((int)R-(int)C+k) : 0;  //number of all-0 rows
        const size_t RX = (k>0) ? 0 : (size_t)(k+1);      //number of all-X rows
        for (size_t r=RX, n=2*RX*C*SH; r<R-R0; ++r)
        {
            cblas_zcopy((int)SH*((int)r+k),z,0,&X[n],1);
            n += 2*C*SH;
        }
        if (R0>0) { cblas_zcopy((int)(R0*C*SH),z,0,&X[2*C*SH*(R-R0)],1); }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
