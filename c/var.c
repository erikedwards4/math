//Gets variance of each row or col of X according to dim.
//For complex case, output is real.
//This works in place.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int var_s (float *Y, float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased);
int var_d (double *Y, double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased);
int var_c (float *Y, float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased);
int var_z (double *Y, double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased);


int var_s (float *Y, float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased)
{
    const float z = 0.0f, o = 1.0f;
    float *x1;
    int r, s, h, l, m, n1 = 0, n2 = 0;
    struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (R<1) { fprintf(stderr,"error in mean_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in mean_s: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in mean_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in mean_s: H (num hyperslices X) must be positive\n"); return 1; }

    //Set ints
    const int RC = R*C, SH = S*H;
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = RC*SH/(M*N);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    const float ni = 1.0f / N;
    const float den = (biased) ? 1.0f/N : 1.0f/(N-1);

    if (N==1) { cblas_scopy(RC*SH,&z,0,Y,1); }
    else if (N==RC*SH)
    {
        Y[0] = 0.0f;
        while (n1<N) { Y[0] += X[n1]; n1++; }
        Y[0] *= ni;
        cblas_saxpy(N,-Y[0],&o,0,X,1);
        Y[0] = 0.0f;
        while (n2<N) { Y[0] += X[n2]*X[n2]; n2++; }
        Y[0] *= den;
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? 1 : N;
        const int ldb = (iscolmajor) ? RC/N : 1;
        const int ldc = (iscolmajor) ? R : C;
        if (!(x1=(float *)malloc((size_t)N*sizeof(float)))) { fprintf(stderr,"error in mean_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(N,&ni,0,x1,1);
        cblas_sgemv(Ord,Tr,R,C,1.0f,X,ldc,x1,1,0.0f,Y,1);
        clock_gettime(CLOCK_REALTIME,&tic);
        for (n1=0, n2=0; n1<RC; n1++, n2++) { if (n2==N) { n2=0; } X[n1] -= Y[n2]; }
        //cblas_scopy(N,&o,0,x1,1);
        //if (dim==0) { cblas_sgemm(Ord,Tr,Tr,R,C,1,-1.0f,x1,lda,Y,ldb,1.0f,X,ldc); }
        //else { cblas_sgemm(Ord,Tr,Tr,R,C,1,-1.0f,Y,ldb,x1,lda,1.0f,X,ldc); }
        clock_gettime(CLOCK_REALTIME,&toc);
        for (n1=0; n1<RC; n1++) { X[n1] *= X[n1]; }
        cblas_sgemv(Ord,Tr,R,C,den,X,ldc,x1,1,0.0f,Y,1);
        free(x1);
    }
    else
    {
        if (iscolmajor && dim<2)
        {
            const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
            if (!(x1=(float *)malloc((size_t)N*sizeof(float)))) { fprintf(stderr,"error in mean_s: problem with malloc. "); perror("malloc"); return 1; }
            cblas_scopy(N,&ni,0,x1,1);
            for (h=0; h<H; h++)
            {
                for (s=0; s<S; s++, n1+=RC, n2+=RC/N)
                {
                    cblas_sgemv(CblasColMajor,Tr,R,C,1.0f,&X[n1],R,x1,1,0.0f,&Y[n2],1);
                }
            }
            free(x1);
        }
        else if (!iscolmajor && dim==2 && H==1)
        {
            if (!(x1=(float *)malloc((size_t)N*sizeof(float)))) { fprintf(stderr,"error in mean_s: problem with malloc. "); perror("malloc"); return 1; }
            cblas_scopy(N,&ni,0,x1,1);
            for (r=0; r<R; r++, n1+=C*SH, n2+=C)
            {
                cblas_sgemv(CblasRowMajor,CblasNoTrans,C,S,1.0f,&X[n1],S,x1,1,0.0f,&Y[n2],1);
            }
            free(x1);
        }
        else
        {
            for (l=0; l<L; l++)
            {
                for (m=0, n1=l*M*N; m<M; m++, n1+=J, n2++)
                {
                    Y[n2] = cblas_sdot(N,&X[n1],I,&ni,0);
                }
            }
        }

        //Subtract mean and get nrm2
        for (l=0, n1=0, n2=0; l<L; l++)
        {
            for (m=0, n1=l*M*N; m<M; m++, n1+=J, n2++)
            {
                cblas_saxpy(N,-Y[n2],&o,0,&X[n1],I);
                Y[n2] = cblas_snrm2(N,&X[n1],I) * den;
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int var_d (double *Y, double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased)
{
    const double o = 1.0;
    double m, den;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in var_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in var_d: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in var_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in var_s: H (num hyperslices X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (biased) { den = (double)R; } else { den = (double)(R-1); }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                m = -cblas_ddot(R,&X[c*R],1,&o,0) / R;
                cblas_daxpy(R,m,&o,0,&X[c*R],1);
                Y[c] = cblas_dnrm2(R,&X[c*R],1) / den;
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                m = -cblas_ddot(R,&X[c],C,&o,0) / R;
                cblas_daxpy(R,m,&o,0,&X[c],C);
                Y[c] = cblas_dnrm2(R,&X[c],C) / den;
            }
        }
    }
    else if (dim==1)
    {
        if (biased) { den = (double)C; } else { den = (double)(C-1); }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                m = -cblas_ddot(C,&X[r],R,&o,0) / C;
                cblas_daxpy(C,m,&o,0,&X[r],R);
                Y[r] = cblas_dnrm2(C,&X[r],R) / den;
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                m = -cblas_ddot(C,&X[r*C],1,&o,0) / C;
                cblas_daxpy(C,m,&o,0,&X[r*C],1);
                Y[r] = cblas_dnrm2(C,&X[r*C],1) / den;
            }
        }
    }
    else
    {
        fprintf(stderr,"error in var_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int var_c (float *Y, float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased)
{
    const float o[2] =  {1.0f,0.0f};
    _Complex float m;
    float den;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in var_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in var_c: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in var_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in var_s: H (num hyperslices X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (biased) { den = (float)R; } else { den = (float)(R-1); }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                m = -cblas_cdotu(R,&X[2*c*R],1,&o[0],0) / R;
                cblas_caxpy(R,(float *)&m,&o[0],0,&X[2*c*R],1);
                Y[c] = cblas_scnrm2(R,&X[2*c*R],1) / den;
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                m = -cblas_cdotu(R,&X[2*c],C,&o[0],0) / R;
                cblas_caxpy(R,(float *)&m,&o[0],0,&X[2*c],C);
                Y[c] = cblas_scnrm2(R,&X[2*c],C) / den;
            }
        }
    }
    else if (dim==1)
    {
        if (biased) { den = (float)C; } else { den = (float)(C-1); }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                m = -cblas_cdotu(C,&X[2*r],R,&o[0],0) / C;
                cblas_caxpy(C,(float *)&m,&o[0],0,&X[2*r],R);
                Y[r] = cblas_scnrm2(C,&X[2*r],R) / den;
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                m = -cblas_cdotu(C,&X[2*r*C],1,&o[0],0) / C;
                cblas_caxpy(C,(float *)&m,&o[0],0,&X[2*r*C],1);
                Y[r] = cblas_scnrm2(C,&X[2*r*C],1) / den;
            }
        }
    }
    else
    {
        fprintf(stderr,"error in var_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int var_z (double *Y, double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased)
{
    const double o[2] =  {1.0,0.0};
    _Complex double m;
    double den;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in var_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in var_z: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in var_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in var_s: H (num hyperslices X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (biased) { den = (double)R; } else { den = (double)(R-1); }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                m = -cblas_zdotu(R,&X[2*c*R],1,&o[0],0) / R;
                cblas_zaxpy(R,(double *)&m,&o[0],0,&X[2*c*R],1);
                Y[c] = cblas_dznrm2(R,&X[2*c*R],1) / den;
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                m = -cblas_zdotu(R,&X[2*c],C,&o[0],0) / R;
                cblas_zaxpy(R,(double *)&m,&o[0],0,&X[2*c],C);
                Y[c] = cblas_dznrm2(R,&X[2*c],C) / den;
            }
        }
    }
    else if (dim==1)
    {
        if (biased) { den = (double)C; } else { den = (double)(C-1); }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                m = -cblas_zdotu(C,&X[2*r],R,&o[0],0) / C;
                cblas_zaxpy(C,(double *)&m,&o[0],0,&X[2*r],R);
                Y[r] = cblas_dznrm2(C,&X[2*r],R) / den;
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                m = -cblas_zdotu(C,&X[2*r*C],1,&o[0],0) / C;
                cblas_zaxpy(C,(double *)&m,&o[0],0,&X[2*r*C],1);
                Y[r] = cblas_dznrm2(C,&X[2*r*C],1) / den;
            }
        }
    }
    else
    {
        fprintf(stderr,"error in var_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
