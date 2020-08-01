//Joins (stacks) 3 inputs X1, X2, X3 into 1 output Y

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int join3_s (float *Y, const float *X1, const float *X2, const float *X3,  const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const size_t R3, const size_t C3, const size_t S3, const size_t H3, const char iscolmajor, const size_t dim);
int join3_d (double *Y, const double *X1, const double *X2, const double *X3,  const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const size_t R3, const size_t C3, const size_t S3, const size_t H3, const char iscolmajor, const size_t dim);
int join3_c (float *Y, const float *X1, const float *X2, const float *X3,  const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const size_t R3, const size_t C3, const size_t S3, const size_t H3, const char iscolmajor, const size_t dim);
int join3_z (double *Y, const double *X1, const double *X2, const double *X3,  const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const size_t R3, const size_t C3, const size_t S3, const size_t H3, const char iscolmajor, const size_t dim);


int join3_s (float *Y, const float *X1, const float *X2, const float *X3,  const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const size_t R3, const size_t C3, const size_t S3, const size_t H3, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in join3_s: dim must be in [0 3]\n"); return 1; }
    if (dim!=0 && (R1!=R2 || R1!=R3)) { fprintf(stderr,"error in join3_s: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1 && (C1!=C2 || C1!=C3)) { fprintf(stderr,"error in join3_s: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2 && (S1!=S2 || S1!=S3)) { fprintf(stderr,"error in join3_s: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3 && (H1!=H2 || H1!=H3)) { fprintf(stderr,"error in join3_s: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0) ? R1 : (dim==1) ? R1*C1 : (dim==2) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0) ? H1*S1*C1*R1 : (dim==1) ? H1*S1*C1 : (dim==2) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0) ? R2 : (dim==1) ? R2*C2 : (dim==2) ? R2*C2*S2 : R2*C2*S2*H2) : ((dim==0) ? H2*S2*C2*R2 : (dim==1) ? H2*S2*C2 : (dim==2) ? H2*S2 : H2);
    const size_t M3 = (iscolmajor) ? ((dim==0) ? R3 : (dim==1) ? R3*C3 : (dim==2) ? R3*C3*S3 : R3*C3*S3*H3) : ((dim==0) ? H3*S3*C3*R3 : (dim==1) ? H3*S3*C3 : (dim==2) ? H3*S3 : H3);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=0; g<G; ++g)
    {
        cblas_scopy((int)B1,X1,1,Y,1); X1 += B1; Y += B1;
        cblas_scopy((int)B2,X2,1,Y,1); X2 += B2; Y += B2;
        cblas_scopy((int)M3,X3,1,Y,1); X3 += M3; Y += M3;
    }

    return 0;
}


int join3_d (double *Y, const double *X1, const double *X2, const double *X3,  const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const size_t R3, const size_t C3, const size_t S3, const size_t H3, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in join3_d: dim must be in [0 3]\n"); return 1; }
    if (dim!=0 && (R1!=R2 || R1!=R3)) { fprintf(stderr,"error in join3_d: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1 && (C1!=C2 || C1!=C3)) { fprintf(stderr,"error in join3_d: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2 && (S1!=S2 || S1!=S3)) { fprintf(stderr,"error in join3_d: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3 && (H1!=H2 || H1!=H3)) { fprintf(stderr,"error in join3_d: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0) ? R1 : (dim==1) ? R1*C1 : (dim==2) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0) ? H1*S1*C1*R1 : (dim==1) ? H1*S1*C1 : (dim==2) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0) ? R2 : (dim==1) ? R2*C2 : (dim==2) ? R2*C2*S2 : R2*C2*S2*H2) : ((dim==0) ? H2*S2*C2*R2 : (dim==1) ? H2*S2*C2 : (dim==2) ? H2*S2 : H2);
    const size_t M3 = (iscolmajor) ? ((dim==0) ? R3 : (dim==1) ? R3*C3 : (dim==2) ? R3*C3*S3 : R3*C3*S3*H3) : ((dim==0) ? H3*S3*C3*R3 : (dim==1) ? H3*S3*C3 : (dim==2) ? H3*S3 : H3);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=0; g<G; ++g)
    {
        cblas_dcopy((int)B1,X1,1,Y,1); X1 += B1; Y += B1;
        cblas_dcopy((int)B2,X2,1,Y,1); X2 += B2; Y += B2;
        cblas_dcopy((int)M3,X3,1,Y,1); X3 += M3; Y += M3;
    }

    return 0;
}


int join3_c (float *Y, const float *X1, const float *X2, const float *X3,  const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const size_t R3, const size_t C3, const size_t S3, const size_t H3, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in join3_c: dim must be in [0 3]\n"); return 1; }
    if (dim!=0 && (R1!=R2 || R1!=R3)) { fprintf(stderr,"error in join3_c: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1 && (C1!=C2 || C1!=C3)) { fprintf(stderr,"error in join3_c: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2 && (S1!=S2 || S1!=S3)) { fprintf(stderr,"error in join3_c: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3 && (H1!=H2 || H1!=H3)) { fprintf(stderr,"error in join3_c: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0) ? R1 : (dim==1) ? R1*C1 : (dim==2) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0) ? H1*S1*C1*R1 : (dim==1) ? H1*S1*C1 : (dim==2) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0) ? R2 : (dim==1) ? R2*C2 : (dim==2) ? R2*C2*S2 : R2*C2*S2*H2) : ((dim==0) ? H2*S2*C2*R2 : (dim==1) ? H2*S2*C2 : (dim==2) ? H2*S2 : H2);
    const size_t M3 = (iscolmajor) ? ((dim==0) ? R3 : (dim==1) ? R3*C3 : (dim==2) ? R3*C3*S3 : R3*C3*S3*H3) : ((dim==0) ? H3*S3*C3*R3 : (dim==1) ? H3*S3*C3 : (dim==2) ? H3*S3 : H3);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=0; g<G; ++g)
    {
        cblas_ccopy((int)B1,X1,1,Y,1); X1 += 2*B1; Y += 2*B1;
        cblas_ccopy((int)B2,X2,1,Y,1); X2 += 2*B2; Y += 2*B2;
        cblas_ccopy((int)M3,X3,1,Y,1); X3 += 2*M3; Y += 2*M3;
    }

    return 0;
}


int join3_z (double *Y, const double *X1, const double *X2, const double *X3,  const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const size_t R3, const size_t C3, const size_t S3, const size_t H3, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in join3_z: dim must be in [0 3]\n"); return 1; }
    if (dim!=0 && (R1!=R2 || R1!=R3)) { fprintf(stderr,"error in join3_z: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1 && (C1!=C2 || C1!=C3)) { fprintf(stderr,"error in join3_z: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2 && (S1!=S2 || S1!=S3)) { fprintf(stderr,"error in join3_z: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3 && (H1!=H2 || H1!=H3)) { fprintf(stderr,"error in join3_z: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0) ? R1 : (dim==1) ? R1*C1 : (dim==2) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0) ? H1*S1*C1*R1 : (dim==1) ? H1*S1*C1 : (dim==2) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0) ? R2 : (dim==1) ? R2*C2 : (dim==2) ? R2*C2*S2 : R2*C2*S2*H2) : ((dim==0) ? H2*S2*C2*R2 : (dim==1) ? H2*S2*C2 : (dim==2) ? H2*S2 : H2);
    const size_t M3 = (iscolmajor) ? ((dim==0) ? R3 : (dim==1) ? R3*C3 : (dim==2) ? R3*C3*S3 : R3*C3*S3*H3) : ((dim==0) ? H3*S3*C3*R3 : (dim==1) ? H3*S3*C3 : (dim==2) ? H3*S3 : H3);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=0; g<G; ++g)
    {
        cblas_zcopy((int)B1,X1,1,Y,1); X1 += 2*B1; Y += 2*B1;
        cblas_zcopy((int)B2,X2,1,Y,1); X2 += 2*B2; Y += 2*B2;
        cblas_zcopy((int)M3,X3,1,Y,1); X3 += 2*M3; Y += 2*M3;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
