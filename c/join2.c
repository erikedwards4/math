//Joins 2 inputs X1, X2 into 1 output Y

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int join2_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int join2_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int join2_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int join2_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);


int join2_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in join2_s: dim must be in [0 3]\n"); return 1; }
    if (dim!=0 && R1!=R2) { fprintf(stderr,"error in join2_s: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1 && C1!=C2) { fprintf(stderr,"error in join2_s: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2 && S1!=S2) { fprintf(stderr,"error in join2_s: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3 && H1!=H2) { fprintf(stderr,"error in join2_s: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0) ? R1 : (dim==1) ? R1*C1 : (dim==2) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0) ? H1*S1*C1*R1 : (dim==1) ? H1*S1*C1 : (dim==2) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0) ? R2 : (dim==1) ? R2*C2 : (dim==2) ? R2*C2*S2 : R2*C2*S2*H2) : ((dim==0) ? H2*S2*C2*R2 : (dim==1) ? H2*S2*C2 : (dim==2) ? H2*S2 : H2);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=0; g<G; ++g)
    {
        cblas_scopy((int)B1,X1,1,Y,1); X1 += B1; Y += B1;
        cblas_scopy((int)B2,X2,1,Y,1); X2 += B2; Y += B2;
    }

    return 0;
}


int join2_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in join2_d: dim must be in [0 3]\n"); return 1; }
    if (dim!=0 && R1!=R2) { fprintf(stderr,"error in join2_d: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1 && C1!=C2) { fprintf(stderr,"error in join2_d: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2 && S1!=S2) { fprintf(stderr,"error in join2_d: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3 && H1!=H2) { fprintf(stderr,"error in join2_d: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0) ? R1 : (dim==1) ? R1*C1 : (dim==2) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0) ? H1*S1*C1*R1 : (dim==1) ? H1*S1*C1 : (dim==2) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0) ? R2 : (dim==1) ? R2*C2 : (dim==2) ? R2*C2*S2 : R2*C2*S2*H2) : ((dim==0) ? H2*S2*C2*R2 : (dim==1) ? H2*S2*C2 : (dim==2) ? H2*S2 : H2);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=0; g<G; ++g)
    {
        cblas_dcopy((int)B1,X1,1,Y,1); X1 += B1; Y += B1;
        cblas_dcopy((int)B2,X2,1,Y,1); X2 += B2; Y += B2;
    }

    return 0;
}


int join2_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in join2_c: dim must be in [0 3]\n"); return 1; }
    if (dim!=0 && R1!=R2) { fprintf(stderr,"error in join2_c: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1 && C1!=C2) { fprintf(stderr,"error in join2_c: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2 && S1!=S2) { fprintf(stderr,"error in join2_c: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3 && H1!=H2) { fprintf(stderr,"error in join2_c: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0) ? R1 : (dim==1) ? R1*C1 : (dim==2) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0) ? H1*S1*C1*R1 : (dim==1) ? H1*S1*C1 : (dim==2) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0) ? R2 : (dim==1) ? R2*C2 : (dim==2) ? R2*C2*S2 : R2*C2*S2*H2) : ((dim==0) ? H2*S2*C2*R2 : (dim==1) ? H2*S2*C2 : (dim==2) ? H2*S2 : H2);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=0; g<G; ++g)
    {
        cblas_ccopy((int)B1,X1,1,Y,1); X1 += 2*B1; Y += 2*B1;
        cblas_ccopy((int)B2,X2,1,Y,1); X2 += 2*B2; Y += 2*B2;
    }

    return 0;
}


int join2_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in join2_z: dim must be in [0 3]\n"); return 1; }
    if (dim!=0 && R1!=R2) { fprintf(stderr,"error in join2_z: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1 && C1!=C2) { fprintf(stderr,"error in join2_z: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2 && S1!=S2) { fprintf(stderr,"error in join2_z: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3 && H1!=H2) { fprintf(stderr,"error in join2_z: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0) ? R1 : (dim==1) ? R1*C1 : (dim==2) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0) ? H1*S1*C1*R1 : (dim==1) ? H1*S1*C1 : (dim==2) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0) ? R2 : (dim==1) ? R2*C2 : (dim==2) ? R2*C2*S2 : R2*C2*S2*H2) : ((dim==0) ? H2*S2*C2*R2 : (dim==1) ? H2*S2*C2 : (dim==2) ? H2*S2 : H2);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=0; g<G; ++g)
    {
        cblas_zcopy((int)B1,X1,1,Y,1); X1 += 2*B1; Y += 2*B1;
        cblas_zcopy((int)B2,X2,1,Y,1); X2 += 2*B2; Y += 2*B2;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
