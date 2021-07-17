//Joins (stacks) 3 inputs X1, X2, X3 into 1 output Y

#include <stdio.h>

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
    if (dim>3u) { fprintf(stderr,"error in join3_s: dim must be in [0 3]\n"); return 1; }
    if (dim!=0u && (R1!=R2 || R1!=R3)) { fprintf(stderr,"error in join3_s: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1u && (C1!=C2 || C1!=C3)) { fprintf(stderr,"error in join3_s: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2u && (S1!=S2 || S1!=S3)) { fprintf(stderr,"error in join3_s: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3u && (H1!=H2 || H1!=H3)) { fprintf(stderr,"error in join3_s: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0u) ? R1 : (dim==1u) ? R1*C1 : (dim==2u) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0u) ? H1*S1*C1*R1 : (dim==1u) ? H1*S1*C1 : (dim==2u) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0u) ? R2 : (dim==1u) ? R2u*C2 : (dim==2u) ? R2u*C2u*S2 : R2u*C2u*S2u*H2) : ((dim==0u) ? H2u*S2u*C2u*R2 : (dim==1u) ? H2u*S2u*C2 : (dim==2u) ? H2u*S2 : H2);
    const size_t B3 = (iscolmajor) ? ((dim==0u) ? R3 : (dim==1u) ? R3u*C3 : (dim==2u) ? R3u*C3u*S3 : R3u*C3u*S3u*H3) : ((dim==0u) ? H3u*S3u*C3u*R3 : (dim==1u) ? H3u*S3u*C3 : (dim==2u) ? H3u*S3 : H3);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=0u; g<G; ++g)
    {
        for (size_t b=0u; b<B1; ++b, ++X1, ++Y) { *Y = *X1; }
        for (size_t b=0u; b<B2; ++b, ++X2, ++Y) { *Y = *X2; }
        for (size_t b=0u; b<B3; ++b, ++X3, ++Y) { *Y = *X3; }
    }

    return 0;
}


int join3_d (double *Y, const double *X1, const double *X2, const double *X3,  const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const size_t R3, const size_t C3, const size_t S3, const size_t H3, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in join3_d: dim must be in [0 3]\n"); return 1; }
    if (dim!=0u && (R1!=R2 || R1!=R3)) { fprintf(stderr,"error in join3_d: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1u && (C1!=C2 || C1!=C3)) { fprintf(stderr,"error in join3_d: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2u && (S1!=S2 || S1!=S3)) { fprintf(stderr,"error in join3_d: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3u && (H1!=H2 || H1!=H3)) { fprintf(stderr,"error in join3_d: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0u) ? R1 : (dim==1u) ? R1*C1 : (dim==2u) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0u) ? H1*S1*C1*R1 : (dim==1u) ? H1*S1*C1 : (dim==2u) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0u) ? R2 : (dim==1u) ? R2u*C2 : (dim==2u) ? R2u*C2u*S2 : R2u*C2u*S2u*H2) : ((dim==0u) ? H2u*S2u*C2u*R2 : (dim==1u) ? H2u*S2u*C2 : (dim==2u) ? H2u*S2 : H2);
    const size_t B3 = (iscolmajor) ? ((dim==0u) ? R3 : (dim==1u) ? R3u*C3 : (dim==2u) ? R3u*C3u*S3 : R3u*C3u*S3u*H3) : ((dim==0u) ? H3u*S3u*C3u*R3 : (dim==1u) ? H3u*S3u*C3 : (dim==2u) ? H3u*S3 : H3);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=0u; g<G; ++g)
    {
        for (size_t b=0u; b<B1; ++b, ++X1, ++Y) { *Y = *X1; }
        for (size_t b=0u; b<B2; ++b, ++X2, ++Y) { *Y = *X2; }
        for (size_t b=0u; b<B3; ++b, ++X3, ++Y) { *Y = *X3; }
    }

    return 0;
}


int join3_c (float *Y, const float *X1, const float *X2, const float *X3,  const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const size_t R3, const size_t C3, const size_t S3, const size_t H3, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in join3_c: dim must be in [0 3]\n"); return 1; }
    if (dim!=0u && (R1!=R2 || R1!=R3)) { fprintf(stderr,"error in join3_c: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1u && (C1!=C2 || C1!=C3)) { fprintf(stderr,"error in join3_c: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2u && (S1!=S2 || S1!=S3)) { fprintf(stderr,"error in join3_c: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3u && (H1!=H2 || H1!=H3)) { fprintf(stderr,"error in join3_c: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0u) ? R1 : (dim==1u) ? R1*C1 : (dim==2u) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0u) ? H1*S1*C1*R1 : (dim==1u) ? H1*S1*C1 : (dim==2u) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0u) ? R2 : (dim==1u) ? R2u*C2 : (dim==2u) ? R2u*C2u*S2 : R2u*C2u*S2u*H2) : ((dim==0u) ? H2u*S2u*C2u*R2 : (dim==1u) ? H2u*S2u*C2 : (dim==2u) ? H2u*S2 : H2);
    const size_t B3 = (iscolmajor) ? ((dim==0u) ? R3 : (dim==1u) ? R3u*C3 : (dim==2u) ? R3u*C3u*S3 : R3u*C3u*S3u*H3) : ((dim==0u) ? H3u*S3u*C3u*R3 : (dim==1u) ? H3u*S3u*C3 : (dim==2u) ? H3u*S3 : H3);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=0u; g<G; ++g)
    {
        for (size_t b=0u; b<B1; ++b, ++X1, ++Y) { *Y = *X1; *++Y = *++X1; }
        for (size_t b=0u; b<B2; ++b, ++X2, ++Y) { *Y = *X2; *++Y = *++X2; }
        for (size_t b=0u; b<B3; ++b, ++X3, ++Y) { *Y = *X3; *++Y = *++X3; }
    }

    return 0;
}


int join3_z (double *Y, const double *X1, const double *X2, const double *X3,  const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const size_t R3, const size_t C3, const size_t S3, const size_t H3, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in join3_z: dim must be in [0 3]\n"); return 1; }
    if (dim!=0u && (R1!=R2 || R1!=R3)) { fprintf(stderr,"error in join3_z: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1u && (C1!=C2 || C1!=C3)) { fprintf(stderr,"error in join3_z: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2u && (S1!=S2 || S1!=S3)) { fprintf(stderr,"error in join3_z: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3u && (H1!=H2 || H1!=H3)) { fprintf(stderr,"error in join3_z: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0u) ? R1 : (dim==1u) ? R1*C1 : (dim==2u) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0u) ? H1*S1*C1*R1 : (dim==1u) ? H1*S1*C1 : (dim==2u) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0u) ? R2 : (dim==1u) ? R2u*C2 : (dim==2u) ? R2u*C2u*S2 : R2u*C2u*S2u*H2) : ((dim==0u) ? H2u*S2u*C2u*R2 : (dim==1u) ? H2u*S2u*C2 : (dim==2u) ? H2u*S2 : H2);
    const size_t B3 = (iscolmajor) ? ((dim==0u) ? R3 : (dim==1u) ? R3u*C3 : (dim==2u) ? R3u*C3u*S3 : R3u*C3u*S3u*H3) : ((dim==0u) ? H3u*S3u*C3u*R3 : (dim==1u) ? H3u*S3u*C3 : (dim==2u) ? H3u*S3 : H3);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=0u; g<G; ++g)
    {
        for (size_t b=0u; b<B1; ++b, ++X1, ++Y) { *Y = *X1; *++Y = *++X1; }
        for (size_t b=0u; b<B2; ++b, ++X2, ++Y) { *Y = *X2; *++Y = *++X2; }
        for (size_t b=0u; b<B3; ++b, ++X3, ++Y) { *Y = *X3; *++Y = *++X3; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
