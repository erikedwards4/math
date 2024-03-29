//Joins 2 inputs X1, X2 into 1 output Y

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int join2_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in join2_s: dim must be in [0 3]\n"); return 1; }
    if (dim!=0u && R1!=R2) { fprintf(stderr,"error in join2_s: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1u && C1!=C2) { fprintf(stderr,"error in join2_s: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2u && S1!=S2) { fprintf(stderr,"error in join2_s: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3u && H1!=H2) { fprintf(stderr,"error in join2_s: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0u) ? R1 : (dim==1u) ? R1*C1 : (dim==2u) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0u) ? H1*S1*C1*R1 : (dim==1u) ? H1*S1*C1 : (dim==2u) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0u) ? R2 : (dim==1u) ? R2*C2 : (dim==2u) ? R2*C2*S2 : R2*C2*S2*H2) : ((dim==0u) ? H2*S2*C2*R2 : (dim==1u) ? H2*S2*C2 : (dim==2u) ? H2*S2 : H2);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=G; g>0u; --g)
    {
        for (size_t b=B1; b>0u; --b, ++X1, ++Y) { *Y = *X1; }
        for (size_t b=B2; b>0u; --b, ++X2, ++Y) { *Y = *X2; }
    }

    return 0;
}


int join2_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in join2_d: dim must be in [0 3]\n"); return 1; }
    if (dim!=0u && R1!=R2) { fprintf(stderr,"error in join2_d: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1u && C1!=C2) { fprintf(stderr,"error in join2_d: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2u && S1!=S2) { fprintf(stderr,"error in join2_d: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3u && H1!=H2) { fprintf(stderr,"error in join2_d: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0u) ? R1 : (dim==1u) ? R1*C1 : (dim==2u) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0u) ? H1*S1*C1*R1 : (dim==1u) ? H1*S1*C1 : (dim==2u) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0u) ? R2 : (dim==1u) ? R2*C2 : (dim==2u) ? R2*C2*S2 : R2*C2*S2*H2) : ((dim==0u) ? H2*S2*C2*R2 : (dim==1u) ? H2*S2*C2 : (dim==2u) ? H2*S2 : H2);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=G; g>0u; --g)
    {
        for (size_t b=B1; b>0u; --b, ++X1, ++Y) { *Y = *X1; }
        for (size_t b=B2; b>0u; --b, ++X2, ++Y) { *Y = *X2; }
    }

    return 0;
}


int join2_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in join2_c: dim must be in [0 3]\n"); return 1; }
    if (dim!=0u && R1!=R2) { fprintf(stderr,"error in join2_c: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1u && C1!=C2) { fprintf(stderr,"error in join2_c: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2u && S1!=S2) { fprintf(stderr,"error in join2_c: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3u && H1!=H2) { fprintf(stderr,"error in join2_c: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0u) ? R1 : (dim==1u) ? R1*C1 : (dim==2u) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0u) ? H1*S1*C1*R1 : (dim==1u) ? H1*S1*C1 : (dim==2u) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0u) ? R2 : (dim==1u) ? R2*C2 : (dim==2u) ? R2*C2*S2 : R2*C2*S2*H2) : ((dim==0u) ? H2*S2*C2*R2 : (dim==1u) ? H2*S2*C2 : (dim==2u) ? H2*S2 : H2);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=G; g>0u; --g)
    {
        for (size_t b=B1; b>0u; --b, ++X1, ++Y) { *Y = *X1; *++Y = *++X1; }
        for (size_t b=B2; b>0u; --b, ++X2, ++Y) { *Y = *X2; *++Y = *++X2; }
    }

    return 0;
}


int join2_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in join2_z: dim must be in [0 3]\n"); return 1; }
    if (dim!=0u && R1!=R2) { fprintf(stderr,"error in join2_z: inputs must have same num rows for dim!=0\n"); return 1; }
    if (dim!=1u && C1!=C2) { fprintf(stderr,"error in join2_z: inputs must have same num cols for dim!=1\n"); return 1; }
    if (dim!=2u && S1!=S2) { fprintf(stderr,"error in join2_z: inputs must have same num slices for dim!=2\n"); return 1; }
    if (dim!=3u && H1!=H2) { fprintf(stderr,"error in join2_z: inputs must have same num hyperslices for dim!=3\n"); return 1; }

    const size_t B1 = (iscolmajor) ? ((dim==0u) ? R1 : (dim==1u) ? R1*C1 : (dim==2u) ? R1*C1*S1 : R1*C1*S1*H1) : ((dim==0u) ? H1*S1*C1*R1 : (dim==1u) ? H1*S1*C1 : (dim==2u) ? H1*S1 : H1);
    const size_t B2 = (iscolmajor) ? ((dim==0u) ? R2 : (dim==1u) ? R2*C2 : (dim==2u) ? R2*C2*S2 : R2*C2*S2*H2) : ((dim==0u) ? H2*S2*C2*R2 : (dim==1u) ? H2*S2*C2 : (dim==2u) ? H2*S2 : H2);
    const size_t G = R1*C1*S1*H1/B1;

    for (size_t g=G; g>0u; --g)
    {
        for (size_t b=B1; b>0u; --b, ++X1, ++Y) { *Y = *X1; *++Y = *++X1; }
        for (size_t b=B2; b>0u; --b, ++X2, ++Y) { *Y = *X2; *++Y = *++X2; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
