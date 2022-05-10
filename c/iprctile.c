//Vec2scalar (reduction) operation.
//Gets index of pth percentile for each vector in X along dim.
//This is the index with value closest to the pth percentile,
//and rounds up for exact ties (0.5 case).

#include <stdio.h>
//#include <stdlib.h>
#include <math.h>
#include "codee_math.h"
//#include "cmpif.c"
//#include "partial_sortif.c"
#include "kselectif.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int iprctile_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p)
{
    if (dim>3u) { fprintf(stderr,"error in iprctile_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>100.0f) { fprintf(stderr,"error in iprctile_s: p must be in [0 100]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    FLT_F *XI;
    if (!(XI=(FLT_F *)malloc(L*sizeof(FLT_F)))) { fprintf(stderr,"error in iprctile_s: problem with malloc. "); perror("malloc"); return 1; }

    //Get index closest to pth prctile after sorting
    const float p1 = (p/100.0f)*(float)(L-1u);
    const size_t i1 = (p<100.0f) ? (size_t)floorf(p1) : L-2u;
    const float w2 = (p<100.0f) ? p1-floorf(p1) : 1.0f;
    const size_t ip = (w2<0.5f) ? i1 : i1+1u;

    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (float)l; }
        //qsort(XI,L,sizeof(FLT_F),cmpif_ascend_s);
        //partial_sortif_s(XI,L,ip,1);
        //*Y = XI[ip].ind;
        *Y = kselectif_s(XI,L-1u,ip,1).ind;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, ++Y)
            {
                for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (float)l; }
                *Y = kselectif_s(XI,L-1u,ip,1).ind;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K) { XI[l].val = *X; XI[l].ind = (float)l; }
                    *Y = kselectif_s(XI,L-1u,ip,1).ind;
                }
            }
        }
    }

    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);

    free(XI);
    return 0;
}


int iprctile_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p)
{
    if (dim>3u) { fprintf(stderr,"error in iprctile_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>100.0) { fprintf(stderr,"error in iprctile_d: p must be in [0 100]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    DBL_D *XI;
    if (!(XI=(DBL_D *)malloc(L*sizeof(DBL_D)))) { fprintf(stderr,"error in iprctile_d: problem with malloc. "); perror("malloc"); return 1; }

    //Get index closest to pth prctile after sorting
    const double p1 = (p/100.0)*(double)(L-1u);
    const size_t i1 = (p<100.0) ? (size_t)floor(p1) : L-2u;
    const double w2 = (p<100.0) ? p1-floor(p1) : 1.0;
    const size_t ip = (w2<0.5) ? i1 : i1+1u;
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (double)l; }
        *Y = kselectif_d(XI,L-1u,ip,1).ind;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, ++Y)
            {
                for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (double)l; }
                *Y = kselectif_d(XI,L-1u,ip,1).ind;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K) { XI[l].val = *X; XI[l].ind = (double)l; }
                    *Y = kselectif_d(XI,L-1u,ip,1).ind;
                }
            }
        }
    }

    free(XI);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
