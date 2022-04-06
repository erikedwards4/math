//Vec2vec operation.
//Gets all percentiles in P along dim of X.
//P is of length Ly, so Y ends up with length Ly along dim.

//The in-place versions still return the same Y, but modify X during processing.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int prctiles_s (float *Y, const float *X, const float *P, const size_t Ly, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in prctiles_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    float *X1;
    if (!(X1=(float *)malloc(Lx*sizeof(float)))) { fprintf(stderr,"error in prctiles_inplace_s: problem with malloc. "); perror("malloc"); return 1; }

    //Prep interpolation
    size_t *i1;
    float *w1, *w2;
    if (!(i1=(size_t *)malloc(Ly*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w1=(float *)malloc(Ly*sizeof(float)))) { fprintf(stderr,"error in prctiles_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w2=(float *)malloc(Ly*sizeof(float)))) { fprintf(stderr,"error in prctiles_s: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t l=Ly; l>0u; --l, ++P, ++i1, ++w1, ++w2)
    {
        if (*P<0.0f || *P>100.0f) { fprintf(stderr,"error in prctiles_s: prctiles must be in [0 100]\n"); return 1; }
        float p1 = (*P/100.0f)*(float)(Lx-1u);
        *i1 = (*P<100.0f) ? (size_t)floorf(p1) : Lx-2u;
        *w2 = (*P<100.0f) ? p1-floorf(p1) : 1.0f;
        *w1 = 1.0f - *w2;
    }
    i1 -= Ly; w1 -= Ly; w2 -= Ly;
    
    if (N==0u || Ly==0u) {}
    else if (Lx==1u && Ly==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (Lx==N)
    {
        for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= Lx;
        if (LAPACKE_slasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_s: problem with LAPACKE function\n"); }
        X1 += *i1++;
        *Y++ = *w1++**X1 + *w2++**(X1+1);
        for (size_t l=Ly; l>1u; --l, ++i1, ++w1, ++w2, ++Y)
        {
            X1 += *i1 - *(i1-1u);
            *Y = *w1**X + *w2**(X1+1);
        }
        X1 -= *(i1-1u); i1 -= Ly; w1 -= Ly; w2 -= Ly;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X1-=*(i1-1u), i1-=Ly, w1-=Ly, w2-=Ly)
            {
                for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= Lx;
                if (LAPACKE_slasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_s: problem with LAPACKE function\n"); }
                X1 += *i1++;
                *Y++ = *w1++**X1 + *w2++**(X1+1);
                for (size_t l=Ly; l>1u; --l, ++i1, ++w1, ++w2, ++Y)
                {
                    X1 += *i1 - *(i1-1u);
                    *Y = *w1**X1 + *w2**(X1+1);
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*Lx-1u, X1-=*(i1-1u), i1-=Ly, w1-=Ly, w2-=Ly, Y-=K*Ly-1u)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    if (LAPACKE_slasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_inplace_s: problem with LAPACKE function\n"); }
                    X1 += *i1++;
                    *Y = *w1++**X1 + *w2++**(X1+1);
                    Y += K;
                    for (size_t l=Ly; l>1u; --l, ++i1, ++w1, ++w2, Y+=K)
                    {
                        X1 += *i1 - *(i1-1u);
                        *Y = *w1**X1 + *w2**(X1+1);
                    }
                }
            }
        }
    }

    free(X1); free(i1); free(w1); free(w2);
    return 0;
}


int prctiles_d (double *Y, const double *X, const double *P, const size_t Ly, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in prctiles_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    double *X1;
    if (!(X1=(double *)malloc(Lx*sizeof(double)))) { fprintf(stderr,"error in prctiles_inplace_d: problem with malloc. "); perror("malloc"); return 1; }

    //Prep interpolation
    size_t *i1;
    double *w1, *w2;
    if (!(i1=(size_t *)malloc(Ly*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w1=(double *)malloc(Ly*sizeof(double)))) { fprintf(stderr,"error in prctiles_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w2=(double *)malloc(Ly*sizeof(double)))) { fprintf(stderr,"error in prctiles_d: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t l=Ly; l>0u; --l, ++P, ++i1, ++w1, ++w2)
    {
        if (*P<0.0 || *P>100.0) { fprintf(stderr,"error in prctiles_d: prctiles must be in [0 100]\n"); return 1; }
        double p1 = (*P/100.0)*(double)(Lx-1u);
        *i1 = (*P<100.0) ? (size_t)floor(p1) : Lx-2u;
        *w2 = (*P<100.0) ? p1-floor(p1) : 1.0;
        *w1 = 1.0 - *w2;
    }
    i1 -= Ly; w1 -= Ly; w2 -= Ly;
    
    if (N==0u || Ly==0u) {}
    else if (Lx==1u && Ly==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (Lx==N)
    {
        for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= Lx;
        if (LAPACKE_dlasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_d: problem with LAPACKE function\n"); }
        X1 += *i1++;
        *Y++ = *w1++**X1 + *w2++**(X1+1);
        for (size_t l=Ly; l>1u; --l, ++i1, ++w1, ++w2, ++Y)
        {
            X1 += *i1 - *(i1-1u);
            *Y = *w1**X + *w2**(X1+1);
        }
        X1 -= *(i1-1u); i1 -= Ly; w1 -= Ly; w2 -= Ly;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X1-=*(i1-1u), i1-=Ly, w1-=Ly, w2-=Ly)
            {
                for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= Lx;
                if (LAPACKE_dlasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_d: problem with LAPACKE function\n"); }
                X1 += *i1++;
                *Y++ = *w1++**X1 + *w2++**(X1+1);
                for (size_t l=Ly; l>1u; --l, ++i1, ++w1, ++w2, ++Y)
                {
                    X1 += *i1 - *(i1-1u);
                    *Y = *w1**X1 + *w2**(X1+1);
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*Lx-1u, X1-=*(i1-1u), i1-=Ly, w1-=Ly, w2-=Ly, Y-=K*Ly-1u)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    if (LAPACKE_dlasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_inplace_d: problem with LAPACKE function\n"); }
                    X1 += *i1++;
                    *Y = *w1++**X1 + *w2++**(X1+1);
                    Y += K;
                    for (size_t l=Ly; l>1u; --l, ++i1, ++w1, ++w2, Y+=K)
                    {
                        X1 += *i1 - *(i1-1u);
                        *Y = *w1**X1 + *w2**(X1+1);
                    }
                }
            }
        }
    }

    free(X1); free(i1); free(w1); free(w2);
    return 0;
}


int prctiles_inplace_s (float *Y, float *X, const float *P, const size_t Ly, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in prctiles_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Prep interpolation
    size_t *i1;
    float *w1, *w2;
    if (!(i1=(size_t *)malloc(Ly*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w1=(float *)malloc(Ly*sizeof(float)))) { fprintf(stderr,"error in prctiles_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w2=(float *)malloc(Ly*sizeof(float)))) { fprintf(stderr,"error in prctiles_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t l=Ly; l>0u; --l, ++P, ++i1, ++w1, ++w2)
    {
        if (*P<0.0f || *P>100.0f) { fprintf(stderr,"error in prctiles_inplace_s: prctiles must be in [0 100]\n"); return 1; }
        float p1 = (*P/100.0f)*(float)(Lx-1u);
        *i1 = (*P<100.0f) ? (size_t)floorf(p1) : Lx-2u;
        *w2 = (*P<100.0f) ? p1-floorf(p1) : 1.0f;
        *w1 = 1.0f - *w2;
    }
    i1 -= Ly; w1 -= Ly; w2 -= Ly;
    
    if (N==0u || Ly==0u) {}
    else if (Lx==1u && Ly==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (Lx==N)
    {
        if (LAPACKE_slasrt_work('I',(int)Lx,X)) { fprintf(stderr,"error in prctiles_inplace_s: problem with LAPACKE function\n"); }
        X += *i1++;
        *Y++ = *w1++**X + *w2++**(X+1);
        for (size_t l=Ly; l>1u; --l, ++i1, ++w1, ++w2, ++Y)
        {
            X += *i1 - *(i1-1u);
            *Y = *w1**X + *w2**(X+1);
        }
        i1 -= Ly; w1 -= Ly; w2 -= Ly;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=Lx-*(i1-1u), i1-=Ly, w1-=Ly, w2-=Ly)
            {
                if (LAPACKE_slasrt_work('I',(int)Lx,X)) { fprintf(stderr,"error in prctiles_inplace_s: problem with LAPACKE function\n"); }
                X += *i1++;
                *Y++ = *w1++**X + *w2++**(X+1);
                for (size_t l=Ly; l>1u; --l, ++i1, ++w1, ++w2, ++Y)
                {
                    X += *i1 - *(i1-1u);
                    *Y = *w1**X + *w2**(X+1);
                }
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(Lx*sizeof(float)))) { fprintf(stderr,"error in prctiles_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*Lx-1u, X1-=*(i1-1u), i1-=Ly, w1-=Ly, w2-=Ly, Y-=K*Ly-1u)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    if (LAPACKE_slasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_inplace_s: problem with LAPACKE function\n"); }
                    X1 += *i1++;
                    *Y = *w1++**X1 + *w2++**(X1+1);
                    Y += K;
                    for (size_t l=Ly; l>1u; --l, ++i1, ++w1, ++w2, Y+=K)
                    {
                        X1 += *i1 - *(i1-1u);
                        *Y = *w1**X1 + *w2**(X1+1);
                    }
                }
            }
            free(X1);
        }
    }

    free(i1); free(w1); free(w2);
    return 0;
}


int prctiles_inplace_d (double *Y, double *X, const double *P, const size_t Ly, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in prctiles_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Prep interpolation
    size_t *i1;
    double *w1, *w2;
    if (!(i1=(size_t *)malloc(Ly*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w1=(double *)malloc(Ly*sizeof(double)))) { fprintf(stderr,"error in prctiles_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w2=(double *)malloc(Ly*sizeof(double)))) { fprintf(stderr,"error in prctiles_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t l=Ly; l>0u; --l, ++P, ++i1, ++w1, ++w2)
    {
        if (*P<0.0 || *P>100.0) { fprintf(stderr,"error in prctiles_inplace_d: prctiles must be in [0 100]\n"); return 1; }
        double p1 = (*P/100.0)*(double)(Lx-1u);
        *i1 = (*P<100.0) ? (size_t)floor(p1) : Lx-2u;
        *w2 = (*P<100.0) ? p1-floor(p1) : 1.0;
        *w1 = 1.0 - *w2;
    }
    i1 -= Ly; w1 -= Ly; w2 -= Ly;
    
    if (N==0u || Ly==0u) {}
    else if (Lx==1u && Ly==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (Lx==N)
    {
        if (LAPACKE_dlasrt_work('I',(int)Lx,X)) { fprintf(stderr,"error in prctiles_inplace_d: problem with LAPACKE function\n"); }
        X += *i1++;
        *Y++ = *w1++**X + *w2++**(X+1);
        for (size_t l=Ly; l>1u; --l, ++i1, ++w1, ++w2, ++Y)
        {
            X += *i1 - *(i1-1u);
            *Y = *w1**X + *w2**(X+1);
        }
        i1 -= Ly; w1 -= Ly; w2 -= Ly;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=Lx-*(i1-1u), i1-=Ly, w1-=Ly, w2-=Ly)
            {
                if (LAPACKE_dlasrt_work('I',(int)Lx,X)) { fprintf(stderr,"error in prctiles_inplace_d: problem with LAPACKE function\n"); }
                X += *i1++;
                *Y++ = *w1++**X + *w2++**(X+1);
                for (size_t l=Ly; l>1u; --l, ++i1, ++w1, ++w2, ++Y)
                {
                    X += *i1 - *(i1-1u);
                    *Y = *w1**X + *w2**(X+1);
                }
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(Lx*sizeof(double)))) { fprintf(stderr,"error in prctiles_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*Lx-1u, X1-=*(i1-1u), i1-=Ly, w1-=Ly, w2-=Ly, Y-=K*Ly-1u)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    if (LAPACKE_dlasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_inplace_d: problem with LAPACKE function\n"); }
                    X1 += *i1++;
                    *Y = *w1++**X1 + *w2++**(X1+1);
                    Y += K;
                    for (size_t l=Ly; l>1u; --l, ++i1, ++w1, ++w2, Y+=K)
                    {
                        X1 += *i1 - *(i1-1u);
                        *Y = *w1**X1 + *w2**(X1+1);
                    }
                }
            }
            free(X1);
        }
    }

    free(i1); free(w1); free(w2);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
