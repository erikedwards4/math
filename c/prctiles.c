//Gets all percentiles in P along dim of X.
//P is of length Ly, so Y ends up with length Ly along dim.

//The in-place versions still return the same Y, but modify X during processing.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int prctiles_s (float *Y, const float *X, const float *P, const size_t Ly, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int prctiles_d (double *Y, const double *X, const double *P, const size_t Ly, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);

int prctiles_inplace_s (float *Y, float *X, const float *P, const size_t Ly, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int prctiles_inplace_d (double *Y, double *X, const double *P, const size_t Ly, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int prctiles_s (float *Y, const float *X, const float *P, const size_t Ly, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in prctiles_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    float *X1;
    if (!(X1=(float *)malloc(Lx*sizeof(float)))) { fprintf(stderr,"error in prctiles_inplace_s: problem with malloc. "); perror("malloc"); return 1; }

    //Prep interpolation
    size_t *i1;
    float *w1, *w2;
    if (!(i1=(size_t *)malloc(Ly*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w1=(float *)malloc(Ly*sizeof(float)))) { fprintf(stderr,"error in prctiles_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w2=(float *)malloc(Ly*sizeof(float)))) { fprintf(stderr,"error in prctiles_s: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t l=0; l<Ly; l++)
    {
        if (P[l]<0.0f || P[l]>100.0f) { fprintf(stderr,"error in prctiles_s: prctiles must be in [0 100]\n"); return 1; }
        float p1 = (P[l]/100.0f)*(Lx-1);
        i1[l] = (P[l]<100.0f) ? (size_t)floorf(p1) : Lx-2;
        w2[l] = (P[l]<100.0f) ? p1-floorf(p1) : 1.0f;
        w1[l] = 1.0f - w2[l];
    }
    
    if (N==0 || Ly==0) {}
    else if (Lx==1 && Ly==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (Lx==N)
    {
        cblas_scopy((int)Lx,X,1,X1,1);
        if (LAPACKE_slasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_s: problem with LAPACKE function\n"); }
        X1 += *i1++;
        *Y++ = *w1++**X1 + *w2++**(X1+1);
        for (size_t l=1; l<Ly; l++, i1++, w1++, w2++)
        {
            X1 += *i1 - *(i1-1);
            *Y++ = *w1**X + *w2**(X1+1);
        }
        X1 -= *(i1-1); i1 -= Ly; w1 -= Ly; w2 -= Ly;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, X+=Lx, X1-=*(i1-1), i1-=Ly, w1-=Ly, w2-=Ly)
            {
                cblas_scopy((int)Lx,X,1,X1,1);
                if (LAPACKE_slasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_s: problem with LAPACKE function\n"); }
                X1 += *i1++;
                *Y++ = *w1++**X1 + *w2++**(X1+1);
                for (size_t l=1; l<Ly; l++, i1++, w1++, w2++)
                {
                    X1 += *i1 - *(i1-1);
                    *Y++ = *w1**X1 + *w2**(X1+1);
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=B*(Lx-1), Y+=B*(Ly-1))
            {
                for (size_t b=0; b<B; b++, X++, Y-=Ly*K-1, X1-=*(i1-1), i1-=Ly, w1-=Ly, w2-=Ly)
                {
                    cblas_scopy((int)Lx,X,(int)K,X1,1);
                    if (LAPACKE_slasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_inplace_s: problem with LAPACKE function\n"); }
                    X1 += *i1++;
                    *Y = *w1++**X1 + *w2++**(X1+1);
                    Y += K;
                    for (size_t l=1; l<Ly; l++, Y+=K, i1++, w1++, w2++)
                    {
                        X1 += *i1 - *(i1-1);
                        *Y = *w1**X1 + *w2**(X1+1);
                    }
                }
            }
        }
    }

    free(X1);
    return 0;
}


int prctiles_d (double *Y, const double *X, const double *P, const size_t Ly, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in prctiles_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    double *X1;
    if (!(X1=(double *)malloc(Lx*sizeof(double)))) { fprintf(stderr,"error in prctiles_inplace_d: problem with malloc. "); perror("malloc"); return 1; }

    //Prep interpolation
    size_t *i1;
    double *w1, *w2;
    if (!(i1=(size_t *)malloc(Ly*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w1=(double *)malloc(Ly*sizeof(double)))) { fprintf(stderr,"error in prctiles_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w2=(double *)malloc(Ly*sizeof(double)))) { fprintf(stderr,"error in prctiles_d: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t l=0; l<Ly; l++)
    {
        if (P[l]<0.0 || P[l]>100.0) { fprintf(stderr,"error in prctiles_d: prctiles must be in [0 100]\n"); return 1; }
        double p1 = (P[l]/100.0)*(Lx-1);
        i1[l] = (P[l]<100.0) ? (size_t)floor(p1) : Lx-2;
        w2[l] = (P[l]<100.0) ? p1-floor(p1) : 1.0;
        w1[l] = 1.0 - w2[l];
    }
    
    if (N==0 || Ly==0) {}
    else if (Lx==1 && Ly==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (Lx==N)
    {
        cblas_dcopy((int)Lx,X,1,X1,1);
        if (LAPACKE_dlasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_d: problem with LAPACKE function\n"); }
        X1 += *i1++;
        *Y++ = *w1++**X1 + *w2++**(X1+1);
        for (size_t l=1; l<Ly; l++, i1++, w1++, w2++)
        {
            X1 += *i1 - *(i1-1);
            *Y++ = *w1**X + *w2**(X1+1);
        }
        X1 -= *(i1-1); i1 -= Ly; w1 -= Ly; w2 -= Ly;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, X+=Lx, X1-=*(i1-1), i1-=Ly, w1-=Ly, w2-=Ly)
            {
                cblas_dcopy((int)Lx,X,1,X1,1);
                if (LAPACKE_dlasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_d: problem with LAPACKE function\n"); }
                X1 += *i1++;
                *Y++ = *w1++**X1 + *w2++**(X1+1);
                for (size_t l=1; l<Ly; l++, i1++, w1++, w2++)
                {
                    X1 += *i1 - *(i1-1);
                    *Y++ = *w1**X1 + *w2**(X1+1);
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=B*(Lx-1), Y+=B*(Ly-1))
            {
                for (size_t b=0; b<B; b++, X++, Y-=Ly*K-1, X1-=*(i1-1), i1-=Ly, w1-=Ly, w2-=Ly)
                {
                    cblas_dcopy((int)Lx,X,(int)K,X1,1);
                    if (LAPACKE_dlasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_inplace_d: problem with LAPACKE function\n"); }
                    X1 += *i1++;
                    *Y = *w1++**X1 + *w2++**(X1+1);
                    Y += K;
                    for (size_t l=1; l<Ly; l++, Y+=K, i1++, w1++, w2++)
                    {
                        X1 += *i1 - *(i1-1);
                        *Y = *w1**X1 + *w2**(X1+1);
                    }
                }
            }
        }
    }

    free(X1);
    return 0;
}


int prctiles_inplace_s (float *Y, float *X, const float *P, const size_t Ly, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in prctiles_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Prep interpolation
    size_t *i1;
    float *w1, *w2;
    if (!(i1=(size_t *)malloc(Ly*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w1=(float *)malloc(Ly*sizeof(float)))) { fprintf(stderr,"error in prctiles_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w2=(float *)malloc(Ly*sizeof(float)))) { fprintf(stderr,"error in prctiles_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t l=0; l<Ly; l++)
    {
        if (P[l]<0.0f || P[l]>100.0f) { fprintf(stderr,"error in prctiles_inplace_s: prctiles must be in [0 100]\n"); return 1; }
        float p1 = (P[l]/100.0f)*(Lx-1);
        i1[l] = (P[l]<100.0f) ? (size_t)floorf(p1) : Lx-2;
        w2[l] = (P[l]<100.0f) ? p1-floorf(p1) : 1.0f;
        w1[l] = 1.0f - w2[l];
    }
    
    if (N==0 || Ly==0) {}
    else if (Lx==1 && Ly==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (Lx==N)
    {
        if (LAPACKE_slasrt_work('I',(int)Lx,X)) { fprintf(stderr,"error in prctiles_inplace_s: problem with LAPACKE function\n"); }
        X += *i1++;
        *Y++ = *w1++**X + *w2++**(X+1);
        for (size_t l=1; l<Ly; l++, i1++, w1++, w2++)
        {
            X += *i1 - *(i1-1);
            *Y++ = *w1**X + *w2**(X+1);
        }
        i1 -= Ly; w1 -= Ly; w2 -= Ly;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, X+=Lx-*(i1-1), i1-=Ly, w1-=Ly, w2-=Ly)
            {
                if (LAPACKE_slasrt_work('I',(int)Lx,X)) { fprintf(stderr,"error in prctiles_inplace_s: problem with LAPACKE function\n"); }
                X += *i1++;
                *Y++ = *w1++**X + *w2++**(X+1);
                for (size_t l=1; l<Ly; l++, i1++, w1++, w2++)
                {
                    X += *i1 - *(i1-1);
                    *Y++ = *w1**X + *w2**(X+1);
                }
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(Lx*sizeof(float)))) { fprintf(stderr,"error in prctiles_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=B*(Lx-1), Y+=B*(Ly-1))
            {
                for (size_t b=0; b<B; b++, X++, Y-=Ly*K-1, X1-=*(i1-1), i1-=Ly, w1-=Ly, w2-=Ly)
                {
                    cblas_scopy((int)Lx,X,(int)K,X1,1);
                    if (LAPACKE_slasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_inplace_s: problem with LAPACKE function\n"); }
                    X1 += *i1++;
                    *Y = *w1++**X1 + *w2++**(X1+1);
                    Y += K;
                    for (size_t l=1; l<Ly; l++, Y+=K, i1++, w1++, w2++)
                    {
                        X1 += *i1 - *(i1-1);
                        *Y = *w1**X1 + *w2**(X1+1);
                    }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int prctiles_inplace_d (double *Y, double *X, const double *P, const size_t Ly, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in prctiles_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Prep interpolation
    size_t *i1;
    double *w1, *w2;
    if (!(i1=(size_t *)malloc(Ly*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w1=(double *)malloc(Ly*sizeof(double)))) { fprintf(stderr,"error in prctiles_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w2=(double *)malloc(Ly*sizeof(double)))) { fprintf(stderr,"error in prctiles_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t l=0; l<Ly; l++)
    {
        if (P[l]<0.0 || P[l]>100.0) { fprintf(stderr,"error in prctiles_inplace_d: prctiles must be in [0 100]\n"); return 1; }
        double p1 = (P[l]/100.0)*(Lx-1);
        i1[l] = (P[l]<100.0) ? (size_t)floor(p1) : Lx-2;
        w2[l] = (P[l]<100.0) ? p1-floor(p1) : 1.0;
        w1[l] = 1.0 - w2[l];
    }
    
    if (N==0 || Ly==0) {}
    else if (Lx==1 && Ly==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (Lx==N)
    {
        if (LAPACKE_dlasrt_work('I',(int)Lx,X)) { fprintf(stderr,"error in prctiles_inplace_d: problem with LAPACKE function\n"); }
        X += *i1++;
        *Y++ = *w1++**X + *w2++**(X+1);
        for (size_t l=1; l<Ly; l++, i1++, w1++, w2++)
        {
            X += *i1 - *(i1-1);
            *Y++ = *w1**X + *w2**(X+1);
        }
        i1 -= Ly; w1 -= Ly; w2 -= Ly;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, X+=Lx-*(i1-1), i1-=Ly, w1-=Ly, w2-=Ly)
            {
                if (LAPACKE_dlasrt_work('I',(int)Lx,X)) { fprintf(stderr,"error in prctiles_inplace_d: problem with LAPACKE function\n"); }
                X += *i1++;
                *Y++ = *w1++**X + *w2++**(X+1);
                for (size_t l=1; l<Ly; l++, i1++, w1++, w2++)
                {
                    X += *i1 - *(i1-1);
                    *Y++ = *w1**X + *w2**(X+1);
                }
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(Lx*sizeof(double)))) { fprintf(stderr,"error in prctiles_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=B*(Lx-1), Y+=B*(Ly-1))
            {
                for (size_t b=0; b<B; b++, X++, Y-=Ly*K-1, X1-=*(i1-1), i1-=Ly, w1-=Ly, w2-=Ly)
                {
                    cblas_dcopy((int)Lx,X,(int)K,X1,1);
                    if (LAPACKE_dlasrt_work('I',(int)Lx,X1)) { fprintf(stderr,"error in prctiles_inplace_d: problem with LAPACKE function\n"); }
                    X1 += *i1++;
                    *Y = *w1++**X1 + *w2++**(X1+1);
                    Y += K;
                    for (size_t l=1; l<Ly; l++, Y+=K, i1++, w1++, w2++)
                    {
                        X1 += *i1 - *(i1-1);
                        *Y = *w1**X1 + *w2**(X1+1);
                    }
                }
            }
            free(X1);
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
