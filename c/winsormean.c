//Vec2scalar (reduction) operation.
//Gets winsorized mean for each vector in X along dim.

//This replaces values outside of the pth and (1-q)th percentiles with
//the min and max remaining values, respectively, and takes mean.

//The inplace version still outputs Y, but modifies X during processing.

#include <stdio.h>
#include <math.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int winsormean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q);
int winsormean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q);

int winsormean_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q);
int winsormean_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q);


int winsormean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q)
{
    if (dim>3) { fprintf(stderr,"error in winsormean_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>50.0f) { fprintf(stderr,"error in winsormean_s: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f / L;

    const float p1 = (p/100.0f)*(L-1), p2 = (1.0f-q/100.0f)*(L-1);
    const size_t i1 = (size_t)ceilf(p1), i2 = (size_t)floorf(p2);
    float mn, mx, sm;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in winsormean_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsormean_s: problem with LAPACKE function\n"); }
        X1 += i2; mx = *X1;
        X1 -= i2-i1; mn = *X1++;
        sm = (i1+1)*mn + (L-i2)*mx;
        for (size_t l=i1+1; l<i2; ++l, ++X1) { sm += *X1; }
        X1 -= i2-i1+1;
        *Y = sm * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, X1-=i2-i1+1, ++Y)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsormean_s: problem with LAPACKE function\n"); }
                X1 += i2; mx = *X1;
                X1 -= i2-i1; mn = *X1++;
                sm = (i1+1)*mn + (L-i2)*mx;
                for (size_t l=i1+1; l<i2; ++l, ++X1) { sm += *X1; }
                *Y = sm * den;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, X1-=i2-i1+1, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsormean_s: problem with LAPACKE function\n"); }
                    X1 += i2; mx = *X1;
                    X1 -= i2-i1; mn = *X1++;
                    sm = (i1+1)*mn + (L-i2)*mx;
                    for (size_t l=i1+1; l<i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm * den;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int winsormean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q)
{
    if (dim>3) { fprintf(stderr,"error in winsormean_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in winsormean_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0 / L;

    const double p1 = (p/100.0)*(L-1), p2 = (1.0-q/100.0)*(L-1);
    const size_t i1 = (size_t)ceil(p1), i2 = (size_t)floor(p2);
    double mn, mx, sm;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in winsormean_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsormean_d: problem with LAPACKE function\n"); }
        X1 += i2; mx = *X1;
        X1 -= i2-i1; mn = *X1++;
        sm = (i1+1)*mn + (L-i2)*mx;
        for (size_t l=i1+1; l<i2; ++l, ++X1) { sm += *X1; }
        X1 -= i2-i1+1;
        *Y = sm * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, X1-=i2-i1+1, ++Y)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsormean_d: problem with LAPACKE function\n"); }
                X1 += i2; mx = *X1;
                X1 -= i2-i1; mn = *X1++;
                sm = (i1+1)*mn + (L-i2)*mx;
                for (size_t l=i1+1; l<i2; ++l, ++X1) { sm += *X1; }
                *Y = sm * den;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, X1-=i2-i1+1, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsormean_d: problem with LAPACKE function\n"); }
                    X1 += i2; mx = *X1;
                    X1 -= i2-i1; mn = *X1++;
                    sm = (i1+1)*mn + (L-i2)*mx;
                    for (size_t l=i1+1; l<i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm * den;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int winsormean_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q)
{
    if (dim>3) { fprintf(stderr,"error in winsormean_inplace_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>50.0f) { fprintf(stderr,"error in winsormean_inplace_s: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f / L;

    const float p1 = (p/100.0f)*(L-1), p2 = (1.0f-q/100.0f)*(L-1);
    const size_t i1 = (size_t)ceilf(p1), i2 = (size_t)floorf(p2);
    float mn, mx, sm;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in winsormean_inplace_s: problem with LAPACKE function\n"); }
        X += i2; mx = *X;
        X -= i2-i1; mn = *X++;
        sm = (i1+1)*mn + (L-i2)*mx;
        for (size_t l=i1+1; l<i2; ++l, ++X) { sm += *X; }
        *Y = sm * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, X+=L-i2, ++Y)
            {
                if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in winsormean_inplace_s: problem with LAPACKE function\n"); }
                X += i2; mx = *X;
                X -= i2-i1; mn = *X++;
                sm = (i1+1)*mn + (L-i2)*mx;
                for (size_t l=i1+1; l<i2; ++l, ++X) { sm += *X; }
                *Y = sm * den;
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in winsormean_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, X1-=i2, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsormean_inplace_s: problem with LAPACKE function\n"); }
                    X1 += i2; mx = *X1;
                    X1 -= i2-i1; mn = *X1++;
                    sm = (i1+1)*mn + (L-i2)*mx;
                    for (size_t l=i1+1; l<i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm * den;
                }
            }
            free(X1);
        }
    }

    return 0;
}


int winsormean_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q)
{
    if (dim>3) { fprintf(stderr,"error in winsormean_inplace_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in winsormean_inplace_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0 / L;

    const double p1 = (p/100.0)*(L-1), p2 = (1.0-q/100.0)*(L-1);
    const size_t i1 = (size_t)ceil(p1), i2 = (size_t)floor(p2);
    double mn, mx, sm;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in winsormean_inplace_d: problem with LAPACKE function\n"); }
        X += i2; mx = *X;
        X -= i2-i1; mn = *X++;
        sm = (i1+1)*mn + (L-i2)*mx;
        for (size_t l=i1+1; l<i2; ++l, ++X) { sm += *X; }
        *Y = sm * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, X+=L-i2, ++Y)
            {
                if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in winsormean_inplace_d: problem with LAPACKE function\n"); }
                X += i2; mx = *X;
                X -= i2-i1; mn = *X++;
                sm = (i1+1)*mn + (L-i2)*mx;
                for (size_t l=i1+1; l<i2; ++l, ++X) { sm += *X; }
                *Y = sm * den;
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in winsormean_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, X1-=i2, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsormean_inplace_d: problem with LAPACKE function\n"); }
                    X1 += i2; mx = *X1;
                    X1 -= i2-i1; mn = *X1++;
                    sm = (i1+1)*mn + (L-i2)*mx;
                    for (size_t l=i1+1; l<i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm * den;
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
