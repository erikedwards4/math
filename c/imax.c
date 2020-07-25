//Vec2scalar (reduction) operation.
//Gets index of maximum of values for each vector in X along dim.
//For complex case, this is the proper absolute value; see iamax for the other one.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int imax_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int imax_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int imax_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int imax_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int imax_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in imax_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float mx;

    if (N==0) {}
    else if (L==1)
    {
        float z = 0.0f;
        cblas_scopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        mx = *X; *Y = 0.0f;
        for (size_t l=1; l<L; l++) { if (X[l]>mx) { mx = X[l]; *Y = l; } }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, Y++)
            {
                mx = *X++; *Y = 0.0f;
                for (size_t l=1; l<L; l++, X++) { if (*X>mx) { mx = *X; *Y = l; } }
            }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, X-=K*L-1, Y++)
                {
                    mx = *X; *Y = 0.0f; X += K;
                    for (size_t l=1; l<L; l++, X+=K) { if (*X>mx) { mx = *X; *Y = l; } }
                }
            }
        }
    }

    return 0;
}


int imax_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in imax_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double mx;

    if (N==0) {}
    else if (L==1)
    {
        double z = 0.0;
        cblas_dcopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        mx = *X; *Y = 0.0;
        for (size_t l=1; l<L; l++) { if (X[l]>mx) { mx = X[l]; *Y = l; } }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, Y++)
            {
                mx = *X++; *Y = 0.0;
                for (size_t l=1; l<L; l++, X++) { if (*X>mx) { mx = *X; *Y = l; } }
            }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, X-=K*L-1, Y++)
                {
                    mx = *X; *Y = 0.0; X += K;
                    for (size_t l=1; l<L; l++, X+=K) { if (*X>mx) { mx = *X; *Y = l; } }
                }
            }
        }
    }

    return 0;
}


int imax_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in imax_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float xx, mx;

    if (N==0) {}
    else if (L==1)
    {
        float z = 0.0f;
        cblas_scopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        mx = *X**X + *(X+1)**(X+1);
        *Y = 0.0f; X += 2;
        for (size_t l=1; l<L; l++, X+=2)
        {
            xx = *X**X + *(X+1)**(X+1);
            if (xx>mx) { mx = xx; *Y = l; }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, Y++)
            {
                mx = *X**X + *(X+1)**(X+1);
                *Y = 0.0f; X += 2;
                for (size_t l=1; l<L; l++, X+=2)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (xx>mx) { mx = xx; *Y = l; }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; b++, X-=2*K*L-2, Y++)
                {
                    mx = *X**X + *(X+1)**(X+1);
                    *Y = 0.0f; X += 2*K;
                    for (size_t l=1; l<L; l++, X+=2*K)
                    {
                        xx = *X**X + *(X+1)**(X+1);
                        if (xx>mx) { mx = xx; *Y = l; }
                    }
                }
            }
        }
    }

    return 0;
}


int imax_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in imax_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double xx, mx;

    if (N==0) {}
    else if (L==1)
    {
        double z = 0.0;
        cblas_dcopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        mx = *X**X + *(X+1)**(X+1);
        *Y = 0.0; X += 2;
        for (size_t l=1; l<L; l++, X+=2)
        {
            xx = *X**X + *(X+1)**(X+1);
            if (xx>mx) { mx = xx; *Y = l; }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, Y++)
            {
                mx = *X**X + *(X+1)**(X+1);
                *Y = 0.0; X += 2;
                for (size_t l=1; l<L; l++, X+=2)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (xx>mx) { mx = xx; *Y = l; }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; b++, X-=2*K*L-2, Y++)
                {
                    mx = *X**X + *(X+1)**(X+1);
                    *Y = 0.0; X += 2*K;
                    for (size_t l=1; l<L; l++, X+=2*K)
                    {
                        xx = *X**X + *(X+1)**(X+1);
                        if (xx>mx) { mx = xx; *Y = l; }
                    }
                }
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
