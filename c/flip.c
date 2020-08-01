//Flips X along dim (reverses order of elements).
//For d=0, this is flipud. For d=1, this is fliplr.
//This has in-place and not-in-place versions.

//This is a vec2vec operation (see readme.txt).

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int flip_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int flip_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int flip_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int flip_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);

int flip_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int flip_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int flip_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int flip_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int flip_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        Y += L - 1;
        for (size_t l=0; l<L; ++l, ++X, --Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            Y += L - 1;
            for (size_t v=0; v<V; ++v, Y+=2*L)
            {
                for (size_t l=0; l<L; ++l, ++X, --Y) { *Y = *X; }
            }
        }
        else if (G==1)
        {
            Y += V*(L-1);
            for (size_t l=0; l<L; ++l, X+=V, Y-=V) { cblas_scopy((int)V,X,1,Y,1); }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, Y+=K+1)
                {
                    Y += K*(L-1);
                    for (size_t l=0; l<L; ++l, Y-=K, X+=K) { *Y = *X; }
                }
            }
        }
    }

    return 0;
}


int flip_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        Y += L - 1;
        for (size_t l=0; l<L; ++l, ++X, --Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            Y += L - 1;
            for (size_t v=0; v<V; ++v, Y+=2*L)
            {
                for (size_t l=0; l<L; ++l, ++X, --Y) { *Y = *X; }
            }
        }
        else if (G==1)
        {
            Y += V*(L-1);
            for (size_t l=0; l<L; ++l, X+=V, Y-=V) { cblas_dcopy((int)V,X,1,Y,1); }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, Y+=K+1)
                {
                    Y += K*(L-1);
                    for (size_t l=0; l<L; ++l, Y-=K, X+=K) { *Y = *X; }
                }
            }
        }
    }

    return 0;
}


int flip_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1) { cblas_ccopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        Y += 2*(L-1);
        for (size_t l=0; l<L; ++l, Y-=2) { *Y = *X++; *(Y+1) = *X++; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            Y += 2*(L-1);
            for (size_t v=0; v<V; ++v, Y+=4*L)
            {
                for (size_t l=0; l<L; ++l, Y-=2) { *Y = *X++; *(Y+1) = *X++; }
            }
        }
        else if (G==1)
        {
            Y += 2*V*(L-1);
            for (size_t l=0; l<L; ++l, X+=2*V, Y-=2*V) { cblas_ccopy((int)V,X,1,Y,1); }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1), Y+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, Y+=2*K+2)
                {
                    Y += 2*K*(L-1);
                    for (size_t l=0; l<L; ++l, Y-=2*K, X+=2*K) { *Y = *X; *(Y+1) = *(X+1); }
                }
            }
        }
    }

    return 0;
}


int flip_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1) { cblas_zcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        Y += 2*(L-1);
        for (size_t l=0; l<L; ++l, Y-=2) { *Y = *X++; *(Y+1) = *X++; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            Y += 2*(L-1);
            for (size_t v=0; v<V; ++v, Y+=4*L)
            {
                for (size_t l=0; l<L; ++l, Y-=2) { *Y = *X++; *(Y+1) = *X++; }
            }
        }
        else if (G==1)
        {
            Y += 2*V*(L-1);
            for (size_t l=0; l<L; ++l, X+=2*V, Y-=2*V) { cblas_zcopy((int)V,X,1,Y,1); }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1), Y+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, Y+=2*K+2)
                {
                    Y += 2*K*(L-1);
                    for (size_t l=0; l<L; ++l, Y-=2*K, X+=2*K) { *Y = *X; *(Y+1) = *(X+1); }
                }
            }
        }
    }

    return 0;
}


int flip_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float x1;

    if (N==0 || L==1) {}
    else if (L==N)
    {   
        for (size_t l=0; l<L/2; ++l, ++X) { x1 = *X; *X = *(X+L-l-1); *(X+L-l-1) = x1; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L-L/2)
            {
                for (size_t l=0; l<L/2; ++l, ++X) { x1 = *X; *X = *(X+L-l-1); *(X+L-l-1) = x1; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X)
                {
                    for (size_t l=0, l2=K*(L-1); l<K*L/2; l+=K, l2-=K) { x1 = X[l]; X[l] = X[l2]; X[l2] = x1; }
                }
            }
        }
    }

    return 0;
}


int flip_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double x1;

    if (N==0 || L==1) {}
    else if (L==N)
    {   
        for (size_t l=0; l<L/2; ++l, ++X) { x1 = *X; *X = *(X+L-l-1); *(X+L-l-1) = x1; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L-L/2)
            {
                for (size_t l=0; l<L/2; ++l, ++X) { x1 = *X; *X = *(X+L-l-1); *(X+L-l-1) = x1; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X)
                {
                    for (size_t l=0, l2=K*(L-1); l<K*L/2; l+=K, l2-=K) { x1 = X[l]; X[l] = X[l2]; X[l2] = x1; }
                }
            }
        }
    }

    return 0;
}


int flip_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_inplace_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float x1[2];

    if (N==0 || L==1) {}
    else if (L==N)
    {   
        for (size_t l=0; l<L; l+=2, X+=2)
        {
            x1[0] = *X; x1[1] = *(X+1);
            *X = *(X+2*(L-l)-2); *(X+1) = *(X+2*(L-l)-1);
            *(X+2*(L-l)-2) = x1[0]; *(X+2*(L-l)-1) = x1[1];
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=2*(L-L/2))
            {
                for (size_t l=0; l<L; l+=2, X+=2)
                {
                    x1[0] = *X; x1[1] = *(X+1);
                    *X = *(X+2*(L-l)-2); *(X+1) = *(X+2*(L-l)-1);
                    *(X+2*(L-l)-2) = x1[0]; *(X+2*(L-l)-1) = x1[1];
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X+=2)
                {
                    for (size_t l=0, l2=2*K*(L-1); l<K*L; l+=2*K, l2-=2*K)
                    {
                        x1[0] = X[l]; x1[1] = X[l+1];
                        X[l] = X[l2]; X[l+1] = X[l2+1];
                        X[l2] = x1[0]; X[l2+1] = x1[1];
                    }
                }
            }
        }
    }

    return 0;
}


int flip_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in flip_inplace_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double x1[2];

    if (N==0 || L==1) {}
    else if (L==N)
    {   
        for (size_t l=0; l<L; l+=2, X+=2)
        {
            x1[0] = *X; x1[1] = *(X+1);
            *X = *(X+2*(L-l)-2); *(X+1) = *(X+2*(L-l)-1);
            *(X+2*(L-l)-2) = x1[0]; *(X+2*(L-l)-1) = x1[1];
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=2*(L-L/2))
            {
                for (size_t l=0; l<L; l+=2, X+=2)
                {
                    x1[0] = *X; x1[1] = *(X+1);
                    *X = *(X+2*(L-l)-2); *(X+1) = *(X+2*(L-l)-1);
                    *(X+2*(L-l)-2) = x1[0]; *(X+2*(L-l)-1) = x1[1];
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X+=2)
                {
                    for (size_t l=0, l2=2*K*(L-1); l<K*L; l+=2*K, l2-=2*K)
                    {
                        x1[0] = X[l]; x1[1] = X[l+1];
                        X[l] = X[l2]; X[l+1] = X[l2+1];
                        X[l2] = x1[0]; X[l2+1] = x1[1];
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
