//Vec2scalar (reduction) operation.
//Gets root-mean-square (RMS) for each vector in X along dim.
//This is the square root of the mean of squares for each vector.
//This is the generalized mean with p=2.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int rms_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int rms_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int rms_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int rms_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int rms_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in rms_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float den = 1.0f/(float)L;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        float sm2 = 0.0f;
        for (size_t l=L; l>0u; --l, ++X) { sm2 += *X * *X; }
        *Y = sqrtf(sm2*den);
        //*Y = cblas_snrm2((int)L,X,1) * den2;  //this is more numerically accurate for large L
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            float sm2;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                sm2 = 0.0f;
                for (size_t l=L; l>0u; --l, ++X) { sm2 += *X * *X; }
                *Y = sqrtf(sm2*den);
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X**X; }
            Y -= V;
            for (size_t l=L; l>1u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += *X**X; }
            }
            for (size_t v=V; v>0u; --v, ++Y) { *Y = sqrtf(*Y*den); }
        }
        else
        {
            float sm2;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    sm2 = 0.0f;
                    for (size_t l=L; l>0u; --l, X+=K) { sm2 += *X * *X; }
                    *Y = sqrtf(sm2*den);
                }
            }
        }
    }

    return 0;
}


int rms_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in rms_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double den = 1.0/(double)L;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        double sm2 = 0.0;
        for (size_t l=L; l>0u; --l, ++X) { sm2 += *X * *X; }
        *Y = sqrt(sm2*den);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            double sm2;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                sm2 = 0.0;
                for (size_t l=L; l>0u; --l, ++X) { sm2 += *X * *X; }
                *Y = sqrt(sm2*den);
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X**X; }
            Y -= V;
            for (size_t l=L; l>1u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += *X**X; }
            }
            for (size_t v=V; v>0u; --v, ++Y) { *Y = sqrt(*Y*den); }
        }
        else
        {
            double sm2;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    sm2 = 0.0;
                    for (size_t l=L; l>0u; --l, X+=K) { sm2 += *X * *X; }
                    *Y = sqrt(sm2*den);
                }
            }
        }
    }

    return 0;
}


int rms_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in rms_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float den = 1.0f/(float)L;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, X+=2, ++Y) { *Y = sqrtf(*X**X + *(X+1)**(X+1)); }
    }
    else if (L==N)
    {
        float sm2 = 0.0f;
        for (size_t l=L; l>0u; --l, X+=2) { sm2 += *X**X + *(X+1)**(X+1); }
        *Y = sqrtf(sm2*den);
        //*Y = cblas_scnrm2((int)L,X,1) * den2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            float sm2;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                sm2 = 0.0f;
                for (size_t l=L; l>0u; --l, X+=2) { sm2 += *X**X + *(X+1)**(X+1); }
                *Y = sqrtf(sm2*den);
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y = *X**X + *(X+1)**(X+1); }
            Y -= V;
            for (size_t l=L; l>1u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y += *X**X + *(X+1)**(X+1); }
            }
            for (size_t v=V; v>0u; --v, ++Y) { *Y = sqrtf(*Y); }
        }
        else
        {
            float sm2;
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, ++Y)
                {
                    sm2 = 0.0f;
                    for (size_t l=L; l>0u; --l, X+=2u*K) { sm2 += *X**X + *(X+1)**(X+1); }
                    *Y = sqrtf(sm2*den);
                }
            }
        }
    }

    return 0;
}


int rms_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in rms_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double den = 1.0/(double)L;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, X+=2, ++Y) { *Y = sqrt(*X**X + *(X+1)**(X+1)); }
    }
    else if (L==N)
    {
        double sm2 = 0.0;
        for (size_t l=L; l>0u; --l, X+=2) { sm2 += *X**X + *(X+1)**(X+1); }
        *Y = sqrt(sm2*den);
        //*Y = cblas_dznrm2((int)L,X,1) * den2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            double sm2;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                sm2 = 0.0;
                for (size_t l=L; l>0u; --l, X+=2) { sm2 += *X**X + *(X+1)**(X+1); }
                *Y = sqrt(sm2*den);
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y = *X**X + *(X+1)**(X+1); }
            Y -= V;
            for (size_t l=L; l>1u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y += *X**X + *(X+1)**(X+1); }
            }
            for (size_t v=V; v>0u; --v, ++Y) { *Y = sqrt(*Y); }
        }
        else
        {
            double sm2;
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, ++Y)
                {
                    sm2 = 0.0;
                    for (size_t l=L; l>0u; --l, X+=2u*K) { sm2 += *X**X + *(X+1)**(X+1); }
                    *Y = sqrt(sm2*den);
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
