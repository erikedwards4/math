//Vec2vec operation.
//Prepads each vector in X with P elements equal to val.
//Thus, Y has the same size as X along all other dims, but has a greater length than X along dim.

//For complex case, val is used for both real and imag parts.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int prepad_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const float val);
int prepad_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const double val);
int prepad_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const float val);
int prepad_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const double val);


int prepad_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const float val)
{
    if (dim>3u) { fprintf(stderr,"error in prepad_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u && Lx>0) {}
    else if (Lx==N)
    {
        for (size_t p=0u; p<P; ++p, ++Y) { *Y = val; }
        for (size_t l=0u; l<Lx; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t p=0u; p<P; ++p, ++Y) { *Y = val; }
                for (size_t l=0u; l<Lx; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1u)
        {
            for (size_t n=0u; n<V*P; ++n, ++Y) { *Y = val; }
            for (size_t n=0u; n<V*Lx; ++n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            const size_t Ly = Lx + P;
            for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*Lx-1u, Y-=K*Ly-1u)
                {
                    for (size_t p=0u; p<P; ++p, Y+=K) { *Y = val; }
                    for (size_t l=0u; l<Lx; ++l, X+=K, Y+=K) { *Y = *X; }
                }
            }
        }
    }

    return 0;
}



int prepad_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const double val)
{
    if (dim>3u) { fprintf(stderr,"error in prepad_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u && Lx>0) {}
    else if (Lx==N)
    {
        for (size_t p=0u; p<P; ++p, ++Y) { *Y = val; }
        for (size_t l=0u; l<Lx; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t p=0u; p<P; ++p, ++Y) { *Y = val; }
                for (size_t l=0u; l<Lx; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1u)
        {
            for (size_t n=0u; n<V*P; ++n, ++Y) { *Y = val; }
            for (size_t n=0u; n<V*Lx; ++n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            const size_t Ly = Lx + P;
            for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*Lx-1u, Y-=K*Ly-1u)
                {
                    for (size_t p=0u; p<P; ++p, Y+=K) { *Y = val; }
                    for (size_t l=0u; l<Lx; ++l, X+=K, Y+=K) { *Y = *X; }
                }
            }
        }
    }

    return 0;
}


int prepad_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const float val)
{
    if (dim>3u) { fprintf(stderr,"error in prepad_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u && Lx>0) {}
    else if (Lx==N)
    {
        for (size_t p=0u; p<2u*P; ++p, ++Y) { *Y = val; }
        for (size_t l=0u; l<2u*Lx; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t p=0u; p<2u*P; ++p, ++Y) { *Y = val; }
                for (size_t l=0u; l<2u*Lx; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1u)
        {
            for (size_t n=0u; n<2u*V*P; ++n, ++Y) { *Y = val; }
            for (size_t n=0u; n<2u*V*Lx; ++n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            const size_t Ly = Lx + P;
            for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*Lx-2u, Y-=2u*K*Ly-2u)
                {
                    for (size_t p=0u; p<P; ++p, Y+=2u*K-1u) { *Y = val; *++Y = val; }
                    for (size_t l=0u; l<Lx; ++l, X+=2u*K-1u, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }
                }
            }
        }
    }

    return 0;
}


int prepad_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const double val)
{
    if (dim>3u) { fprintf(stderr,"error in prepad_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u && Lx>0) {}
    else if (Lx==N)
    {
        for (size_t p=0u; p<2u*P; ++p, ++Y) { *Y = val; }
        for (size_t l=0u; l<2u*Lx; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t p=0u; p<2u*P; ++p, ++Y) { *Y = val; }
                for (size_t l=0u; l<2u*Lx; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1u)
        {
            for (size_t n=0u; n<2u*V*P; ++n, ++Y) { *Y = val; }
            for (size_t n=0u; n<2u*V*Lx; ++n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            const size_t Ly = Lx + P;
            for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*Lx-2u, Y-=2u*K*Ly-2u)
                {
                    for (size_t p=0u; p<P; ++p, Y+=2u*K-1u) { *Y = val; *++Y = val; }
                    for (size_t l=0u; l<Lx; ++l, X+=2u*K-1u, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }
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
