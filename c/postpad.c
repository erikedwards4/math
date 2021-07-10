//Vec2vec operation.
//Prepads each vector in X with P elements equal to val.
//Thus, Y has the same size as X along all other dims, but has a greater length than X along dim.

//For complex case, val is used for both real and imag parts.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int postpad_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const float val);
int postpad_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const double val);
int postpad_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const float val);
int postpad_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const double val);


int postpad_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const float val)
{
    if (dim>3) { fprintf(stderr,"error in postpad_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0 && Lx>0) {}
    else if (Lx==N)
    {
        for (size_t l=0u; l<Lx; ++l, ++X, ++Y) { *Y = *X; }
        for (size_t p=0u; p<P; ++p, ++Y) { *Y = val; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<Lx; ++l, ++X, ++Y) { *Y = *X; }
                for (size_t p=0u; p<P; ++p, ++Y) { *Y = val; }
            }
        }
        else if (G==1)
        {
            for (size_t n=0u; n<V*Lx; ++n, ++X, ++Y) { *Y = *X; }
            for (size_t n=0u; n<V*P; ++n, ++Y) { *Y = val; }
        }
        else
        {
            const size_t Ly = Lx + P;
            for (size_t g=0u; g<G; ++g, X+=B*(Lx-1), Y+=B*(Ly-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*Lx-1, Y-=K*Ly-1)
                {
                    for (size_t l=0u; l<Lx; ++l, X+=K, Y+=K) { *Y = *X; }
                    for (size_t p=0u; p<P; ++p, Y+=K) { *Y = val; }
                }
            }
        }
    }

    return 0;
}



int postpad_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const double val)
{
    if (dim>3) { fprintf(stderr,"error in postpad_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0 && Lx>0) {}
    else if (Lx==N)
    {
        for (size_t l=0u; l<Lx; ++l, ++X, ++Y) { *Y = *X; }
        for (size_t p=0u; p<P; ++p, ++Y) { *Y = val; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<Lx; ++l, ++X, ++Y) { *Y = *X; }
                for (size_t p=0u; p<P; ++p, ++Y) { *Y = val; }
            }
        }
        else if (G==1)
        {
            for (size_t n=0u; n<V*Lx; ++n, ++X, ++Y) { *Y = *X; }
            for (size_t n=0u; n<V*P; ++n, ++Y) { *Y = val; }
        }
        else
        {
            const size_t Ly = Lx + P;
            for (size_t g=0u; g<G; ++g, X+=B*(Lx-1), Y+=B*(Ly-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*Lx-1, Y-=K*Ly-1)
                {
                    for (size_t l=0u; l<Lx; ++l, X+=K, Y+=K) { *Y = *X; }
                    for (size_t p=0u; p<P; ++p, Y+=K) { *Y = val; }
                }
            }
        }
    }

    return 0;
}


int postpad_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const float val)
{
    if (dim>3) { fprintf(stderr,"error in postpad_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0 && Lx>0) {}
    else if (Lx==N)
    {
        for (size_t l=0u; l<2*Lx; ++l, ++X, ++Y) { *Y = *X; }
        for (size_t p=0u; p<2*P; ++p, ++Y) { *Y = val; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<2*Lx; ++l, ++X, ++Y) { *Y = *X; }
                for (size_t p=0u; p<2*P; ++p, ++Y) { *Y = val; }
            }
        }
        else if (G==1)
        {
            for (size_t n=0u; n<2*V*Lx; ++n, ++X, ++Y) { *Y = *X; }
            for (size_t n=0u; n<2*V*P; ++n, ++Y) { *Y = val; }
        }
        else
        {
            const size_t Ly = Lx + P;
            for (size_t g=0u; g<G; ++g, X+=2*B*(Lx-1), Y+=2*B*(Ly-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*Lx-2, Y-=2*K*Ly-2)
                {
                    for (size_t l=0u; l<Lx; ++l, X+=2*K-1, Y+=2*K-1) { *Y = *X; *++Y = *++X; }
                    for (size_t p=0u; p<P; ++p, Y+=2*K-1) { *Y = val; *++Y = val; }
                }
            }
        }
    }

    return 0;
}


int postpad_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const double val)
{
    if (dim>3) { fprintf(stderr,"error in postpad_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0 && Lx>0) {}
    else if (Lx==N)
    {
        for (size_t l=0u; l<2*Lx; ++l, ++X, ++Y) { *Y = *X; }
        for (size_t p=0u; p<2*P; ++p, ++Y) { *Y = val; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<2*Lx; ++l, ++X, ++Y) { *Y = *X; }
                for (size_t p=0u; p<2*P; ++p, ++Y) { *Y = val; }
            }
        }
        else if (G==1)
        {
            for (size_t n=0u; n<2*V*Lx; ++n, ++X, ++Y) { *Y = *X; }
            for (size_t n=0u; n<2*V*P; ++n, ++Y) { *Y = val; }
        }
        else
        {
            const size_t Ly = Lx + P;
            for (size_t g=0u; g<G; ++g, X+=2*B*(Lx-1), Y+=2*B*(Ly-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*Lx-2, Y-=2*K*Ly-2)
                {
                    for (size_t l=0u; l<Lx; ++l, X+=2*K-1, Y+=2*K-1) { *Y = *X; *++Y = *++X; }
                    for (size_t p=0u; p<P; ++p, Y+=2*K-1) { *Y = val; *++Y = val; }
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
