//Takes input tensor X and repeats it Nr x Nc x Ns x Nh times,
//so that the output Y has size R*Nr x C*Nc x S*Ns x H*Nh.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int repmat_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const int iscolmajor);
int repmat_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const int iscolmajor);
int repmat_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const int iscolmajor);
int repmat_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const int iscolmajor);


int repmat_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const int iscolmajor)
{
    const size_t N = R*C*S*H, NN = Nr*Nc*Ns*Nh;

    if (N==0u || NN==0u) {}
    else if (NN==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (iscolmajor)
    {
        const size_t RC = R*C, RCS = RC*S;
        for (size_t nh=0u; nh<Nh; ++nh)
        {
            for (size_t h=0u; h<H; ++h)
            {
                for (size_t ns=0u; ns<Ns; ++ns)
                {
                    for (size_t s=0u; s<S; ++s)
                    {
                        for (size_t nc=0u; nc<Nc; ++nc)
                        {
                            for (size_t c=0u; c<C; ++c)
                            {
                                for (size_t nr=0u; nr<Nr; ++nr)
                                {
                                    for (size_t r=0u; r<R; ++r, ++X, ++Y)
                                    {
                                        *Y = *X;
                                    }
                                    if (nr<Nr-1u) { X -= R; }
                                }
                            }
                            if (nc<Nc-1u) { X -= RC; }
                        }
                    }
                    if (ns<Ns-1u) { X -= RCS; }
                }
            }
            if (nh<Nh-1u) { X -= N; }
        }
    }
    else
    {
        const size_t SH = S*H, CSH = C*SH;
        for (size_t nr=0u; nr<Nr; ++nr)
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t nc=0u; nc<Nc; ++nc)
                {
                    for (size_t c=0u; c<C; ++c)
                    {
                        for (size_t ns=0u; ns<Ns; ++ns)
                        {
                            for (size_t s=0u; s<S; ++s)
                            {
                                for (size_t nh=0u; nh<Nh; ++nh)
                                {
                                    for (size_t h=0u; h<H; ++h, ++X, ++Y)
                                    {
                                        *Y = *X;
                                    }
                                    if (nh<Nh-1u) { X -= H; }
                                }
                            }
                            if (ns<Ns-1u) { X -= SH; }
                        }
                    }
                    if (nc<Nc-1u) { X -= CSH; }
                }
            }
            if (nr<Nr-1u) { X -= N; }
        }
    }

    return 0;
}


int repmat_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const int iscolmajor)
{
    const size_t N = R*C*S*H, NN = Nr*Nc*Ns*Nh;

    if (N==0u || NN==0u) {}
    else if (NN==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (iscolmajor)
    {
        const size_t RC = R*C, RCS = RC*S;
        for (size_t nh=0u; nh<Nh; ++nh)
        {
            for (size_t h=0u; h<H; ++h)
            {
                for (size_t ns=0u; ns<Ns; ++ns)
                {
                    for (size_t s=0u; s<S; ++s)
                    {
                        for (size_t nc=0u; nc<Nc; ++nc)
                        {
                            for (size_t c=0u; c<C; ++c)
                            {
                                for (size_t nr=0u; nr<Nr; ++nr)
                                {
                                    for (size_t r=0u; r<R; ++r, ++X, ++Y)
                                    {
                                        *Y = *X;
                                    }
                                    if (nr<Nr-1u) { X -= R; }
                                }
                            }
                            if (nc<Nc-1u) { X -= RC; }
                        }
                    }
                    if (ns<Ns-1u) { X -= RCS; }
                }
            }
            if (nh<Nh-1u) { X -= N; }
        }
    }
    else
    {
        const size_t SH = S*H, CSH = C*SH;
        for (size_t nr=0u; nr<Nr; ++nr)
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t nc=0u; nc<Nc; ++nc)
                {
                    for (size_t c=0u; c<C; ++c)
                    {
                        for (size_t ns=0u; ns<Ns; ++ns)
                        {
                            for (size_t s=0u; s<S; ++s)
                            {
                                for (size_t nh=0u; nh<Nh; ++nh)
                                {
                                    for (size_t h=0u; h<H; ++h, ++X, ++Y)
                                    {
                                        *Y = *X;
                                    }
                                    if (nh<Nh-1u) { X -= H; }
                                }
                            }
                            if (ns<Ns-1u) { X -= SH; }
                        }
                    }
                    if (nc<Nc-1u) { X -= CSH; }
                }
            }
            if (nr<Nr-1u) { X -= N; }
        }
    }

    return 0;
}


int repmat_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const int iscolmajor)
{
    const size_t N = R*C*S*H, NN = Nr*Nc*Ns*Nh;

    if (N==0u || NN==0u) {}
    else if (NN==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (iscolmajor)
    {
        const size_t RC = R*C, RCS = RC*S;
        for (size_t nh=0u; nh<Nh; ++nh)
        {
            for (size_t h=0u; h<H; ++h)
            {
                for (size_t ns=0u; ns<Ns; ++ns)
                {
                    for (size_t s=0u; s<S; ++s)
                    {
                        for (size_t nc=0u; nc<Nc; ++nc)
                        {
                            for (size_t c=0u; c<C; ++c)
                            {
                                for (size_t nr=0u; nr<Nr; ++nr)
                                {
                                    for (size_t r=0u; r<R; ++r, ++X, ++Y)
                                    {
                                        *Y = *X; *++Y = *++X;
                                    }
                                    if (nr<Nr-1u) { X -= 2u*R; }
                                }
                            }
                            if (nc<Nc-1u) { X -= 2u*RC; }
                        }
                    }
                    if (ns<Ns-1u) { X -= 2u*RCS; }
                }
            }
            if (nh<Nh-1u) { X -= 2u*N; }
        }
    }
    else
    {
        const size_t SH = S*H, CSH = C*SH;
        for (size_t nr=0u; nr<Nr; ++nr)
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t nc=0u; nc<Nc; ++nc)
                {
                    for (size_t c=0u; c<C; ++c)
                    {
                        for (size_t ns=0u; ns<Ns; ++ns)
                        {
                            for (size_t s=0u; s<S; ++s)
                            {
                                for (size_t nh=0u; nh<Nh; ++nh)
                                {
                                    for (size_t h=0u; h<H; ++h, ++X, ++Y)
                                    {
                                        *Y = *X; *++Y = *++X;
                                    }
                                    if (nh<Nh-1u) { X -= 2u*H; }
                                }
                            }
                            if (ns<Ns-1u) { X -= 2u*SH; }
                        }
                    }
                    if (nc<Nc-1u) { X -= 2u*CSH; }
                }
            }
            if (nr<Nr-1u) { X -= 2u*N; }
        }
    }

    return 0;
}


int repmat_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const int iscolmajor)
{
    const size_t N = R*C*S*H, NN = Nr*Nc*Ns*Nh;

    if (N==0u || NN==0u) {}
    else if (NN==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (iscolmajor)
    {
        const size_t RC = R*C, RCS = RC*S;
        for (size_t nh=0u; nh<Nh; ++nh)
        {
            for (size_t h=0u; h<H; ++h)
            {
                for (size_t ns=0u; ns<Ns; ++ns)
                {
                    for (size_t s=0u; s<S; ++s)
                    {
                        for (size_t nc=0u; nc<Nc; ++nc)
                        {
                            for (size_t c=0u; c<C; ++c)
                            {
                                for (size_t nr=0u; nr<Nr; ++nr)
                                {
                                    for (size_t r=0u; r<R; ++r, ++X, ++Y)
                                    {
                                        *Y = *X; *++Y = *++X;
                                    }
                                    if (nr<Nr-1u) { X -= 2u*R; }
                                }
                            }
                            if (nc<Nc-1u) { X -= 2u*RC; }
                        }
                    }
                    if (ns<Ns-1u) { X -= 2u*RCS; }
                }
            }
            if (nh<Nh-1u) { X -= 2u*N; }
        }
    }
    else
    {
        const size_t SH = S*H, CSH = C*SH;
        for (size_t nr=0u; nr<Nr; ++nr)
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t nc=0u; nc<Nc; ++nc)
                {
                    for (size_t c=0u; c<C; ++c)
                    {
                        for (size_t ns=0u; ns<Ns; ++ns)
                        {
                            for (size_t s=0u; s<S; ++s)
                            {
                                for (size_t nh=0u; nh<Nh; ++nh)
                                {
                                    for (size_t h=0u; h<H; ++h, ++X, ++Y)
                                    {
                                        *Y = *X; *++Y = *++X;
                                    }
                                    if (nh<Nh-1u) { X -= 2u*H; }
                                }
                            }
                            if (ns<Ns-1u) { X -= 2u*SH; }
                        }
                    }
                    if (nc<Nc-1u) { X -= 2u*CSH; }
                }
            }
            if (nr<Nr-1u) { X -= 2u*N; }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
