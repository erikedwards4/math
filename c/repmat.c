//Takes input tensor X and repeats it Nr x Nc x Ns x Nh times,
//so that the output Y has size R*Nr x C*Nc x S*Ns x H*Nh.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int repmat_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const char iscolmajor);
int repmat_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const char iscolmajor);
int repmat_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const char iscolmajor);
int repmat_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const char iscolmajor);


int repmat_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const char iscolmajor)
{
    const size_t N = R*C*S*H, NN = Nr*Nc*Ns*Nh;

    if (N==0 || NN==0) {}
    else if (NN==1)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (iscolmajor)
    {
        const size_t RC = R*C, RCS = RC*S;
        for (size_t nh=0; nh<Nh; ++nh)
        {
            for (size_t h=0u; h<H; ++h)
            {
                for (size_t ns=0; ns<Ns; ++ns)
                {
                    for (size_t s=0u; s<S; ++s)
                    {
                        for (size_t nc=0; nc<Nc; ++nc)
                        {
                            for (size_t c=0u; c<C; ++c)
                            {
                                for (size_t nr=0; nr<Nr; ++nr)
                                {
                                    for (size_t r=0u; r<R; ++r, ++X, ++Y)
                                    {
                                        *Y = *X;
                                    }
                                    if (nr<Nr-1) { X -= R; }
                                }
                            }
                            if (nc<Nc-1) { X -= RC; }
                        }
                    }
                    if (ns<Ns-1) { X -= RCS; }
                }
            }
            if (nh<Nh-1) { X -= N; }
        }
    }
    else
    {
        const size_t SH = S*H, CSH = C*SH;
        for (size_t nr=0; nr<Nr; ++nr)
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t nc=0; nc<Nc; ++nc)
                {
                    for (size_t c=0u; c<C; ++c)
                    {
                        for (size_t ns=0; ns<Ns; ++ns)
                        {
                            for (size_t s=0u; s<S; ++s)
                            {
                                for (size_t nh=0; nh<Nh; ++nh)
                                {
                                    for (size_t h=0u; h<H; ++h, ++X, ++Y)
                                    {
                                        *Y = *X;
                                    }
                                    if (nh<Nh-1) { X -= H; }
                                }
                            }
                            if (ns<Ns-1) { X -= SH; }
                        }
                    }
                    if (nc<Nc-1) { X -= CSH; }
                }
            }
            if (nr<Nr-1) { X -= N; }
        }
    }

    return 0;
}


int repmat_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const char iscolmajor)
{
    const size_t N = R*C*S*H, NN = Nr*Nc*Ns*Nh;

    if (N==0 || NN==0) {}
    else if (NN==1)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (iscolmajor)
    {
        const size_t RC = R*C, RCS = RC*S;
        for (size_t nh=0; nh<Nh; ++nh)
        {
            for (size_t h=0u; h<H; ++h)
            {
                for (size_t ns=0; ns<Ns; ++ns)
                {
                    for (size_t s=0u; s<S; ++s)
                    {
                        for (size_t nc=0; nc<Nc; ++nc)
                        {
                            for (size_t c=0u; c<C; ++c)
                            {
                                for (size_t nr=0; nr<Nr; ++nr)
                                {
                                    for (size_t r=0u; r<R; ++r, ++X, ++Y)
                                    {
                                        *Y = *X;
                                    }
                                    if (nr<Nr-1) { X -= R; }
                                }
                            }
                            if (nc<Nc-1) { X -= RC; }
                        }
                    }
                    if (ns<Ns-1) { X -= RCS; }
                }
            }
            if (nh<Nh-1) { X -= N; }
        }
    }
    else
    {
        const size_t SH = S*H, CSH = C*SH;
        for (size_t nr=0; nr<Nr; ++nr)
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t nc=0; nc<Nc; ++nc)
                {
                    for (size_t c=0u; c<C; ++c)
                    {
                        for (size_t ns=0; ns<Ns; ++ns)
                        {
                            for (size_t s=0u; s<S; ++s)
                            {
                                for (size_t nh=0; nh<Nh; ++nh)
                                {
                                    for (size_t h=0u; h<H; ++h, ++X, ++Y)
                                    {
                                        *Y = *X;
                                    }
                                    if (nh<Nh-1) { X -= H; }
                                }
                            }
                            if (ns<Ns-1) { X -= SH; }
                        }
                    }
                    if (nc<Nc-1) { X -= CSH; }
                }
            }
            if (nr<Nr-1) { X -= N; }
        }
    }

    return 0;
}


int repmat_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const char iscolmajor)
{
    const size_t N = R*C*S*H, NN = Nr*Nc*Ns*Nh;

    if (N==0 || NN==0) {}
    else if (NN==1)
    {
        for (size_t n=0u; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (iscolmajor)
    {
        const size_t RC = R*C, RCS = RC*S;
        for (size_t nh=0; nh<Nh; ++nh)
        {
            for (size_t h=0u; h<H; ++h)
            {
                for (size_t ns=0; ns<Ns; ++ns)
                {
                    for (size_t s=0u; s<S; ++s)
                    {
                        for (size_t nc=0; nc<Nc; ++nc)
                        {
                            for (size_t c=0u; c<C; ++c)
                            {
                                for (size_t nr=0; nr<Nr; ++nr)
                                {
                                    for (size_t r=0u; r<R; ++r, ++X, ++Y)
                                    {
                                        *Y = *X; *++Y = *++X;
                                    }
                                    if (nr<Nr-1) { X -= 2*R; }
                                }
                            }
                            if (nc<Nc-1) { X -= 2*RC; }
                        }
                    }
                    if (ns<Ns-1) { X -= 2*RCS; }
                }
            }
            if (nh<Nh-1) { X -= 2*N; }
        }
    }
    else
    {
        const size_t SH = S*H, CSH = C*SH;
        for (size_t nr=0; nr<Nr; ++nr)
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t nc=0; nc<Nc; ++nc)
                {
                    for (size_t c=0u; c<C; ++c)
                    {
                        for (size_t ns=0; ns<Ns; ++ns)
                        {
                            for (size_t s=0u; s<S; ++s)
                            {
                                for (size_t nh=0; nh<Nh; ++nh)
                                {
                                    for (size_t h=0u; h<H; ++h, ++X, ++Y)
                                    {
                                        *Y = *X; *++Y = *++X;
                                    }
                                    if (nh<Nh-1) { X -= 2*H; }
                                }
                            }
                            if (ns<Ns-1) { X -= 2*SH; }
                        }
                    }
                    if (nc<Nc-1) { X -= 2*CSH; }
                }
            }
            if (nr<Nr-1) { X -= 2*N; }
        }
    }

    return 0;
}


int repmat_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const char iscolmajor)
{
    const size_t N = R*C*S*H, NN = Nr*Nc*Ns*Nh;

    if (N==0 || NN==0) {}
    else if (NN==1)
    {
        for (size_t n=0u; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (iscolmajor)
    {
        const size_t RC = R*C, RCS = RC*S;
        for (size_t nh=0; nh<Nh; ++nh)
        {
            for (size_t h=0u; h<H; ++h)
            {
                for (size_t ns=0; ns<Ns; ++ns)
                {
                    for (size_t s=0u; s<S; ++s)
                    {
                        for (size_t nc=0; nc<Nc; ++nc)
                        {
                            for (size_t c=0u; c<C; ++c)
                            {
                                for (size_t nr=0; nr<Nr; ++nr)
                                {
                                    for (size_t r=0u; r<R; ++r, ++X, ++Y)
                                    {
                                        *Y = *X; *++Y = *++X;
                                    }
                                    if (nr<Nr-1) { X -= 2*R; }
                                }
                            }
                            if (nc<Nc-1) { X -= 2*RC; }
                        }
                    }
                    if (ns<Ns-1) { X -= 2*RCS; }
                }
            }
            if (nh<Nh-1) { X -= 2*N; }
        }
    }
    else
    {
        const size_t SH = S*H, CSH = C*SH;
        for (size_t nr=0; nr<Nr; ++nr)
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t nc=0; nc<Nc; ++nc)
                {
                    for (size_t c=0u; c<C; ++c)
                    {
                        for (size_t ns=0; ns<Ns; ++ns)
                        {
                            for (size_t s=0u; s<S; ++s)
                            {
                                for (size_t nh=0; nh<Nh; ++nh)
                                {
                                    for (size_t h=0u; h<H; ++h, ++X, ++Y)
                                    {
                                        *Y = *X; *++Y = *++X;
                                    }
                                    if (nh<Nh-1) { X -= 2*H; }
                                }
                            }
                            if (ns<Ns-1) { X -= 2*SH; }
                        }
                    }
                    if (nc<Nc-1) { X -= 2*CSH; }
                }
            }
            if (nr<Nr-1) { X -= 2*N; }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
