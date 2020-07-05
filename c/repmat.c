//Takes input tensor X and repeats it Nr x Nc x Ns x Nh times,
//so that the output Y has size R*Nr x C*Nc x S*Ns x H*Nh.

#include <stdio.h>
#include <cblas.h>

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
    const size_t CR = C*R, SCR = S*CR;
    const size_t SH = S*H, CSH = C*SH;

    if (Nr*Nc*Ns*Nh==1)
    {
        cblas_scopy((int)(CR*SH),X,1,Y,1);
    }
    else if (iscolmajor)
    {
        if (R*Nr==1)
        {
            for (size_t h=0, n2=0; h<H*Nh; h++)
            {
                for (size_t s=0; s<S*Ns; s++)
                {
                    size_t n = SCR*h%H + CR*s%S;
                    for (size_t c=0; c<Nc; c++, n2+=C) { cblas_scopy((int)C,&X[n],1,&Y[n2],1); }
                }
            }
        }
        else
        {
            for (size_t h=0, n2=0; h<H*Nh; h++)
            {
                for (size_t s=0; s<S*Ns; s++)
                {
                    for (size_t c=0; c<C*Nc; c++)
                    {
                        size_t n = SCR*h%H + CR*s%S + R*c%C;
                        for (size_t r=0; r<Nr; r++, n2+=R) { cblas_scopy((int)R,&X[n],1,&Y[n2],1); }
                    }
                }
            }
        }
    }
    else
    {
        if (CSH*Nc*Ns*Nh==1)
        {
            for (size_t r=0, n2=0; r<Nr; r++, n2+=R) { cblas_scopy((int)R,X,1,&Y[n2],1); }
        }
        else if (SH*Ns*Nh==1)
        {
            for (size_t r=0, n2=0; r<R*Nr; r++)
            {
                size_t n = CSH*r%R;
                for (size_t c=0; c<Nc; c++, n2+=C) { cblas_scopy((int)C,&X[n],1,&Y[n2],1); }
            }
        }
        else
        {
            for (size_t r=0, n2=0; r<R*Nr; r++)
            {
                for (size_t c=0; c<C*Nc; c++)
                {
                    for (size_t s=0; s<S*Ns; s++)
                    {
                        size_t n = CSH*r%R + SH*c%C + H*s%S;
                        for (size_t h=0; h<Nh; h++, n2+=H) { cblas_scopy((int)H,&X[n],1,&Y[n2],1); }
                    }
                }
            }
        }
    }

    return 0;
}


int repmat_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const char iscolmajor)
{
    const size_t CR = C*R, SCR = S*CR;
    const size_t SH = S*H, CSH = C*SH;

    if (Nr*Nc*Ns*Nh==1)
    {
        cblas_dcopy((int)(CR*SH),X,1,Y,1);
    }
    else if (iscolmajor)
    {
        if (R*Nr==1)
        {
            for (size_t h=0, n2=0; h<H*Nh; h++)
            {
                for (size_t s=0; s<S*Ns; s++)
                {
                    size_t n = SCR*h%H + CR*s%S;
                    for (size_t c=0; c<Nc; c++, n2+=C) { cblas_dcopy((int)C,&X[n],1,&Y[n2],1); }
                }
            }
        }
        else
        {
            for (size_t h=0, n2=0; h<H*Nh; h++)
            {
                for (size_t s=0; s<S*Ns; s++)
                {
                    for (size_t c=0; c<C*Nc; c++)
                    {
                        size_t n = SCR*h%H + CR*s%S + R*c%C;
                        for (size_t r=0; r<Nr; r++, n2+=R) { cblas_dcopy((int)R,&X[n],1,&Y[n2],1); }
                    }
                }
            }
        }
    }
    else
    {
        if (CSH*Nc*Ns*Nh==1)
        {
            for (size_t r=0, n2=0; r<Nr; r++, n2+=R) { cblas_dcopy((int)R,X,1,&Y[n2],1); }
        }
        else if (SH*Ns*Nh==1)
        {
            for (size_t r=0, n2=0; r<R*Nr; r++)
            {
                size_t n = CSH*r%R;
                for (size_t c=0; c<Nc; c++, n2+=C) { cblas_dcopy((int)C,&X[n],1,&Y[n2],1); }
            }
        }
        else
        {
            for (size_t r=0, n2=0; r<R*Nr; r++)
            {
                for (size_t c=0; c<C*Nc; c++)
                {
                    for (size_t s=0; s<S*Ns; s++)
                    {
                        size_t n = CSH*r%R + SH*c%C + H*s%S;
                        for (size_t h=0; h<Nh; h++, n2+=H) { cblas_dcopy((int)H,&X[n],1,&Y[n2],1); }
                    }
                }
            }
        }
    }

    return 0;
}


int repmat_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const char iscolmajor)
{
    const size_t CR = C*R, SCR = S*CR;
    const size_t SH = S*H, CSH = C*SH;

    if (Nr*Nc*Ns*Nh==1)
    {
        cblas_ccopy((int)(CR*SH),X,1,Y,1);
    }
    else if (iscolmajor)
    {
        if (R*Nr==1)
        {
            for (size_t h=0, n2=0; h<H*Nh; h++)
            {
                for (size_t s=0; s<S*Ns; s++)
                {
                    size_t n = 2*(SCR*h%H + CR*s%S);
                    for (size_t c=0; c<Nc; c++, n2+=2*C) { cblas_ccopy((int)C,&X[n],1,&Y[n2],1); }
                }
            }
        }
        else
        {
            for (size_t h=0, n2=0; h<H*Nh; h++)
            {
                for (size_t s=0; s<S*Ns; s++)
                {
                    for (size_t c=0; c<C*Nc; c++)
                    {
                        size_t n = 2*(SCR*h%H + CR*s%S + R*c%C);
                        for (size_t r=0; r<Nr; r++, n2+=2*R) { cblas_ccopy((int)R,&X[n],1,&Y[n2],1); }
                    }
                }
            }
        }
    }
    else
    {
        if (CSH*Nc*Ns*Nh==1)
        {
            for (size_t r=0, n2=0; r<Nr; r++, n2+=2*R) { cblas_ccopy((int)R,X,1,&Y[n2],1); }
        }
        else if (SH*Ns*Nh==1)
        {
            for (size_t r=0, n2=0; r<R*Nr; r++)
            {
                size_t n = 2*(CSH*r%R);
                for (size_t c=0; c<Nc; c++, n2+=2*C) { cblas_ccopy((int)C,&X[n],1,&Y[n2],1); }
            }
        }
        else
        {
            for (size_t r=0, n2=0; r<R*Nr; r++)
            {
                for (size_t c=0; c<C*Nc; c++)
                {
                    for (size_t s=0; s<S*Ns; s++)
                    {
                        size_t n = 2*(CSH*r%R + SH*c%C + H*s%S);
                        for (size_t h=0; h<Nh; h++, n2+=2*H) { cblas_ccopy((int)H,&X[n],1,&Y[n2],1); }
                    }
                }
            }
        }
    }

    return 0;
}


int repmat_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Nr, const size_t Nc, const size_t Ns, const size_t Nh, const char iscolmajor)
{
    const size_t CR = C*R, SCR = S*CR;
    const size_t SH = S*H, CSH = C*SH;

    if (Nr*Nc*Ns*Nh==1)
    {
        cblas_zcopy((int)(CR*SH),X,1,Y,1);
    }
    else if (iscolmajor)
    {
        if (R*Nr==1)
        {
            for (size_t h=0, n2=0; h<H*Nh; h++)
            {
                for (size_t s=0; s<S*Ns; s++)
                {
                    size_t n = 2*(SCR*h%H + CR*s%S);
                    for (size_t c=0; c<Nc; c++, n2+=2*C) { cblas_zcopy((int)C,&X[n],1,&Y[n2],1); }
                }
            }
        }
        else
        {
            for (size_t h=0, n2=0; h<H*Nh; h++)
            {
                for (size_t s=0; s<S*Ns; s++)
                {
                    for (size_t c=0; c<C*Nc; c++)
                    {
                        size_t n = 2*(SCR*h%H + CR*s%S + R*c%C);
                        for (size_t r=0; r<Nr; r++, n2+=2*R) { cblas_zcopy((int)R,&X[n],1,&Y[n2],1); }
                    }
                }
            }
        }
    }
    else
    {
        if (CSH*Nc*Ns*Nh==1)
        {
            for (size_t r=0, n2=0; r<Nr; r++, n2+=2*R) { cblas_zcopy((int)R,X,1,&Y[n2],1); }
        }
        else if (SH*Ns*Nh==1)
        {
            for (size_t r=0, n2=0; r<R*Nr; r++)
            {
                size_t n = 2*(CSH*r%R);
                for (size_t c=0; c<Nc; c++, n2+=2*C) { cblas_zcopy((int)C,&X[n],1,&Y[n2],1); }
            }
        }
        else
        {
            for (size_t r=0, n2=0; r<R*Nr; r++)
            {
                for (size_t c=0; c<C*Nc; c++)
                {
                    for (size_t s=0; s<S*Ns; s++)
                    {
                        size_t n = 2*(CSH*r%R + SH*c%C + H*s%S);
                        for (size_t h=0; h<Nh; h++, n2+=2*H) { cblas_zcopy((int)H,&X[n],1,&Y[n2],1); }
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
