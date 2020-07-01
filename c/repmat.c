//Takes input tensor X and repeats it Nr x Nc x Ns x Nh times,
//so that the output Y has size R*Nr x C*Nc x S*Ns x H*Nh.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int repmat_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const int Nr, const int Nc, const int Ns, const int Nh, const char iscolmajor);
int repmat_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const int Nr, const int Nc, const int Ns, const int Nh, const char iscolmajor);
int repmat_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const int Nr, const int Nc, const int Ns, const int Nh, const char iscolmajor);
int repmat_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const int Nr, const int Nc, const int Ns, const int Nh, const char iscolmajor);


int repmat_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const int Nr, const int Nc, const int Ns, const int Nh, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in repmat_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in repmat_s: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in repmat_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in repmat_s: H (num hyperslices X) must be positive\n"); return 1; }
    if (Nr<1) { fprintf(stderr,"error in repmat_s: Nr (num row reps) must be positive\n"); return 1; }
    if (Nc<1) { fprintf(stderr,"error in repmat_s: Nc (num col reps) must be positive\n"); return 1; }
    if (Ns<1) { fprintf(stderr,"error in repmat_s: Ns (num slice reps) must be positive\n"); return 1; }
    if (Nh<1) { fprintf(stderr,"error in repmat_s: Nh (num hyperslice reps) must be positive\n"); return 1; }

    const int CR = C*R, SCR = S*CR;
    const int SH = S*H, CSH = C*SH;
    int n, n2 = 0;

    if (iscolmajor)
    {
        if (R*Nr==1)
        {
            for (int h=0; h<H*Nh; h++)
            {
                for (int s=0; s<S*Ns; s++)
                {
                    n = SCR*h%H + CR*s%S;
                    for (int c=0; c<Nc; c++, n2+=C) { cblas_scopy(C,&X[n],1,&Y[n2],1); }
                }
            }
        }
        else
        {
            for (int h=0; h<H*Nh; h++)
            {
                for (int s=0; s<S*Ns; s++)
                {
                    for (int c=0; c<C*Nc; c++)
                    {
                        n = SCR*h%H + CR*s%S + R*c%C;
                        for (int r=0; r<Nr; r++, n2+=R) { cblas_scopy(R,&X[n],1,&Y[n2],1); }
                    }
                }
            }
        }
    }
    else
    {
        if (CSH*Nc*Ns*Nh==1)
        {
            for (int r=0; r<Nr; r++, n2+=R) { cblas_scopy(R,X,1,&Y[n2],1); }
        }
        else if (SH*Ns*Nh==1)
        {
            for (int r=0; r<R*Nr; r++)
            {
                n = CSH*r%R;
                for (int c=0; c<Nc; c++, n2+=C) { cblas_scopy(C,&X[n],1,&Y[n2],1); }
            }
        }
        else
        {
            for (int r=0; r<R*Nr; r++)
            {
                for (int c=0; c<C*Nc; c++)
                {
                    for (int s=0; s<S*Ns; s++)
                    {
                        n = CSH*r%R + SH*c%C + H*s%S;
                        for (int h=0; h<Nh; h++, n2+=H) { cblas_scopy(H,&X[n],1,&Y[n2],1); }
                    }
                }
            }
        }
    }

    return 0;
}


int repmat_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const int Nr, const int Nc, const int Ns, const int Nh, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in repmat_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in repmat_d: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in repmat_d: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in repmat_d: H (num hyperslices X) must be positive\n"); return 1; }
    if (Nr<1) { fprintf(stderr,"error in repmat_d: Nr (num row reps) must be positive\n"); return 1; }
    if (Nc<1) { fprintf(stderr,"error in repmat_d: Nc (num col reps) must be positive\n"); return 1; }
    if (Ns<1) { fprintf(stderr,"error in repmat_d: Ns (num slice reps) must be positive\n"); return 1; }
    if (Nh<1) { fprintf(stderr,"error in repmat_d: Nh (num hyperslice reps) must be positive\n"); return 1; }

    const int CR = C*R, SCR = S*CR;
    const int SH = S*H, CSH = C*SH;
    int n, n2 = 0;

    if (iscolmajor)
    {
        if (R*Nr==1)
        {
            for (int h=0; h<H*Nh; h++)
            {
                for (int s=0; s<S*Ns; s++)
                {
                    n = SCR*h%H + CR*s%S;
                    for (int c=0; c<Nc; c++, n2+=C) { cblas_dcopy(C,&X[n],1,&Y[n2],1); }
                }
            }
        }
        else
        {
            for (int h=0; h<H*Nh; h++)
            {
                for (int s=0; s<S*Ns; s++)
                {
                    for (int c=0; c<C*Nc; c++)
                    {
                        n = SCR*h%H + CR*s%S + R*c%C;
                        for (int r=0; r<Nr; r++, n2+=R) { cblas_dcopy(R,&X[n],1,&Y[n2],1); }
                    }
                }
            }
        }
    }
    else
    {
        if (CSH*Nc*Ns*Nh==1)
        {
            for (int r=0; r<Nr; r++, n2+=R) { cblas_dcopy(R,X,1,&Y[n2],1); }
        }
        else if (SH*Ns*Nh==1)
        {
            for (int r=0; r<R*Nr; r++)
            {
                n = CSH*r%R;
                for (int c=0; c<Nc; c++, n2+=C) { cblas_dcopy(C,&X[n],1,&Y[n2],1); }
            }
        }
        else
        {
            for (int r=0; r<R*Nr; r++)
            {
                for (int c=0; c<C*Nc; c++)
                {
                    for (int s=0; s<S*Ns; s++)
                    {
                        n = CSH*r%R + SH*c%C + H*s%S;
                        for (int h=0; h<Nh; h++, n2+=H) { cblas_dcopy(H,&X[n],1,&Y[n2],1); }
                    }
                }
            }
        }
    }

    return 0;
}


int repmat_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const int Nr, const int Nc, const int Ns, const int Nh, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in repmat_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in repmat_c: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in repmat_c: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in repmat_c: H (num hyperslices X) must be positive\n"); return 1; }
    if (Nr<1) { fprintf(stderr,"error in repmat_c: Nr (num row reps) must be positive\n"); return 1; }
    if (Nc<1) { fprintf(stderr,"error in repmat_c: Nc (num col reps) must be positive\n"); return 1; }
    if (Ns<1) { fprintf(stderr,"error in repmat_c: Ns (num slice reps) must be positive\n"); return 1; }
    if (Nh<1) { fprintf(stderr,"error in repmat_c: Nh (num hyperslice reps) must be positive\n"); return 1; }

    const int CR = C*R, SCR = S*CR;
    const int SH = S*H, CSH = C*SH;
    int n, n2 = 0;

    if (iscolmajor)
    {
        if (R*Nr==1)
        {
            for (int h=0; h<H*Nh; h++)
            {
                for (int s=0; s<S*Ns; s++)
                {
                    n = 2*(SCR*h%H + CR*s%S);
                    for (int c=0; c<Nc; c++, n2+=2*C) { cblas_ccopy(C,&X[n],1,&Y[n2],1); }
                }
            }
        }
        else
        {
            for (int h=0; h<H*Nh; h++)
            {
                for (int s=0; s<S*Ns; s++)
                {
                    for (int c=0; c<C*Nc; c++)
                    {
                        n = 2*(SCR*h%H + CR*s%S + R*c%C);
                        for (int r=0; r<Nr; r++, n2+=2*R) { cblas_ccopy(R,&X[n],1,&Y[n2],1); }
                    }
                }
            }
        }
    }
    else
    {
        if (CSH*Nc*Ns*Nh==1)
        {
            for (int r=0; r<Nr; r++, n2+=2*R) { cblas_ccopy(R,X,1,&Y[n2],1); }
        }
        else if (SH*Ns*Nh==1)
        {
            for (int r=0; r<R*Nr; r++)
            {
                n = 2*(CSH*r%R);
                for (int c=0; c<Nc; c++, n2+=2*C) { cblas_ccopy(C,&X[n],1,&Y[n2],1); }
            }
        }
        else
        {
            for (int r=0; r<R*Nr; r++)
            {
                for (int c=0; c<C*Nc; c++)
                {
                    for (int s=0; s<S*Ns; s++)
                    {
                        n = 2*(CSH*r%R + SH*c%C + H*s%S);
                        for (int h=0; h<Nh; h++, n2+=2*H) { cblas_ccopy(H,&X[n],1,&Y[n2],1); }
                    }
                }
            }
        }
    }

    return 0;
}


int repmat_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const int Nr, const int Nc, const int Ns, const int Nh, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in repmat_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in repmat_z: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in repmat_z: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in repmat_z: H (num hyperslices X) must be positive\n"); return 1; }
    if (Nr<1) { fprintf(stderr,"error in repmat_z: Nr (num row reps) must be positive\n"); return 1; }
    if (Nc<1) { fprintf(stderr,"error in repmat_z: Nc (num col reps) must be positive\n"); return 1; }
    if (Ns<1) { fprintf(stderr,"error in repmat_z: Ns (num slice reps) must be positive\n"); return 1; }
    if (Nh<1) { fprintf(stderr,"error in repmat_z: Nh (num hyperslice reps) must be positive\n"); return 1; }

    const int CR = C*R, SCR = S*CR;
    const int SH = S*H, CSH = C*SH;
    int n, n2 = 0;

    if (iscolmajor)
    {
        if (R*Nr==1)
        {
            for (int h=0; h<H*Nh; h++)
            {
                for (int s=0; s<S*Ns; s++)
                {
                    n = 2*(SCR*h%H + CR*s%S);
                    for (int c=0; c<Nc; c++, n2+=2*C) { cblas_zcopy(C,&X[n],1,&Y[n2],1); }
                }
            }
        }
        else
        {
            for (int h=0; h<H*Nh; h++)
            {
                for (int s=0; s<S*Ns; s++)
                {
                    for (int c=0; c<C*Nc; c++)
                    {
                        n = 2*(SCR*h%H + CR*s%S + R*c%C);
                        for (int r=0; r<Nr; r++, n2+=2*R) { cblas_zcopy(R,&X[n],1,&Y[n2],1); }
                    }
                }
            }
        }
    }
    else
    {
        if (CSH*Nc*Ns*Nh==1)
        {
            for (int r=0; r<Nr; r++, n2+=2*R) { cblas_zcopy(R,X,1,&Y[n2],1); }
        }
        else if (SH*Ns*Nh==1)
        {
            for (int r=0; r<R*Nr; r++)
            {
                n = 2*(CSH*r%R);
                for (int c=0; c<Nc; c++, n2+=2*C) { cblas_zcopy(C,&X[n],1,&Y[n2],1); }
            }
        }
        else
        {
            for (int r=0; r<R*Nr; r++)
            {
                for (int c=0; c<C*Nc; c++)
                {
                    for (int s=0; s<S*Ns; s++)
                    {
                        n = 2*(CSH*r%R + SH*c%C + H*s%S);
                        for (int h=0; h<Nh; h++, n2+=2*H) { cblas_zcopy(H,&X[n],1,&Y[n2],1); }
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
