//Linear algebra function.
//Singular value decomposition of general matrix X.
//Singular values and vectors are returned in descending order.

//Since only the first K are kept, the first K*R (or K*C) elements hold the result.
//The input matrices must be allocated to their full size!
//That is, when allocating, assume that K is at K_full = min(R,C).

//For complex case, singular values in S are real-valued.

//The inplace version destroys X during processing.
//The not-inplace version keeps X const, but requires full copy of X!

#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int svd_s (float *U, float *S, float *Vt, const float *X, const size_t R, const size_t C, const char iscolmajor, const size_t K);
int svd_d (double *U, double *S, double *Vt, const double *X, const size_t R, const size_t C, const char iscolmajor, const size_t K);
int svd_c (float *U, float *S, float *Vt, const float *X, const size_t R, const size_t C, const char iscolmajor, const size_t K);
int svd_z (double *U, double *S, double *Vt, const double *X, const size_t R, const size_t C, const char iscolmajor, const size_t K);

int svd_inplace_s (float *U, float *S, float *Vt, float *X, const size_t R, const size_t C, const char iscolmajor, const size_t K);
int svd_inplace_d (double *U, double *S, double *Vt, double *X, const size_t R, const size_t C, const char iscolmajor, const size_t K);
int svd_inplace_c (float *U, float *S, float *Vt, float *X, const size_t R, const size_t C, const char iscolmajor, const size_t K);
int svd_inplace_z (double *U, double *S, double *Vt, double *X, const size_t R, const size_t C, const char iscolmajor, const size_t K);


int svd_s (float *U, float *S, float *Vt, const float *X, const size_t R, const size_t C, const char iscolmajor, const size_t K)
{
    const size_t N = R*C;

    if (N==0u) {}
    else
    {
        const size_t Kmax = (R<C) ? R : C;
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        const int ldu = (iscolmajor) ? (int)R : (int)Kmax;
        const int ldv = (iscolmajor) ? (int)Kmax : (int)C;

        float *tmp, *Xtmp;
        if (!(tmp=(float *)malloc(Kmax*sizeof(float)))) { fprintf(stderr,"error in svd_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Xtmp=(float *)malloc(N*sizeof(float)))) { fprintf(stderr,"error in svd_s: problem with malloc. "); perror("malloc"); return 1; }

        for (size_t n=0u; n<N; ++n, ++X, ++Xtmp) { *Xtmp = *X; }
        Xtmp -= N;

        if (LAPACKE_sgesvd(Ord,'S','S',(int)R,(int)C,Xtmp,lda,S,U,ldu,Vt,ldv,tmp))
        { fprintf(stderr,"error in svd_s: problem with LAPACKE_sgesvd function\n"); }

        //Keep K
        if (K>0 && K<Kmax)
        {
            if (iscolmajor)
            {
                for (size_t c=0u; c<C; ++c)
                {
                    for (size_t r=0u; r<K; ++r, ++Vt) { *Vt = *(Vt+c*(Kmax-K)); }
                }
            }
            else
            {
                for (size_t r=0u; r<R; ++r)
                {
                    for (size_t c=0u; c<K; ++c, ++U) { *U = *(U+r*(Kmax-K)); }
                }
            }
        }

        free(tmp); free(Xtmp);
    }

    return 0;
}


int svd_d (double *U, double *S, double *Vt, const double *X, const size_t R, const size_t C, const char iscolmajor, const size_t K)
{
    const size_t N = R*C;

    if (N==0u) {}
    else
    {
        const size_t Kmax = (R<C) ? R : C;
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        const int ldu = (iscolmajor) ? (int)R : (int)Kmax;
        const int ldv = (iscolmajor) ? (int)Kmax : (int)C;

        double *tmp, *Xtmp;
        if (!(tmp=(double *)malloc(Kmax*sizeof(double)))) { fprintf(stderr,"error in svd_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Xtmp=(double *)malloc(N*sizeof(double)))) { fprintf(stderr,"error in svd_s: problem with malloc. "); perror("malloc"); return 1; }

        for (size_t n=0u; n<N; ++n, ++X, ++Xtmp) { *Xtmp = *X; }
        Xtmp -= N;

        if (LAPACKE_dgesvd(Ord,'S','S',(int)R,(int)C,Xtmp,lda,S,U,ldu,Vt,ldv,tmp))
        { fprintf(stderr,"error in svd_d: problem with LAPACKE_sgesvd function\n"); }

        //Keep K
        if (K>0 && K<Kmax)
        {
            if (iscolmajor)
            {
                for (size_t c=0u; c<C; ++c)
                {
                    for (size_t r=0u; r<K; ++r, ++Vt) { *Vt = *(Vt+c*(Kmax-K)); }
                }
            }
            else
            {
                for (size_t r=0u; r<R; ++r)
                {
                    for (size_t c=0u; c<K; ++c, ++U) { *U = *(U+r*(Kmax-K)); }
                }
            }
        }

        free(tmp); free(Xtmp);
    }

    return 0;
}


int svd_c (float *U, float *S, float *Vt, const float *X, const size_t R, const size_t C, const char iscolmajor, const size_t K)
{
    const size_t N = R*C;

    if (N==0u) {}
    else
    {
        const size_t Kmax = (R<C) ? R : C;
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        const int ldu = (iscolmajor) ? (int)R : (int)Kmax;
        const int ldv = (iscolmajor) ? (int)Kmax : (int)C;

        float *tmp, *Xtmp;
        if (!(tmp=(float *)malloc(2*Kmax*sizeof(float)))) { fprintf(stderr,"error in svd_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Xtmp=(float *)malloc(2*N*sizeof(float)))) { fprintf(stderr,"error in svd_s: problem with malloc. "); perror("malloc"); return 1; }

        for (size_t n=0u; n<2*N; ++n, ++X, ++Xtmp) { *Xtmp = *X; }
        Xtmp -= 2*N;

        if (LAPACKE_cgesvd(Ord,'S','S',(int)R,(int)C,(_Complex float *)Xtmp,lda,S,(_Complex float *)U,ldu,(_Complex float *)Vt,ldv,tmp))
        { fprintf(stderr,"error in svd_c: problem with LAPACKE_cgesvd function\n"); }

        //Keep K
        if (K>0 && K<Kmax)
        {
            if (iscolmajor)
            {
                for (size_t c=0u; c<C; ++c)
                {
                    for (size_t r=0u; r<K; ++r, ++Vt) { *Vt = *(Vt+2*c*(Kmax-K)); ++Vt; *Vt = *(Vt+2*c*(Kmax-K)); }
                }
            }
            else
            {
                for (size_t r=0u; r<R; ++r)
                {
                    for (size_t c=0u; c<K; ++c, ++U) { *U = *(U+2*r*(Kmax-K)); ++U; *U = *(U+2*r*(Kmax-K)); }
                }
            }
        }

        free(tmp); free(Xtmp);
    }

    return 0;
}


int svd_z (double *U, double *S, double *Vt, const double *X, const size_t R, const size_t C, const char iscolmajor, const size_t K)
{
    const size_t N = R*C;

    if (N==0u) {}
    else
    {
        const size_t Kmax = (R<C) ? R : C;
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        const int ldu = (iscolmajor) ? (int)R : (int)Kmax;
        const int ldv = (iscolmajor) ? (int)Kmax : (int)C;

        double *tmp, *Xtmp;
        if (!(tmp=(double *)malloc(2*Kmax*sizeof(double)))) { fprintf(stderr,"error in svd_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Xtmp=(double *)malloc(2*N*sizeof(double)))) { fprintf(stderr,"error in svd_s: problem with malloc. "); perror("malloc"); return 1; }

        for (size_t n=0u; n<2*N; ++n, ++X, ++Xtmp) { *Xtmp = *X; }
        Xtmp -= 2*N;

        if (LAPACKE_zgesvd(Ord,'S','S',(int)R,(int)C,(_Complex double *)Xtmp,lda,S,(_Complex double *)U,ldu,(_Complex double *)Vt,ldv,tmp))
        { fprintf(stderr,"error in svd_z: problem with LAPACKE_zgesvd function\n"); }

        //Keep K
        if (K>0 && K<Kmax)
        {
            if (iscolmajor)
            {
                for (size_t c=0u; c<C; ++c)
                {
                    for (size_t r=0u; r<K; ++r, ++Vt) { *Vt = *(Vt+2*c*(Kmax-K)); ++Vt; *Vt = *(Vt+2*c*(Kmax-K)); }
                }
            }
            else
            {
                for (size_t r=0u; r<R; ++r)
                {
                    for (size_t c=0u; c<K; ++c, ++U) { *U = *(U+2*r*(Kmax-K)); ++U; *U = *(U+2*r*(Kmax-K)); }
                }
            }
        }

        free(tmp); free(Xtmp);
    }

    return 0;
}


int svd_inplace_s (float *U, float *S, float *Vt, float *X, const size_t R, const size_t C, const char iscolmajor, const size_t K)
{
    if (R*C==0) {}
    else
    {
        const size_t Kmax = (R<C) ? R : C;
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        const int ldu = (iscolmajor) ? (int)R : (int)Kmax;
        const int ldv = (iscolmajor) ? (int)Kmax : (int)C;

        float *tmp;
        if (!(tmp=(float *)malloc(Kmax*sizeof(float)))) { fprintf(stderr,"error in svd_s: problem with malloc. "); perror("malloc"); return 1; }

        if (LAPACKE_sgesvd(Ord,'S','S',(int)R,(int)C,X,lda,S,U,ldu,Vt,ldv,tmp))
        { fprintf(stderr,"error in svd_s: problem with LAPACKE_sgesvd function\n"); }

        //Keep K
        if (K>0 && K<Kmax)
        {
            if (iscolmajor)
            {
                for (size_t c=0u; c<C; ++c)
                {
                    for (size_t r=0u; r<K; ++r, ++Vt) { *Vt = *(Vt+c*(Kmax-K)); }
                }
            }
            else
            {
                for (size_t r=0u; r<R; ++r)
                {
                    for (size_t c=0u; c<K; ++c, ++U) { *U = *(U+r*(Kmax-K)); }
                }
            }
        }

        free(tmp);
    }

    return 0;
}


int svd_inplace_d (double *U, double *S, double *Vt, double *X, const size_t R, const size_t C, const char iscolmajor, const size_t K)
{
    if (R*C==0) {}
    else
    {
        const size_t Kmax = (R<C) ? R : C;
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        const int ldu = (iscolmajor) ? (int)R : (int)Kmax;
        const int ldv = (iscolmajor) ? (int)Kmax : (int)C;

        double *tmp;
        if (!(tmp=(double *)malloc(Kmax*sizeof(double)))) { fprintf(stderr,"error in svd_d: problem with malloc. "); perror("malloc"); return 1; }

        if (LAPACKE_dgesvd(Ord,'S','S',(int)R,(int)C,X,lda,S,U,ldu,Vt,ldv,tmp))
        { fprintf(stderr,"error in svd_d: problem with LAPACKE_sgesvd function\n"); }

        //Keep K
        if (K>0 && K<Kmax)
        {
            if (iscolmajor)
            {
                for (size_t c=0u; c<C; ++c)
                {
                    for (size_t r=0u; r<K; ++r, ++Vt) { *Vt = *(Vt+c*(Kmax-K)); }
                }
            }
            else
            {
                for (size_t r=0u; r<R; ++r)
                {
                    for (size_t c=0u; c<K; ++c, ++U) { *U = *(U+r*(Kmax-K)); }
                }
            }
        }

        free(tmp);
    }

    return 0;
}


int svd_inplace_c (float *U, float *S, float *Vt, float *X, const size_t R, const size_t C, const char iscolmajor, const size_t K)
{
    if (R*C==0) {}
    else
    {
        const size_t Kmax = (R<C) ? R : C;
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        const int ldu = (iscolmajor) ? (int)R : (int)Kmax;
        const int ldv = (iscolmajor) ? (int)Kmax : (int)C;

        float *tmp;
        if (!(tmp=(float *)malloc(2*Kmax*sizeof(float)))) { fprintf(stderr,"error in svd_c: problem with malloc. "); perror("malloc"); return 1; }

        if (LAPACKE_cgesvd(Ord,'S','S',(int)R,(int)C,(_Complex float *)X,lda,S,(_Complex float *)U,ldu,(_Complex float *)Vt,ldv,tmp))
        { fprintf(stderr,"error in svd_c: problem with LAPACKE_cgesvd function\n"); }

        //Keep K
        if (K>0 && K<Kmax)
        {
            if (iscolmajor)
            {
                for (size_t c=0u; c<C; ++c)
                {
                    for (size_t r=0u; r<K; ++r, ++Vt) { *Vt = *(Vt+2*c*(Kmax-K)); ++Vt; *Vt = *(Vt+2*c*(Kmax-K)); }
                }
            }
            else
            {
                for (size_t r=0u; r<R; ++r)
                {
                    for (size_t c=0u; c<K; ++c, ++U) { *U = *(U+2*r*(Kmax-K)); ++U; *U = *(U+2*r*(Kmax-K)); }
                }
            }
        }

        free(tmp);
    }

    return 0;
}


int svd_inplace_z (double *U, double *S, double *Vt, double *X, const size_t R, const size_t C, const char iscolmajor, const size_t K)
{
    if (R*C==0) {}
    else
    {
        const size_t Kmax = (R<C) ? R : C;
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        const int ldu = (iscolmajor) ? (int)R : (int)Kmax;
        const int ldv = (iscolmajor) ? (int)Kmax : (int)C;
        
        double *tmp;
        if (!(tmp=(double *)malloc(2*Kmax*sizeof(double)))) { fprintf(stderr,"error in svd_z: problem with malloc. "); perror("malloc"); return 1; }

        if (LAPACKE_zgesvd(Ord,'S','S',(int)R,(int)C,(_Complex double *)X,lda,S,(_Complex double *)U,ldu,(_Complex double *)Vt,ldv,tmp))
        { fprintf(stderr,"error in svd_z: problem with LAPACKE_zgesvd function\n"); }

        //Keep K
        if (K>0 && K<Kmax)
        {
            if (iscolmajor)
            {
                for (size_t c=0u; c<C; ++c)
                {
                    for (size_t r=0u; r<K; ++r, ++Vt) { *Vt = *(Vt+2*c*(Kmax-K)); ++Vt; *Vt = *(Vt+2*c*(Kmax-K)); }
                }
            }
            else
            {
                for (size_t r=0u; r<R; ++r)
                {
                    for (size_t c=0u; c<K; ++c, ++U) { *U = *(U+2*r*(Kmax-K)); ++U; *U = *(U+2*r*(Kmax-K)); }
                }
            }
        }

        free(tmp);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
