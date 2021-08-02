//Linear algebra function.
//Cholesky decomposition of square, Hermitian, positive-definite matrix X.
//This has in-place and not-in-place versions.

//The _work version can be used just as well; it is slightly faster by skipping NaN check,
//but I have decided to use the the NaN-check versions for eig and svd also.

#include <stdio.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int chol_s (float *Y, const float *X, const size_t R, const int iscolmajor, const int upper);
int chol_d (double *Y, const double *X, const size_t R, const int iscolmajor, const int upper);
int chol_c (float *Y, const float *X, const size_t R, const int iscolmajor, const int upper);
int chol_z (double *Y, const double *X, const size_t R, const int iscolmajor, const int upper);

int chol_inplace_s (float *X, const size_t R, const int iscolmajor, const int upper);
int chol_inplace_d (double *X, const size_t R, const int iscolmajor, const int upper);
int chol_inplace_c (float *X, const size_t R, const int iscolmajor, const int upper);
int chol_inplace_z (double *X, const size_t R, const int iscolmajor, const int upper);


int chol_s (float *Y, const float *X, const size_t R, const int iscolmajor, const int upper)
{
    const size_t N = R*R;

    if (N==0u) {}
    else if (N==1u) { *Y = *X; }
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int Uplo = (upper) ? 'U' : 'L';
        
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
        Y -= N;
        
        if (LAPACKE_spotrf(Ord,Uplo,(int)R,Y,(int)R))
        { fprintf(stderr,"error in chol_s: problem with LAPACKE_spotrf function\n"); }

        if ((iscolmajor && upper) || (!iscolmajor && !upper))
        {
            for (size_t r1=0u; r1<R; ++r1)
            {
                Y += r1 + 1u;
                for (size_t r2=r1+1u; r2<R; ++r2, ++Y) { *Y = 0.0f; }
            }
        }
        else
        {
            for (size_t r1=0u; r1<R; ++r1, Y+=R-r1+1u)
            {
                for (size_t r2=0u; r2<r1; ++r2, ++Y) { *Y = 0.0f; }
            }
        }
    }

    return 0;
}


int chol_d (double *Y, const double *X, const size_t R, const int iscolmajor, const int upper)
{
    const size_t N = R*R;

    if (N==0u) {}
    else if (N==1u) { *Y = *X; }
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int Uplo = (upper) ? 'U' : 'L';
        
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
        Y -= N;
        
        if (LAPACKE_dpotrf(Ord,Uplo,(int)R,Y,(int)R))
        { fprintf(stderr,"error in chol_d: problem with LAPACKE_dpotrf function\n"); }

        if ((iscolmajor && upper) || (!iscolmajor && !upper))
        {
            for (size_t r1=0u; r1<R; ++r1)
            {
                Y += r1 + 1u;
                for (size_t r2=r1+1u; r2<R; ++r2, ++Y) { *Y = 0.0; }
            }
        }
        else
        {
            for (size_t r1=0u; r1<R; ++r1, Y+=R-r1+1u)
            {
                for (size_t r2=0u; r2<r1; ++r2, ++Y) { *Y = 0.0; }
            }
        }
    }

    return 0;
}


int chol_c (float *Y, const float *X, const size_t R, const int iscolmajor, const int upper)
{
    const size_t N = R*R;

    if (N==0u) {}
    else if (N==1u) { *Y = *X; *++Y = *++X; }
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int Uplo = (upper) ? 'U' : 'L';
        
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
        Y -= 2u*N;
        
        if (LAPACKE_cpotrf(Ord,Uplo,(int)R,(_Complex float *)Y,(int)R))
        { fprintf(stderr,"error in chol_c: problem with LAPACKE_cpotrf function\n"); }

        if ((iscolmajor && upper) || (!iscolmajor && !upper))
        {
            for (size_t r1=0u; r1<R; ++r1)
            {
                Y += 2u*r1 + 2u;
                for (size_t r2=r1+1u; r2<R; ++r2, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
            }
        }
        else
        {
            for (size_t r1=0u; r1<R; ++r1, Y+=2u*(R-r1+1u))
            {
                for (size_t r2=0u; r2<r1; ++r2, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
            }
        }
    }

    return 0;
}


int chol_z (double *Y, const double *X, const size_t R, const int iscolmajor, const int upper)
{
    const size_t N = R*R;

    if (N==0u) {}
    else if (N==1u) { *Y = *X; *++Y = *++X; }
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int Uplo = (upper) ? 'U' : 'L';
        
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
        Y -= 2u*N;
        
        if (LAPACKE_zpotrf(Ord,Uplo,(int)R,(_Complex double *)Y,(int)R))
        { fprintf(stderr,"error in chol_z: problem with LAPACKE_zpotrf function\n"); }

        if ((iscolmajor && upper) || (!iscolmajor && !upper))
        {
            for (size_t r1=0u; r1<R; ++r1)
            {
                Y += 2u*r1 + 2u;
                for (size_t r2=r1+1u; r2<R; ++r2, ++Y) { *Y = 0.0; *++Y = 0.0; }
            }
        }
        else
        {
            for (size_t r1=0u; r1<R; ++r1, Y+=2u*(R-r1+1u))
            {
                for (size_t r2=0u; r2<r1; ++r2, ++Y) { *Y = 0.0; *++Y = 0.0; }
            }
        }
    }

    return 0;
}


int chol_inplace_s (float *X, const size_t R, const int iscolmajor, const int upper)
{
    if (R<2) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int Uplo = (upper) ? 'U' : 'L';
        
        if (LAPACKE_spotrf(Ord,Uplo,(int)R,X,(int)R))
        { fprintf(stderr,"error in chol_inplace_s: problem with LAPACKE_spotrf function\n"); }

        if ((iscolmajor && upper) || (!iscolmajor && !upper))
        {
            for (size_t r1=0u; r1<R; ++r1)
            {
                X += r1 + 1u;
                for (size_t r2=r1+1u; r2<R; ++r2, ++X) { *X = 0.0f; }
            }
        }
        else
        {
            for (size_t r1=0u; r1<R; ++r1, X+=R-r1+1u)
            {
                for (size_t r2=0u; r2<r1; ++r2, ++X) { *X = 0.0f; }
            }
        }
    }

    return 0;
}


int chol_inplace_d (double *X, const size_t R, const int iscolmajor, const int upper)
{
    if (R<2) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int Uplo = (upper) ? 'U' : 'L';
        
        if (LAPACKE_dpotrf(Ord,Uplo,(int)R,X,(int)R))
        { fprintf(stderr,"error in chol_inplace_d: problem with LAPACKE_dpotrf function\n"); }

        if ((iscolmajor && upper) || (!iscolmajor && !upper))
        {
            for (size_t r1=0u; r1<R; ++r1)
            {
                X += r1 + 1u;
                for (size_t r2=r1+1u; r2<R; ++r2, ++X) { *X = 0.0; }
            }
        }
        else
        {
            for (size_t r1=0u; r1<R; ++r1, X+=R-r1+1u)
            {
                for (size_t r2=0u; r2<r1; ++r2, ++X) { *X = 0.0; }
            }
        }
    }

    return 0;
}


int chol_inplace_c (float *X, const size_t R, const int iscolmajor, const int upper)
{
    if (R<2) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int Uplo = (upper) ? 'U' : 'L';
        
        if (LAPACKE_cpotrf(Ord,Uplo,(int)R,(_Complex float *)X,(int)R))
        { fprintf(stderr,"error in chol_inplace_c: problem with LAPACKE_cpotrf function\n"); }

        if ((iscolmajor && upper) || (!iscolmajor && !upper))
        {
            for (size_t r1=0u; r1<R; ++r1)
            {
                X += 2u*r1 + 2u;
                for (size_t r2=r1+1u; r2<R; ++r2, ++X) { *X = 0.0f; *++X = 0.0f; }
            }
        }
        else
        {
            for (size_t r1=0u; r1<R; ++r1, X+=2u*(R-r1+1u))
            {
                for (size_t r2=0u; r2<r1; ++r2, ++X) { *X = 0.0f; *++X = 0.0f; }
            }
        }
    }

    return 0;
}


int chol_inplace_z (double *X, const size_t R, const int iscolmajor, const int upper)
{
    if (R<2) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const int Uplo = (upper) ? 'U' : 'L';
        
        if (LAPACKE_zpotrf(Ord,Uplo,(int)R,(_Complex double *)X,(int)R))
        { fprintf(stderr,"error in chol_inplace_z: problem with LAPACKE_zpotrf function\n"); }

        if ((iscolmajor && upper) || (!iscolmajor && !upper))
        {
            for (size_t r1=0u; r1<R; ++r1)
            {
                X += 2u*r1 + 2u;
                for (size_t r2=r1+1u; r2<R; ++r2, ++X) { *X = 0.0; *++X = 0.0; }
            }
        }
        else
        {
            for (size_t r1=0u; r1<R; ++r1, X+=2u*(R-r1+1u))
            {
                for (size_t r2=0u; r2<r1; ++r2, ++X) { *X = 0.0; *++X = 0.0; }
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
