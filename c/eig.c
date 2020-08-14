//Linear algebra function.
//Eigendecomposition of square, Hermitian-symmetric matrix X.
//Eigenvalues and vectors are returned in descending order (fliplr from initial order).

//Since only the first K are kept, the first K*R elements of U (or X for inplace) hold the result.
//The input matrices must be allocated to their full size!
//That is, when allocating, assume that K is at K_full = R.

//For complex case, eigenvalues in V are real-valued.

#include <stdio.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int eig_s (float *U, float *V, const float *X, const size_t R, const char iscolmajor, const size_t K);
int eig_d (double *U, double *V, const double *X, const size_t R, const char iscolmajor, const size_t K);
int eig_c (float *U, float *V, const float *X, const size_t R, const char iscolmajor, const size_t K);
int eig_z (double *U, double *V, const double *X, const size_t R, const char iscolmajor, const size_t K);

int eig_inplace_s (float *X, float *V, const size_t R, const char iscolmajor, const size_t K);
int eig_inplace_d (double *X, double *V, const size_t R, const char iscolmajor, const size_t K);
int eig_inplace_c (float *X, float *V, const size_t R, const char iscolmajor, const size_t K);
int eig_inplace_z (double *X, double *V, const size_t R, const char iscolmajor, const size_t K);


int eig_s (float *U, float *V, const float *X, const size_t R, const char iscolmajor, const size_t K)
{
    if (R==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        
        //Only lower half needs to be copied, but difference is negligible with much more code
        for (size_t n=0; n<R*R; ++n, ++X, ++U) { *U = *X; }
        U -= R*R;
        
        if (LAPACKE_ssyev(Ord,'V','L',(int)R,U,(int)R,V))
        { fprintf(stderr,"error in eig_s: problem with LAPACKE_ssyev function\n"); }

        //Flip to descending order
        float v1;
        for (size_t r=0; r<R/2; ++r, ++V) { v1 = *V; *V = *(V+R-2*r-1); *(V+R-2*r-1) = v1; }
        if (iscolmajor)
        {
            for (size_t c=0; c<R/2; ++c)
            {
                for (size_t r=0; r<R; ++r, ++U) { v1 = *U; *U = *(U+(R-2*c-1)*R); *(U+(R-2*c-1)*R) = v1; }
            }
        }
        else
        {
            for (size_t r=0; r<R; ++r, U+=R-R/2)
            {
                for (size_t c=0; c<R/2; ++c, ++U) { v1 = *U; *U = *(U+R-2*c-1); *(U+R-2*c-1) = v1; }
            }
            if (K>0 && K<R)
            {
                U -= R*R;
                for (size_t r=0; r<R; ++r)
                {
                    for (size_t k=0; k<K; ++k, ++U) { *U = *(U+r*(R-K)); }
                }
            }
        }
    }

    return 0;
}


int eig_d (double *U, double *V, const double *X, const size_t R, const char iscolmajor, const size_t K)
{
    if (R==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        
        for (size_t n=0; n<R*R; ++n, ++X, ++U) { *U = *X; }
        U -= R*R;
        
        if (LAPACKE_dsyev(Ord,'V','L',(int)R,U,(int)R,V))
        { fprintf(stderr,"error in eig_d: problem with LAPACKE_dsyev function\n"); }

        //Flip to descending order
        double v1;
        for (size_t r=0; r<R/2; ++r, ++V) { v1 = *V; *V = *(V+R-2*r-1); *(V+R-2*r-1) = v1; }
        if (iscolmajor)
        {
            for (size_t c=0; c<R/2; ++c)
            {
                for (size_t r=0; r<R; ++r, ++U) { v1 = *U; *U = *(U+(R-2*c-1)*R); *(U+(R-2*c-1)*R) = v1; }
            }
        }
        else
        {
            for (size_t r=0; r<R; ++r, U+=R-R/2)
            {
                for (size_t c=0; c<R/2; ++c, ++U) { v1 = *U; *U = *(U+R-2*c-1); *(U+R-2*c-1) = v1; }
            }
            if (K>0 && K<R)
            {
                U -= R*R;
                for (size_t r=0; r<R; ++r)
                {
                    for (size_t k=0; k<K; ++k, ++U) { *U = *(U+r*(R-K)); }
                }
            }
        }
    }

    return 0;
}


int eig_c (float *U, float *V, const float *X, const size_t R, const char iscolmajor, const size_t K)
{
    if (R==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        
        for (size_t n=0; n<2*R*R; ++n, ++X, ++U) { *U = *X; }
        U -= 2*R*R;
        
        if (LAPACKE_cheev(Ord,'V','L',(int)R,(_Complex float *)U,(int)R,V))
        { fprintf(stderr,"error in eig_c: problem with LAPACKE_csyev function\n"); }

        //Flip to descending order
        float v1;
        for (size_t r=0; r<R/2; ++r, ++V) { v1 = *V; *V = *(V+R-2*r-1); *(V+R-2*r-1) = v1; }
        if (iscolmajor)
        {
            for (size_t c=0; c<R/2; ++c)
            {
                for (size_t r=0; r<2*R; ++r, ++U) { v1 = *U; *U = *(U+2*(R-2*c-1)*R); *(U+2*(R-2*c-1)*R) = v1; }
            }
        }
        else
        {
            float v1r, v1i;
            for (size_t r=0; r<R; ++r, U+=2*(R-R/2))
            {
                for (size_t c=0; c<R/2; ++c)
                {
                    v1r = *U; v1i = *(U+1);
                    *U = *(U+2*(R-2*c)-2); ++U;
                    *U = *(U+2*(R-2*c)-2); ++U;
                    *(U+2*(R-2*c)-4) = v1r;
                    *(U+2*(R-2*c)-3) = v1i;
                }
            }
            if (K>0 && K<R)
            {
                U -= 2*R*R;
                for (size_t r=0; r<R; ++r)
                {
                    for (size_t k=0; k<2*K; ++k, ++U) { *U = *(U+2*r*(R-K)); }
                }
            }
        }
    }

    return 0;
}


int eig_z (double *U, double *V, const double *X, const size_t R, const char iscolmajor, const size_t K)
{
    if (R==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        
        for (size_t n=0; n<2*R*R; ++n, ++X, ++U) { *U = *X; }
        U -= 2*R*R;
        
        if (LAPACKE_zheev(Ord,'V','L',(int)R,(_Complex double *)U,(int)R,V))
        { fprintf(stderr,"error in eig_z: problem with LAPACKE_zsyev function\n"); }

        //Flip to descending order
        double v1;
        for (size_t r=0; r<R/2; ++r, ++V) { v1 = *V; *V = *(V+R-2*r-1); *(V+R-2*r-1) = v1; }
        if (iscolmajor)
        {
            for (size_t c=0; c<R/2; ++c)
            {
                for (size_t r=0; r<2*R; ++r, ++U) { v1 = *U; *U = *(U+2*(R-2*c-1)*R); *(U+2*(R-2*c-1)*R) = v1; }
            }
        }
        else
        {
            double v1r, v1i;
            for (size_t r=0; r<R; ++r, U+=2*(R-R/2))
            {
                for (size_t c=0; c<R/2; ++c)
                {
                    v1r = *U; v1i = *(U+1);
                    *U = *(U+2*(R-2*c)-2); ++U;
                    *U = *(U+2*(R-2*c)-2); ++U;
                    *(U+2*(R-2*c)-4) = v1r;
                    *(U+2*(R-2*c)-3) = v1i;
                }
            }
            if (K>0 && K<R)
            {
                U -= 2*R*R;
                for (size_t r=0; r<R; ++r)
                {
                    for (size_t k=0; k<2*K; ++k, ++U) { *U = *(U+2*r*(R-K)); }
                }
            }
        }
    }

    return 0;
}


int eig_inplace_s (float *X, float *V, const size_t R, const char iscolmajor, const size_t K)
{
    if (R==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        
        if (LAPACKE_ssyev(Ord,'V','L',(int)R,X,(int)R,V))
        { fprintf(stderr,"error in eig_inplace_s: problem with LAPACKE_ssyev function\n"); }

        //Flip to descending order
        float v1;
        for (size_t r=0; r<R/2; ++r, ++V) { v1 = *V; *V = *(V+R-2*r-1); *(V+R-2*r-1) = v1; }
        if (iscolmajor)
        {
            for (size_t c=0; c<R/2; ++c)
            {
                for (size_t r=0; r<R; ++r, ++X) { v1 = *X; *X = *(X+(R-2*c-1)*R); *(X+(R-2*c-1)*R) = v1; }
            }
        }
        else
        {
            for (size_t r=0; r<R; ++r, X+=R-R/2)
            {
                for (size_t c=0; c<R/2; ++c, ++X) { v1 = *X; *X = *(X+R-2*c-1); *(X+R-2*c-1) = v1; }
            }
            if (K>0 && K<R)
            {
                X -= R*R;
                for (size_t r=0; r<R; ++r)
                {
                    for (size_t k=0; k<K; ++k, ++X) { *X = *(X+r*(R-K)); }
                }
            }
        }
    }

    return 0;
}


int eig_inplace_d (double *X, double *V, const size_t R, const char iscolmajor, const size_t K)
{
    if (R==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        
        if (LAPACKE_dsyev(Ord,'V','L',(int)R,X,(int)R,V))
        { fprintf(stderr,"error in eig_inplace_d: problem with LAPACKE_dsyev function\n"); }

        //Flip to descending order
        double v1;
        for (size_t r=0; r<R/2; ++r, ++V) { v1 = *V; *V = *(V+R-2*r-1); *(V+R-2*r-1) = v1; }
        if (iscolmajor)
        {
            for (size_t c=0; c<R/2; ++c)
            {
                for (size_t r=0; r<R; ++r, ++X) { v1 = *X; *X = *(X+(R-2*c-1)*R); *(X+(R-2*c-1)*R) = v1; }
            }
        }
        else
        {
            for (size_t r=0; r<R; ++r, X+=R-R/2)
            {
                for (size_t c=0; c<R/2; ++c, ++X) { v1 = *X; *X = *(X+R-2*c-1); *(X+R-2*c-1) = v1; }
            }
            if (K>0 && K<R)
            {
                X -= R*R;
                for (size_t r=0; r<R; ++r)
                {
                    for (size_t k=0; k<K; ++k, ++X) { *X = *(X+r*(R-K)); }
                }
            }
        }
    }

    return 0;
}


int eig_inplace_c (float *X, float *V, const size_t R, const char iscolmajor, const size_t K)
{
    if (R==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        
        if (LAPACKE_cheev(Ord,'V','L',(int)R,(_Complex float *)X,(int)R,V))
        { fprintf(stderr,"error in eig_inplace_c: problem with LAPACKE_csyev function\n"); }

        //Flip to descending order
        float v1;
        for (size_t r=0; r<R/2; ++r, ++V) { v1 = *V; *V = *(V+R-2*r-1); *(V+R-2*r-1) = v1; }
        if (iscolmajor)
        {
            for (size_t c=0; c<R/2; ++c)
            {
                for (size_t r=0; r<2*R; ++r, ++X) { v1 = *X; *X = *(X+2*(R-2*c-1)*R); *(X+2*(R-2*c-1)*R) = v1; }
            }
        }
        else
        {
            float v1r, v1i;
            for (size_t r=0; r<R; ++r, X+=2*(R-R/2))
            {
                for (size_t c=0; c<R/2; ++c)
                {
                    v1r = *X; v1i = *(X+1);
                    *X = *(X+2*(R-2*c)-2); ++X;
                    *X = *(X+2*(R-2*c)-2); ++X;
                    *(X+2*(R-2*c)-4) = v1r;
                    *(X+2*(R-2*c)-3) = v1i;
                }
            }
            if (K>0 && K<R)
            {
                X -= 2*R*R;
                for (size_t r=0; r<R; ++r)
                {
                    for (size_t k=0; k<2*K; ++k, ++X) { *X = *(X+2*r*(R-K)); }
                }
            }
        }
    }

    return 0;
}


int eig_inplace_z (double *X, double *V, const size_t R, const char iscolmajor, const size_t K)
{
    if (R==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        
        if (LAPACKE_zheev(Ord,'V','L',(int)R,(_Complex double *)X,(int)R,V))
        { fprintf(stderr,"error in eig_inplace_z: problem with LAPACKE_zsyev function\n"); }

        //Flip to descending order
        double v1;
        for (size_t r=0; r<R/2; ++r, ++V) { v1 = *V; *V = *(V+R-2*r-1); *(V+R-2*r-1) = v1; }
        if (iscolmajor)
        {
            for (size_t c=0; c<R/2; ++c)
            {
                for (size_t r=0; r<2*R; ++r, ++X) { v1 = *X; *X = *(X+2*(R-2*c-1)*R); *(X+2*(R-2*c-1)*R) = v1; }
            }
        }
        else
        {
            double v1r, v1i;
            for (size_t r=0; r<R; ++r, X+=2*(R-R/2))
            {
                for (size_t c=0; c<R/2; ++c)
                {
                    v1r = *X; v1i = *(X+1);
                    *X = *(X+2*(R-2*c)-2); ++X;
                    *X = *(X+2*(R-2*c)-2); ++X;
                    *(X+2*(R-2*c)-4) = v1r;
                    *(X+2*(R-2*c)-3) = v1i;
                }
            }
            if (K>0 && K<R)
            {
                X -= 2*R*R;
                for (size_t r=0; r<R; ++r)
                {
                    for (size_t k=0; k<2*K; ++k, ++X) { *X = *(X+2*r*(R-K)); }
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
