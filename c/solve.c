//Linear algebra function.
//Least-squares solution of over- or under-determined linear system: A*X = B.
//The system can be over- or under-determined.

//This is similar to the \\ (backslash) operator in Octave: X = A\B.

//The inplace version modifies A and B during processing, with X found in B.
//The not-inplace version leaves them constant, but requires tmp copies.

#include <stdio.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int solve_s (float *X, const float *A, const float *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr);
int solve_d (double *X, const double *A, const double *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr);
int solve_c (float *X, const float *A, const float *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr);
int solve_z (double *X, const double *A, const double *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr);

int solve_inplace_s (float *A, float *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr);
int solve_inplace_d (double *A, double *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr);
int solve_inplace_c (float *A, float *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr);
int solve_inplace_z (double *A, double *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr);


int solve_s (float *X, const float *A, const float *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr)
{
    const size_t R2 = (tr) ? C1 : R1;
    const size_t R = (tr) ? R1 : C1, C = C2;
    const size_t N1 = R1*C1, N2 = R2*C2;

    if (R*C==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const char Trans = (tr) ? 'T' : 'N';
        const int lda = (iscolmajor) ? (int)R1 : (int)C1;
        const int ldb = (iscolmajor) ? (int)R2 : (int)C2;

        //tmp copies
        float *Atmp, *Btmp;
        if (!(Atmp=(float *)malloc(N1*sizeof(float)))) { fprintf(stderr,"error in solve_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Btmp=(float *)malloc(N2*sizeof(float)))) { fprintf(stderr,"error in solve_s: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t n=0u; n<N1; ++n, ++A, ++Atmp) { *Atmp = *A; }
        for (size_t n=0u; n<N2; ++n, ++B, ++Btmp) { *Btmp = *B; }
        Atmp -= N1; Btmp -= N2;
        
        //Solve
        if (LAPACKE_sgels(Ord,Trans,(int)R1,(int)C1,(int)C2,Atmp,lda,Btmp,ldb))
        { fprintf(stderr,"error in solve_s: problem with LAPACKE_sgels function\n"); }

        //Copy B to output X
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c, Btmp+=R2-R)
            {
                for (size_t r=0u; r<R; ++r, ++Btmp, ++X) { *X = *Btmp; }
            }
            Btmp -= R2*C;
        }
        else
        {
            for (size_t r=0u; r<R; ++r, Btmp+=C2-C)
            {
                for (size_t c=0u; c<C; ++c, ++Btmp, ++X) { *X = *Btmp; }
            }
            Btmp -= R*C2;
        }

        //Free
        free(Atmp); free(Btmp);
    }

    return 0;
}


int solve_d (double *X, const double *A, const double *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr)
{
    const size_t R2 = (tr) ? C1 : R1;
    const size_t R = (tr) ? R1 : C1, C = C2;
    const size_t N1 = R1*C1, N2 = R2*C2;

    if (R*C==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const char Trans = (tr) ? 'T' : 'N';
        const int lda = (iscolmajor) ? (int)R1 : (int)C1;
        const int ldb = (iscolmajor) ? (int)R2 : (int)C2;

        //tmp copies
        double *Atmp, *Btmp;
        if (!(Atmp=(double *)malloc(N1*sizeof(double)))) { fprintf(stderr,"error in solve_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Btmp=(double *)malloc(N2*sizeof(double)))) { fprintf(stderr,"error in solve_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t n=0u; n<N1; ++n, ++A, ++Atmp) { *Atmp = *A; }
        for (size_t n=0u; n<N2; ++n, ++B, ++Btmp) { *Btmp = *B; }
        Atmp -= N1; Btmp -= N2;
        
        //Solve
        if (LAPACKE_dgels(Ord,Trans,(int)R1,(int)C1,(int)C2,Atmp,lda,Btmp,ldb))
        { fprintf(stderr,"error in solve_d: problem with LAPACKE_dgels function\n"); }

        //Copy B to output X
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c, Btmp+=R2-R)
            {
                for (size_t r=0u; r<R; ++r, ++Btmp, ++X) { *X = *Btmp; }
            }
            Btmp -= R2*C;
        }
        else
        {
            for (size_t r=0u; r<R; ++r, Btmp+=C2-C)
            {
                for (size_t c=0u; c<C; ++c, ++Btmp, ++X) { *X = *Btmp; }
            }
            Btmp -= R*C2;
        }

        //Free
        free(Atmp); free(Btmp);
    }

    return 0;
}


int solve_c (float *X, const float *A, const float *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr)
{
    const size_t R2 = (tr) ? C1 : R1;
    const size_t R = (tr) ? R1 : C1, C = C2;
    const size_t N1 = R1*C1, N2 = R2*C2;

    if (R*C==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const char Trans = (tr) ? 'T' : 'N';
        const int lda = (iscolmajor) ? (int)R1 : (int)C1;
        const int ldb = (iscolmajor) ? (int)R2 : (int)C2;

        //tmp copies
        float *Atmp, *Btmp;
        if (!(Atmp=(float *)malloc(2*N1*sizeof(float)))) { fprintf(stderr,"error in solve_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Btmp=(float *)malloc(2*N2*sizeof(float)))) { fprintf(stderr,"error in solve_c: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t n=0u; n<2*N1; ++n, ++A, ++Atmp) { *Atmp = *A; }
        for (size_t n=0u; n<2*N2; ++n, ++B, ++Btmp) { *Btmp = *B; }
        Atmp -= 2*N1; Btmp -= 2*N2;
        
        //Solve
        if (LAPACKE_cgels(Ord,Trans,(int)R1,(int)C1,(int)C2,(_Complex float *)Atmp,lda,(_Complex float *)Btmp,ldb))
        { fprintf(stderr,"error in solve_c: problem with LAPACKE_cgels function\n"); }

        //Copy B to output X
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c, Btmp+=2*(R2-R))
            {
                for (size_t r=0u; r<2*R; ++r, ++Btmp, ++X) { *X = *Btmp; }
            }
            Btmp -= 2*R2*C;
        }
        else
        {
            for (size_t r=0u; r<R; ++r, Btmp+=2*(C2-C))
            {
                for (size_t c=0u; c<2*C; ++c, ++Btmp, ++X) { *X = *Btmp; }
            }
            Btmp -= 2*R*C2;
        }

        //Free
        free(Atmp); free(Btmp);
    }

    return 0;
}


int solve_z (double *X, const double *A, const double *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr)
{
    const size_t R2 = (tr) ? C1 : R1;
    const size_t R = (tr) ? R1 : C1, C = C2;
    const size_t N1 = R1*C1, N2 = R2*C2;

    if (R*C==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const char Trans = (tr) ? 'T' : 'N';
        const int lda = (iscolmajor) ? (int)R1 : (int)C1;
        const int ldb = (iscolmajor) ? (int)R2 : (int)C2;

        //tmp copies
        double *Atmp, *Btmp;
        if (!(Atmp=(double *)malloc(2*N1*sizeof(double)))) { fprintf(stderr,"error in solve_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Btmp=(double *)malloc(2*N2*sizeof(double)))) { fprintf(stderr,"error in solve_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t n=0u; n<2*N1; ++n, ++A, ++Atmp) { *Atmp = *A; }
        for (size_t n=0u; n<2*N2; ++n, ++B, ++Btmp) { *Btmp = *B; }
        Atmp -= 2*N1; Btmp -= 2*N2;
        
        //Solve
        if (LAPACKE_zgels(Ord,Trans,(int)R1,(int)C1,(int)C2,(_Complex double *)Atmp,lda,(_Complex double *)Btmp,ldb))
        { fprintf(stderr,"error in solve_z: problem with LAPACKE_zgels function\n"); }

        //Copy B to output X
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c, Btmp+=2*(R2-R))
            {
                for (size_t r=0u; r<2*R; ++r, ++Btmp, ++X) { *X = *Btmp; }
            }
            Btmp -= 2*R2*C;
        }
        else
        {
            for (size_t r=0u; r<R; ++r, Btmp+=2*(C2-C))
            {
                for (size_t c=0u; c<2*C; ++c, ++Btmp, ++X) { *X = *Btmp; }
            }
            Btmp -= 2*R*C2;
        }

        //Free
        free(Atmp); free(Btmp);
    }

    return 0;
}


int solve_inplace_s (float *A, float *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr)
{
    const size_t R2 = (tr) ? C1 : R1;
    const size_t R = (tr) ? R1 : C1, C = C2;

    if (R*C==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const char Trans = (tr) ? 'T' : 'N';
        const int lda = (iscolmajor) ? (int)R1 : (int)C1;
        const int ldb = (iscolmajor) ? (int)R2 : (int)C2;

        //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
        
        if (LAPACKE_sgels(Ord,Trans,(int)R1,(int)C1,(int)C2,A,lda,B,ldb))
        { fprintf(stderr,"error in solve_inplace_s: problem with LAPACKE_sgels function\n"); }

        //This works, but is not any faster (skips Nan check, but not significant)
        // float work_query, *work;
        // int lwork = -1;
        // if (LAPACKE_sgels_work(Ord,Trans,(int)R1,(int)C1,(int)C2,A,lda,B,ldb,&work_query,lwork))
        // { fprintf(stderr,"error in solve_inplace_s: problem with LAPACKE_sgels_work function\n"); }
        // lwork = (int)work_query;
        // if (!(work=(float *)malloc((size_t)lwork*sizeof(float)))) { fprintf(stderr,"error in solve_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        // if (LAPACKE_sgels_work(Ord,Trans,(int)R1,(int)C1,(int)C2,A,lda,B,ldb,work,lwork))
        // { fprintf(stderr,"error in solve_inplace_s: problem with LAPACKE_sgels_work function\n"); }
        // free(work);

        //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                for (size_t r=0u; r<R; ++r, ++B) { *B = *(B+c*(R2-R)); }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t c=0u; c<C; ++c, ++B) { *B = *(B+r*(C2-C)); }
            }
        }
    }

    return 0;
}


int solve_inplace_d (double *A, double *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr)
{
    const size_t R2 = (tr) ? C1 : R1;
    const size_t R = (tr) ? R1 : C1, C = C2;

    if (R*C==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const char Trans = (tr) ? 'T' : 'N';
        const int lda = (iscolmajor) ? (int)R1 : (int)C1;
        const int ldb = (iscolmajor) ? (int)R2 : (int)C2;
        
        if (LAPACKE_dgels(Ord,Trans,(int)R1,(int)C1,(int)C2,A,lda,B,ldb))
        { fprintf(stderr,"error in solve_inplace_d: problem with LAPACKE_dgels function\n"); }

        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                for (size_t r=0u; r<R; ++r, ++B) { *B = *(B+c*(R2-R)); }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t c=0u; c<C; ++c, ++B) { *B = *(B+r*(C2-C)); }
            }
        }
    }

    return 0;
}


int solve_inplace_c (float *A, float *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr)
{
    const size_t R2 = (tr) ? C1 : R1;
    const size_t R = (tr) ? R1 : C1, C = C2;

    if (R*C==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const char Trans = (tr) ? 'T' : 'N';
        const int lda = (iscolmajor) ? (int)R1 : (int)C1;
        const int ldb = (iscolmajor) ? (int)R2 : (int)C2;
        
        if (LAPACKE_cgels(Ord,Trans,(int)R1,(int)C1,(int)C2,(_Complex float *)A,lda,(_Complex float *)B,ldb))
        { fprintf(stderr,"error in solve_inplace_c: problem with LAPACKE_cgels function\n"); }

        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                for (size_t r=0u; r<2*R; ++r, ++B) { *B = *(B+2*c*(R2-R)); }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t c=0u; c<2*C; ++c, ++B) { *B = *(B+2*r*(C2-C)); }
            }
        }
    }

    return 0;
}


int solve_inplace_z (double *A, double *B, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor, const char tr)
{
    const size_t R2 = (tr) ? C1 : R1;
    const size_t R = (tr) ? R1 : C1, C = C2;

    if (R*C==0) {}
    else
    {
        const int Ord = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        const char Trans = (tr) ? 'T' : 'N';
        const int lda = (iscolmajor) ? (int)R1 : (int)C1;
        const int ldb = (iscolmajor) ? (int)R2 : (int)C2;
        
        if (LAPACKE_zgels(Ord,Trans,(int)R1,(int)C1,(int)C2,(_Complex double *)A,lda,(_Complex double *)B,ldb))
        { fprintf(stderr,"error in solve_inplace_z: problem with LAPACKE_zgels function\n"); }

        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                for (size_t r=0u; r<2*R; ++r, ++B) { *B = *(B+2*c*(R2-R)); }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t c=0u; c<2*C; ++c, ++B) { *B = *(B+2*r*(C2-C)); }
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
