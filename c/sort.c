//Vec2vec operation.
//Sorts each vector in X along dim.
//This has in-place and not-in-place versions.
//Too bad there is no LAPACKE_?laord for the in-place version.

#include <stdio.h>
#include <math.h>
#include <lapacke.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

// static int cmp_ascend_s (const void *a, const void *b);
// static int cmp_ascend_d (const void *a, const void *b);
static int cmp_ascend_c (const void *a, const void *b);
static int cmp_ascend_z (const void *a, const void *b);

// static int cmp_descend_s (const void *a, const void *b);
// static int cmp_descend_d (const void *a, const void *b);
static int cmp_descend_c (const void *a, const void *b);
static int cmp_descend_z (const void *a, const void *b);

int sort_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sort_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sort_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sort_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);

int sort_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sort_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sort_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sort_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);


// static int cmp_ascend_s (const void *a, const void *b)
// {
// 	const float x1 = *(const float*)a, x2 = *(const float*)b;
// 	if (x1!=x1) { return 1; }
//     else if (x2!=x2) { return -1; }
//     else if (x1>x2) { return 1; }
//     else if (x2>x1) { return -1; }
//     else { return 0; }
// }


// static int cmp_ascend_d (const void *a, const void *b)
// {
// 	const double x1 = *(const double*)a, x2 = *(const double*)b;
// 	if (x1!=x1) { return 1; }
//     else if (x2!=x2) { return -1; }
//     else if (x1>x2) { return 1; }
//     else if (x2>x1) { return -1; }
//     else { return 0; }
// }


static int cmp_ascend_c (const void *a, const void *b)
{
	const float sq1 = *(const float *)a**(const float *)a + *((const float *)(a)+1)**((const float *)(a)+1);
    const float sq2 = *(const float *)b**(const float *)b + *((const float *)(b)+1)**((const float *)(b)+1);
	if (sq1!=sq1) { return 1; }
    else if (sq2!=sq2) { return -1; }
    else if (sq1>sq2) { return 1; }
    else if (sq2>sq1) { return -1; }
	else if (atan2f(*((const float *)(a)+1),*(const float *)a) > atan2f(*((const float *)(b)+1),*(const float *)b)) { return 1; }
    else if (atan2f(*((const float *)(b)+1),*(const float *)b) > atan2f(*((const float *)(a)+1),*(const float *)a)) { return -1; }
    else { return 0; }
}


static int cmp_ascend_z (const void *a, const void *b)
{
	const double sq1 = *(const double *)a**(const double *)a + *((const double *)(a)+1)**((const double *)(a)+1);
    const double sq2 = *(const double *)b**(const double *)b + *((const double *)(b)+1)**((const double *)(b)+1);
	if (sq1!=sq1) { return 1; }
    else if (sq2!=sq2) { return -1; }
    else if (sq1>sq2) { return 1; }
    else if (sq2>sq1) { return -1; }
	else if (atan2(*((const double *)(a)+1),*(const double *)a) > atan2(*((const double *)(b)+1),*(const double *)b)) { return 1; }
    else if (atan2(*((const double *)(b)+1),*(const double *)b) > atan2(*((const double *)(a)+1),*(const double *)a)) { return -1; }
    else { return 0; }
}


// static int cmp_descend_s (const void *a, const void *b)
// {
// 	const float x1 = *(const float*)a, x2 = *(const float*)b;
// 	if (x1!=x1) { return -1; }
//     else if (x2!=x2) { return 1; }
//     else if (x1<x2) { return 1; }
//     else if (x2<x1) { return -1; }
//     else { return 0; }
// }


// static int cmp_descend_d (const void *a, const void *b)
// {
// 	const double x1 = *(const double*)a, x2 = *(const double*)b;
// 	if (x1!=x1) { return -1; }
//     else if (x2!=x2) { return 1; }
//     else if (x1<x2) { return 1; }
//     else if (x2<x1) { return -1; }
//     else { return 0; }
// }


static int cmp_descend_c (const void *a, const void *b)
{
	const float sq1 = *(const float *)a**(const float *)a + *((const float *)(a)+1)**((const float *)(a)+1);
    const float sq2 = *(const float *)b**(const float *)b + *((const float *)(b)+1)**((const float *)(b)+1);
	if (sq1!=sq1) { return -1; }
    else if (sq2!=sq2) { return 1; }
    else if (sq1<sq2) { return 1; }
    else if (sq2<sq1) { return -1; }
	else if (atan2f(*((const float *)(a)+1),*(const float *)a) < atan2f(*((const float *)(b)+1),*(const float *)b)) { return 1; }
    else if (atan2f(*((const float *)(b)+1),*(const float *)b) < atan2f(*((const float *)(a)+1),*(const float *)a)) { return -1; }
    else { return 0; }
}


static int cmp_descend_z (const void *a, const void *b)
{
	const double sq1 = *(const double *)a**(const double *)a + *((const double *)(a)+1)**((const double *)(a)+1);
    const double sq2 = *(const double *)b**(const double *)b + *((const double *)(b)+1)**((const double *)(b)+1);
	if (sq1!=sq1) { return -1; }
    else if (sq2!=sq2) { return 1; }
    else if (sq1<sq2) { return 1; }
    else if (sq2<sq1) { return -1; }
	else if (atan2(*((const double *)(a)+1),*(const double *)a) < atan2(*((const double *)(b)+1),*(const double *)b)) { return 1; }
    else if (atan2(*((const double *)(b)+1),*(const double *)b) < atan2(*((const double *)(a)+1),*(const double *)a)) { return -1; }
    else { return 0; }
}


int sort_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const char id = (ascend) ? 'I' : 'D';
    //int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_s : cmp_descend_s;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = *X; }
        Y -= L;
        if (LAPACKE_slasrt_work(id,(int)L,Y)) { fprintf(stderr,"error in sort_s: problem with LAPACKE function\n"); }
        //qsort(Y,L,sizeof(float),comp);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            //for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
            //Y -= N; //surprisingly, only same speed
            for (size_t v=0u; v<V; ++v, Y+=L)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = *X; }
                Y -= L;
                if (LAPACKE_slasrt_work(id,(int)L,Y)) { fprintf(stderr,"error in sort_s: problem with LAPACKE function\n"); }
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in sort_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, X1-=L, Y-=K*L-1u)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work(id,(int)L,X1)) { fprintf(stderr,"error in sort_s: problem with LAPACKE function\n"); }
                    for (size_t l=0u; l<L; ++l, ++X1, Y+=K) { *Y = *X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int sort_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const char id = (ascend) ? 'I' : 'D';

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = *X; }
        if (LAPACKE_dlasrt_work(id,(int)L,Y)) { fprintf(stderr,"error in sort_d: problem with LAPACKE function\n"); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, Y+=L)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = *X; }
                Y -= L;
                if (LAPACKE_dlasrt_work(id,(int)L,Y)) { fprintf(stderr,"error in sort_d: problem with LAPACKE function\n"); }
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in sort_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, X1-=L, Y-=K*L-1u)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work(id,(int)L,X1)) { fprintf(stderr,"error in sort_d: problem with LAPACKE function\n"); }
                    for (size_t l=0u; l<L; ++l, ++X1, Y+=K) { *Y = *X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int sort_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_c : cmp_descend_c;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<2u*L; ++l, ++X, ++Y) { *Y = *X; }
        qsort(Y,L,2*sizeof(float),comp);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, Y+=2u*L)
            {
                for (size_t l=0u; l<2u*L; ++l, ++X, ++Y) { *Y = *X; }
                Y -= 2u*L;
                qsort(Y,L,2*sizeof(float),comp);
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(2u*L*sizeof(float)))) { fprintf(stderr,"error in sort_c: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X1-=2u*L, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=0u; l<L; ++l, X+=2u*K-1u, ++X1) { *X1 = *X; *++X1 = *++X; }
                    X1 -= 2u*L;
                    qsort(X1,L,2*sizeof(float),comp);
                    X1 -= 2u*L;
                    for (size_t l=0u; l<L; ++l, ++X1, Y+=2u*K-1u) { *Y = *X1; *++Y = *++X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int sort_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_z : cmp_descend_z;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<2u*L; ++l, ++X, ++Y) { *Y = *X; }
        qsort(Y,L,2*sizeof(double),comp);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, Y+=2u*L)
            {
                for (size_t l=0u; l<2u*L; ++l, ++X, ++Y) { *Y = *X; }
                Y -= 2u*L;
                qsort(Y,L,2*sizeof(double),comp);
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(2u*L*sizeof(double)))) { fprintf(stderr,"error in sort_z: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X1-=2u*L, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=0u; l<L; ++l, X+=2u*K-1u, ++X1) { *X1 = *X; *++X1 = *++X; }
                    X1 -= 2u*L;
                    qsort(X1,L,2*sizeof(double),comp);
                    X1 -= 2u*L;
                    for (size_t l=0u; l<L; ++l, ++X1, Y+=2u*K-1u) { *Y = *X1; *++Y = *++X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int sort_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const char id = (ascend) ? 'I' : 'D';
    //int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_s : cmp_descend_s;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        if (LAPACKE_slasrt_work(id,(int)L,X)) { fprintf(stderr,"error in sort_inplace_s: problem with LAPACKE function\n"); }
        //qsort(X,L,sizeof(float),comp);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, X+=L)
            {
                if (LAPACKE_slasrt_work(id,(int)L,X)) { fprintf(stderr,"error in sort_inplace_s: problem with LAPACKE function\n"); }
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in sort_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, ++X1, X+=K+1u)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work(id,(int)L,X1)) { fprintf(stderr,"error in sort_inplace_s: problem with LAPACKE function\n"); }
                    X1 += L-1u; X -= K;
                    for (size_t l=0u; l<L; ++l, --X1, X-=K) { *X = *X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int sort_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const char id = (ascend) ? 'I' : 'D';

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        if (LAPACKE_dlasrt_work(id,(int)L,X)) { fprintf(stderr,"error in sort_inplace_d: problem with LAPACKE function\n"); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, X+=L)
            {
                if (LAPACKE_dlasrt_work(id,(int)L,X)) { fprintf(stderr,"error in sort_inplace_d: problem with LAPACKE function\n"); }
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in sort_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, ++X1, X+=K+1u)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work(id,(int)L,X1)) { fprintf(stderr,"error in sort_inplace_d: problem with LAPACKE function\n"); }
                    X1 += L-1u; X -= K;
                    for (size_t l=0u; l<L; ++l, --X1, X-=K) { *X = *X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int sort_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_inplace_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_c : cmp_descend_c;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        qsort(X,L,2*sizeof(float),comp);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, X+=2u*L)
            {
                qsort(X,L,2*sizeof(float),comp);
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(2u*L*sizeof(float)))) { fprintf(stderr,"error in sort_inplace_c: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X1+=2, X+=2u*K+2u)
                {
                    for (size_t l=0u; l<L; ++l, X+=2u*K-1u, ++X1) { *X1 = *X; *++X1 = *++X; }
                    X1 -= 2u*L;
                    qsort(X1,L,2*sizeof(float),comp);
                    X1 += 2u*L-2u; X -= 2u*K;
                    for (size_t l=0u; l<L; ++l, X1-=2, X-=2u*K) { *X = *X1; *(X+1) = *(X1+1); }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int sort_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_inplace_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_z : cmp_descend_z;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        qsort(X,L,2*sizeof(double),comp);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, X+=2u*L)
            {
                qsort(X,L,2*sizeof(double),comp);
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(2u*L*sizeof(double)))) { fprintf(stderr,"error in sort_inplace_z: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X1+=2, X+=2u*K+2u)
                {
                    for (size_t l=0u; l<L; ++l, X+=2u*K-1u, ++X1) { *X1 = *X; *++X1 = *++X; }
                    X1 -= 2u*L;
                    qsort(X1,L,2*sizeof(double),comp);
                    X1 += 2u*L-2u; X -= 2u*K;
                    for (size_t l=0u; l<L; ++l, X1-=2, X-=2u*K) { *X = *X1; *(X+1) = *(X1+1); }
                }
            }
            free(X1);
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
