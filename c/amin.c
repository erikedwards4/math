//Gets maximum of absolute values for each row or col of X according to dim.
//This is also the -Inf-norm (or min-norm) of each vector in X.

//For complex case, finds max absolute value and outputs the complex number.
//For complex case, this is |Xr|+|Xi|; see max for the usual sqrt(Xr*Xr+Xl*Xi).


#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int amin_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int amin_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int amin_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int amin_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int amin_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in amin_s: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (L==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        size_t l = cblas_isamin((int)N,X,1);
        *Y = *(X+l);
    }
    else if (SH==1)
    {
        const size_t V = N/L;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t v=0; v<V; v++, Y++)
            {
                size_t l = cblas_isamin((int)L,X,1);
                X += l; *Y = *X; X += L-l;
            }
        }
        else
        {
            float *mns;
            if (!(mns=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in amin_s: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0; l<L; l++, Y-=V, mns-=V)
            {
                for (size_t v=0; v<V; v++, X++, Y++, mns++)
                {
                    if (*X**X<*mns) { *Y = *X; *mns = *X**X; }
                }
            }
            free(mns);
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y++)
            {
                size_t l = cblas_isamin((int)L,X,1);
                X += l; *Y = *X; X += L-l;
            }
        }
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, Y++)
            {
                size_t l = cblas_isamin((int)L,X,(int)K);
                X += l*K; *Y = *X;
                X -= (int)((L-l)*K)-(int)J;
            }
        }
    }

    return 0;
}


int amin_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in amin_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (L==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        size_t l = cblas_idamin((int)N,X,1);
        *Y = *(X+l);
    }
    else if (SH==1)
    {
        const size_t V = N/L;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t v=0; v<V; v++, Y++)
            {
                size_t l = cblas_idamin((int)L,X,1);
                X += l; *Y = *X; X += L-l;
            }
        }
        else
        {
            double *mns;
            if (!(mns=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in amin_d: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0; l<L; l++, Y-=V, mns-=V)
            {
                for (size_t v=0; v<V; v++, X++, Y++, mns++)
                {
                    if (*X**X<*mns) { *Y = *X; *mns = *X**X; }
                }
            }
            free(mns);
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y++)
            {
                size_t l = cblas_idamin((int)L,X,1);
                X += l; *Y = *X; X += L-l;
            }
        }
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, Y++)
            {
                size_t l = cblas_idamin((int)L,X,(int)K);
                X += l*K; *Y = *X;
                X -= (int)((L-l)*K)-(int)J;
            }
        }
    }

    return 0;
}


int amin_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in amin_c: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (L==1) { cblas_ccopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        size_t l = cblas_icamin((int)N,X,1);
        X += 2*l; *Y = *X; *(Y+1) = *(X+1);
    }
    else if (SH==1)
    {
        const size_t V = N/L;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t v=0; v<V; v++, Y+=2)
            {
                size_t l = cblas_icamin((int)L,X,1);
                X += 2*l; *Y = *X; *(Y+1) = *(X+1);
                X += 2*(L-l);
            }
        }
        else
        {
            float *mns;
            if (!(mns=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in amin_c: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0; l<L; l++, Y-=2*V, mns-=V)
            {
                for (size_t v=0; v<V; v++, X+=2, Y+=2, mns++)
                {
                    if (sqrt(*X**X)+sqrt(*(X+1)**(X+1))<*mns)
                    {
                        *Y = *X; *(Y+1) = *(X+1);
                        *mns = sqrt(*X**X)+sqrt(*(X+1)**(X+1));
                    }
                }
            }
            free(mns);
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y+=2)
            {
                size_t l = cblas_icamin((int)L,X,1);
                X += 2*l; *Y = *X; *(Y+1) = *(X+1);
                X += 2*(L-l);
            }
        }
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t g=0; g<G; g++, X+=2*B*(L-J))
        {
            for (size_t b=0; b<B; b++, Y+=2)
            {
                size_t l = cblas_icamin((int)L,X,(int)K);
                X += 2*l*K; *Y = *X; *(Y+1) = *(X+1);
                X -= 2*((int)((L-l)*K)-(int)J);
            }
        }
    }

    return 0;
}


int amin_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in amin_z: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (L==1) { cblas_zcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        size_t l = cblas_izamin((int)N,X,1);
        X += 2*l; *Y = *X; *(Y+1) = *(X+1);
    }
    else if (SH==1)
    {
        const size_t V = N/L;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t v=0; v<V; v++, Y+=2)
            {
                size_t l = cblas_izamin((int)L,X,1);
                X += 2*l; *Y = *X; *(Y+1) = *(X+1);
                X += 2*(L-l);
            }
        }
        else
        {
            double *mns;
            if (!(mns=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in amin_z: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0; l<L; l++, Y-=2*V, mns-=V)
            {
                for (size_t v=0; v<V; v++, X+=2, Y+=2, mns++)
                {
                    if (sqrt(*X**X)+sqrt(*(X+1)**(X+1))<*mns)
                    {
                        *Y = *X; *(Y+1) = *(X+1);
                        *mns = sqrt(*X**X)+sqrt(*(X+1)**(X+1));
                    }
                }
            }
            free(mns);
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, Y+=2)
            {
                size_t l = cblas_izamin((int)L,X,1);
                X += 2*l; *Y = *X; *(Y+1) = *(X+1);
                X += 2*(L-l);
            }
        }
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t g=0; g<G; g++, X+=2*B*(L-J))
        {
            for (size_t b=0; b<B; b++, Y+=2)
            {
                size_t l = cblas_izamin((int)L,X,(int)K);
                X += 2*l*K; *Y = *X; *(Y+1) = *(X+1);
                X -= 2*((int)((L-l)*K)-(int)J);
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
