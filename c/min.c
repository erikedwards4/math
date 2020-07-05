//Gets minimum of values for each row or col of X according to dim.
//For complex case, gets min absolute value (output is real).

//The inplace versions replace the first R or C elements of X,
//without regard to iscolmajor since output is a vector.

//Unfortunately, compiler doesn't find cblas_i?amin (but finds cblas_i?amax),
//so I do an inefficient solution for complex-valued case for now.


#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int min_s (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int min_d (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int min_c (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int min_z (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);

int min_inplace_s (float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int min_inplace_d (double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int min_inplace_c (float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int min_inplace_z (double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);


int min_s (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1) { cblas_scopy(N,X,1,Y,1); }
    else if (N1==N)
    {
        Y[0] = X[0];
        for (size_t n=1; n<N; n++) { if (X[n]<Y[0]) { Y[0] = X[n]; } }
    }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
            for (size_t n2=0, n=0; n2<N2; n2++)
            {
                Y[n2] = X[n]; n++;
                for (size_t n1=1; n1<N1; n1++, n++) { if (X[n]<Y[n2]) { Y[n2] = X[n]; } }
            }
            // for (size_t n2=0; n2<N2; n2++, Y++)
            // {
            //     *Y = *X; X++;
            //     for (size_t n1=1; n1<N1; n1++, X++) { if (*X<*Y) { *Y = *X; } }
            // }
            clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
        }
        else
        {
            for (size_t n2=0; n2<N2; n2++, Y++)
            {
                //Y[n2] = X[n2];  //so strange, this is slower here, but faster above and below
                //for (size_t n1=1; n1<N1; n1++) { if (X[n1*N2+n2]<Y[n2]) { Y[n2] = X[n1*N2+n2]; } }
                //float mn = X[n2];
                //for (size_t n1=1; n1<N1; n1++) { if (X[n1*N2+n2]<mn) { mn = X[n1*N2+n2]; } }
                *Y = X[n2];
                for (size_t n1=1; n1<N1; n1++) { if (X[n1*N2+n2]<*Y) { *Y = X[n1*N2+n2]; } }
            }
        }
    }
    else if (iscolmajor && dim==0)
    {
        for (size_t h=0, n=0, n2=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                for (size_t c=0; c<C; c++, n2++)
                {
                    Y[n2] = X[n];
                    for (size_t n1=0; n1<N1; n1++, n++) { if (X[n]<Y[n2]) { Y[n2] = X[n]; } }
                }
            }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0, n=0, n2=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, n2++)
            {
                Y[n2] = X[n];
                for (size_t n1=0; n1<N1; n1++, n++) { if (X[n]<Y[n2]) { Y[n2] = X[n]; } }
            }
        }
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, n+=J, n2++)
            {
                Y[n2] = X[n];
                for (size_t n1=0; n1<N1; n1+=K) { if (X[n+n1]<Y[n2]) { Y[n2] = X[n+n1]; } }
            }
        }
    }

    return 0;
}


int min_d (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1) { cblas_dcopy(N,X,1,Y,1); }
    else if (N1==N)
    {
        Y[0] = X[0];
        for (size_t n=1; n<N; n++) { if (X[n]<Y[0]) { Y[0] = X[n]; } }
    }
    else if (SH==1)
    {
        const size_t N2 = N/N1;
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0, n=0; n2<N2; n2++)
            {
                Y[n2] = X[n]; n++;
                for (size_t n1=1; n1<N1; n1++, n++) { if (X[n]<Y[n2]) { Y[n2] = X[n]; } }
            }
        }
        else
        {
            for (size_t n2=0; n2<N2; n2++)
            {
                double mn = X[n2];
                for (size_t n1=1; n1<N1; n1++) { if (X[n1*N2+n2]<mn) { mn = X[n1*N2+n2]; } }
                Y[n2] = mn;
            }
        }
    }
    else if (iscolmajor && dim==0)
    {
        for (size_t h=0, n=0, n2=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                for (size_t c=0; c<C; c++, n2++)
                {
                    Y[n2] = X[n];
                    for (size_t n1=0; n1<N1; n1++, n++) { if (X[n]<Y[n2]) { Y[n2] = X[n]; } }
                }
            }
        }
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        for (size_t r=0, n=0, n2=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, n2++)
            {
                Y[n2] = X[n];
                for (size_t n1=0; n1<N1; n1++, n++) { if (X[n]<Y[n2]) { Y[n2] = X[n]; } }
            }
        }
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, n+=J, n2++)
            {
                Y[n2] = X[n];
                for (size_t n1=0; n1<N1; n1+=K) { if (X[n+n1]<Y[n2]) { Y[n2] = X[n+n1]; } }
            }
        }
    }
    
    return 0;
}


// int min_c (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
// {
//     if (dim==0)
//     {
//         if (iscolmajor)
//         {
//             for (c=0; c<C; c++)
//             {
//                 r = (int)cblas_icamin((int)R,&X[2*c*R],1);
//                 Y[c] = sqrtf(X[2*(c*R+r)]*X[2*(c*R+r)] + X[2*(c*R+r)+1]*X[2*(c*R+r)+1]);
//             }
//         }
//         else
//         {
//             for (c=0; c<C; c++)
//             {
//                 r = (int)cblas_icamin((int)R,&X[2*c],(int)C);
//                 Y[c] = sqrtf(X[2*(r*C+c)]*X[2*(r*C+c)] + X[2*(r*C+c)+1]*X[2*(r*C+c)+1]);
//             }
//         }
//     }
//     else if (dim==1)
//     {
//         if (iscolmajor)
//         {
//             for (r=0; r<R; r++)
//             {
//                 c = (int)cblas_icamin((int)C,&X[2*r],(int)R);
//                 Y[r] = sqrtf(X[2*(c*R+r)]*X[2*(c*R+r)] + X[2*(c*R+r)+1]*X[2*(c*R+r)+1]);
//             }
//         }
//         else
//         {
//             for (r=0; r<R; r++)
//             {
//                 c = (int)cblas_icamin((int)C,&X[2*r*C],1);
//                 Y[r] = sqrtf(X[2*(r*C+c)]*X[2*(r*C+c)] + X[2*(r*C+c)+1]*X[2*(r*C+c)+1]);
//             }
//         }
//     }
//     else
//     {
//         fprintf(stderr,"error in min_c: dim must be 0 or 1.\n"); return 1;
//     }

//     return 0;
// }


// int min_z (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
// {
//     int r, c;

//     if (R<1) { fprintf(stderr,"error in min_z: R (nrows X) must be positive\n"); return 1; }
//     if (C<1) { fprintf(stderr,"error in min_z: C (ncols X) must be positive\n"); return 1; }

//     if (dim==0)
//     {
//         if (iscolmajor)
//         {
//             for (c=0; c<C; c++)
//             {
//                 r = (int)cblas_izamin((int)R,&X[2*c*R],1);
//                 Y[c] = sqrt(X[2*(c*R+r)]*X[2*(c*R+r)] + X[2*(c*R+r)+1]*X[2*(c*R+r)+1]);
//             }
//         }
//         else
//         {
//             for (c=0; c<C; c++)
//             {
//                 r = (int)cblas_izamin((int)R,&X[2*c],(int)C);
//                 Y[c] = sqrt(X[2*(r*C+c)]*X[2*(r*C+c)] + X[2*(r*C+c)+1]*X[2*(r*C+c)+1]);
//             }
//         }
//     }
//     else if (dim==1)
//     {
//         if (iscolmajor)
//         {
//             for (r=0; r<R; r++)
//             {
//                 c = (int)cblas_izamin((int)C,&X[2*r],(int)R);
//                 Y[r] = sqrt(X[2*(c*R+r)]*X[2*(c*R+r)] + X[2*(c*R+r)+1]*X[2*(c*R+r)+1]);
//             }
//         }
//         else
//         {
//             for (r=0; r<R; r++)
//             {
//                 c = (int)cblas_izamin((int)C,&X[2*r*C],1);
//                 Y[r] = sqrt(X[2*(r*C+c)]*X[2*(r*C+c)] + X[2*(r*C+c)+1]*X[2*(r*C+c)+1]);
//             }
//         }
//     }
//     else
//     {
//         fprintf(stderr,"error in min_z: dim must be 0 or 1.\n"); return 1;
//     }

//     return 0;
// }


int min_inplace_s (float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    int r, c;
    float mn;
    struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    
    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                r = (int)cblas_isamax((int)R,&X[c*R],1);
                if (X[c*R+r]>0.0f)
                {
                    mn = X[c*R+r];
                    for (r=0; r<R; r++) { X[c*R+r] -= mn; }
                    r = (int)cblas_isamax((int)R,&X[c*R],1);
                    X[c] = X[c*R+r] + mn;
                }
                else { X[c] = X[c*R+r]; }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                mn = X[c];
                for (r=1; r<R; r++) { if (X[r*C+c]<mn) { mn = X[r*C+c]; } }
                X[c] = mn;
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                mn = X[r];
                for (c=1; c<C; c++) { if (X[c*R+r]<mn) { mn = X[c*R+r]; } }
                X[r] = mn;
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                c = (int)cblas_isamax((int)C,&X[r*C],1);
                if (X[r*C+c]>0.0f)
                {
                    mn = X[r*C+c];
                    for (c=0; c<C; c++) { X[r*C+c] -= mn; }
                    c = (int)cblas_isamax((int)C,&X[r*C],1);
                    X[r] = X[r*C+c] + mn;
                }
                else { X[r] = X[r*C+c]; }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in min_inplace_s: dim must be 0 or 1.\n"); return 1;
    }
    clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int min_inplace_d (double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    int r, c;
    double mn;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                r = (int)cblas_idamax((int)R,&X[c*R],1);
                if (X[c*R+r]>0.0)
                {
                    mn = X[c*R+r];
                    for (r=0; r<R; r++) { X[c*R+r] -= mn; }
                    r = (int)cblas_idamax((int)R,&X[c*R],1);
                    X[c] = X[c*R+r] + mn;
                }
                else { X[c] = X[c*R+r]; }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                mn = X[c];
                for (r=1; r<R; r++) { if (X[r*C+c]<mn) { mn = X[r*C+c]; } }
                X[c] = mn;
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                mn = X[r];
                for (c=1; c<C; c++) { if (X[c*R+r]<mn) { mn = X[c*R+r]; } }
                X[r] = mn;
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                c = (int)cblas_idamax((int)C,&X[r*C],1);
                if (X[r*C+c]>0.0)
                {
                    mn = X[r*C+c];
                    for (c=0; c<C; c++) { X[r*C+c] -= mn; }
                    c = (int)cblas_idamax((int)C,&X[r*C],1);
                    X[r] = X[r*C+c] + mn;
                }
                else { X[r] = X[r*C+c]; }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in min_inplace_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int min_inplace_c (float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    float mn, a2;
    int r, c;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                mn = X[2*c*R]*X[2*c*R] + X[2*c*R+1]*X[2*c*R+1];
                for (r=1; r<R; r++)
                {
                    a2 = X[2*(c*R+r)]*X[2*(c*R+r)] + X[2*(c*R+r)+1]*X[2*(c*R+r)+1];
                    if (a2<mn) { mn = a2; }
                }
                X[c] = sqrtf(mn);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                mn = X[2*c]*X[2*c] + X[2*c+1]*X[2*c+1];
                for (r=1; r<R; r++)
                {
                    a2 = X[2*(r*C+c)]*X[2*(r*C+c)] + X[2*(r*C+c)+1]*X[2*(r*C+c)+1];
                    if (a2<mn) { mn = a2; }
                }
                X[c] = sqrtf(mn);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                mn = X[2*r]*X[2*r] + X[2*r+1]*X[2*r+1];
                for (c=1; c<C; c++)
                {
                    a2 = X[2*(c*R+r)]*X[2*(c*R+r)] + X[2*(c*R+r)+1]*X[2*(c*R+r)+1];
                    if (a2<mn) { mn = a2; }
                }
                X[r] = sqrtf(mn);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                mn = X[2*r*C]*X[2*r*C] + X[2*r*C+1]*X[2*r*C+1];
                for (c=1; c<C; c++)
                {
                    a2 = X[2*(r*C+c)]*X[2*(r*C+c)] + X[2*(r*C+c)+1]*X[2*(r*C+c)+1];
                    if (a2<mn) { mn = a2; }
                }
                X[r] = sqrtf(mn);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in min_inplace_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int min_inplace_z (double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    double mn, a2;
    int r, c;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                mn = X[2*c*R]*X[2*c*R] + X[2*c*R+1]*X[2*c*R+1];
                for (r=1; r<R; r++)
                {
                    a2 = X[2*(c*R+r)]*X[2*(c*R+r)] + X[2*(c*R+r)+1]*X[2*(c*R+r)+1];
                    if (a2<mn) { mn = a2; }
                }
                X[c] = sqrt(mn);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                mn = X[2*c]*X[2*c] + X[2*c+1]*X[2*c+1];
                for (r=1; r<R; r++)
                {
                    a2 = X[2*(r*C+c)]*X[2*(r*C+c)] + X[2*(r*C+c)+1]*X[2*(r*C+c)+1];
                    if (a2<mn) { mn = a2; }
                }
                X[c] = sqrt(mn);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                mn = X[2*r]*X[2*r] + X[2*r+1]*X[2*r+1];
                for (c=1; c<C; c++)
                {
                    a2 = X[2*(c*R+r)]*X[2*(c*R+r)] + X[2*(c*R+r)+1]*X[2*(c*R+r)+1];
                    if (a2<mn) { mn = a2; }
                }
                X[r] = sqrt(mn);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                mn = X[2*r*C]*X[2*r*C] + X[2*r*C+1]*X[2*r*C+1];
                for (c=1; c<C; c++)
                {
                    a2 = X[2*(r*C+c)]*X[2*(r*C+c)] + X[2*(r*C+c)+1]*X[2*(r*C+c)+1];
                    if (a2<mn) { mn = a2; }
                }
                X[r] = sqrt(mn);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in min_inplace_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
