//2-input elementwise function.
//Does elementwise multiplication: Y = X1*X2
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int times_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int times_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int times_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int times_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);

int times_inplace_s (float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int times_inplace_d (double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int times_inplace_c (float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int times_inplace_z (double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);


int times_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    
    if (N1==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X1[0] * X2[n]; }
    }    
    else if (N2==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X1[n] * X2[0]; }
    }
    else if (N1==N2)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X1[n] * X2[n]; }
    }
    else if (iscolmajor)
    {
        const int r1i = (int)(R1>1), r2i = (int)(R2>1);
        const int c1i = (int)R1*((int)(C1>1)-(int)(R1>1)), c2i = (int)R2*((int)(C2>1)-(int)(R2>1));
        const int s1i = (int)(R1*C1)*((int)(S1>1)-(int)(C1>1)), s2i = (int)(R2*C2)*((int)(S2>1)-(int)(C2>1));
        const int h1i = (int)(R1*C1*S1)*((int)(H1>1)-(int)(S1>1)), h2i = (int)(R2*C2*S2)*((int)(H2>1)-(int)(S2>1));
        for (size_t h=0; h<H; h++, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0; s<S; s++, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0; c<C; c++, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0; r<R; r++, X1+=r1i, X2+=r2i)
                    {
                        *Y++ = *X1 * *X2;
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = (int)(H1>1), h2i = (int)(H2>1);
        const int s1i = (int)H1*((int)(S1>1)-(int)(H1>1)), s2i = (int)H2*((int)(S2>1)-(int)(H2>1));
        const int c1i = (int)(H1*S1)*((int)(C1>1)-(int)(S1>1)), c2i = (int)(H2*S2)*((int)(C2>1)-(int)(S2>1));
        const int r1i = (int)(H1*S1*C1)*((int)(R1>1)-(int)(C1>1)), r2i = (int)(H2*S2*C2)*((int)(R2>1)-(int)(C2>1));
        for (size_t r=0; r<R; r++, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0; c<C; c++, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0; s<S; s++, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0; h<H; h++, X1+=h1i, X2+=h2i)
                    {
                        *Y++ = *X1 * *X2;
                    }
                }
            }
        }
    }

    return 0;
}


int times_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    
    if (N1==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X1[0] * X2[n]; }
    }    
    else if (N2==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X1[n] * X2[0]; }
    }
    else if (N1==N2)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X1[n] * X2[n]; }
    }
    else if (iscolmajor)
    {
        const int r1i = (int)(R1>1), r2i = (int)(R2>1);
        const int c1i = (int)R1*((int)(C1>1)-(int)(R1>1)), c2i = (int)R2*((int)(C2>1)-(int)(R2>1));
        const int s1i = (int)(R1*C1)*((int)(S1>1)-(int)(C1>1)), s2i = (int)(R2*C2)*((int)(S2>1)-(int)(C2>1));
        const int h1i = (int)(R1*C1*S1)*((int)(H1>1)-(int)(S1>1)), h2i = (int)(R2*C2*S2)*((int)(H2>1)-(int)(S2>1));
        for (size_t h=0; h<H; h++, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0; s<S; s++, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0; c<C; c++, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0; r<R; r++, X1+=r1i, X2+=r2i)
                    {
                        *Y++ = *X1 * *X2;
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = (int)(H1>1), h2i = (int)(H2>1);
        const int s1i = (int)H1*((int)(S1>1)-(int)(H1>1)), s2i = (int)H2*((int)(S2>1)-(int)(H2>1));
        const int c1i = (int)(H1*S1)*((int)(C1>1)-(int)(S1>1)), c2i = (int)(H2*S2)*((int)(C2>1)-(int)(S2>1));
        const int r1i = (int)(H1*S1*C1)*((int)(R1>1)-(int)(C1>1)), r2i = (int)(H2*S2*C2)*((int)(R2>1)-(int)(C2>1));
        for (size_t r=0; r<R; r++, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0; c<C; c++, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0; s<S; s++, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0; h<H; h++, X1+=h1i, X2+=h2i)
                    {
                        *Y++ = *X1 * *X2;
                    }
                }
            }
        }
    }

    return 0;
}


int times_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;

    if (N1==1)
    {
        for (size_t n2=0; n2<2*N; n2+=2) { Y[n2] = X1[0]+X2[n2]; Y[n2+1] = X1[1]+X2[n2+1]; }
    }    
    else if (N2==1)
    {
        for (size_t n2=0; n2<2*N; n2+=2) { Y[n2] = X1[n2]+X2[0]; Y[n2+1] = X1[n2+1]+X2[0]; }
    }
    else if (N1==N2)
    {
        for (size_t n=0; n<2*N; n++) { Y[n] = X1[n] * X2[n]; }
    }
    else if (iscolmajor)
    {
        const int r1i = 2*(int)(R1>1), r2i = 2*(int)(R2>1);
        const int c1i = 2*(int)R1*((int)(C1>1)-(int)(R1>1)), c2i = 2*(int)R2*((int)(C2>1)-(int)(R2>1));
        const int s1i = 2*(int)(R1*C1)*((int)(S1>1)-(int)(C1>1)), s2i = 2*(int)(R2*C2)*((int)(S2>1)-(int)(C2>1));
        const int h1i = 2*(int)(R1*C1*S1)*((int)(H1>1)-(int)(S1>1)), h2i = 2*(int)(R2*C2*S2)*((int)(H2>1)-(int)(S2>1));
        for (size_t h=0; h<H; h++, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0; s<S; s++, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0; c<C; c++, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0; r<R; r++, X1+=r1i, X2+=r2i)
                    {
                        *Y++ = *X1 * *X2;
                        *Y++ = *(X1+1) * *(X2+1);
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = 2*(int)(H1>1), h2i = 2*(int)(H2>1);
        const int s1i = 2*(int)H1*((int)(S1>1)-(int)(H1>1)), s2i = 2*(int)H2*((int)(S2>1)-(int)(H2>1));
        const int c1i = 2*(int)(H1*S1)*((int)(C1>1)-(int)(S1>1)), c2i = 2*(int)(H2*S2)*((int)(C2>1)-(int)(S2>1));
        const int r1i = 2*(int)(H1*S1*C1)*((int)(R1>1)-(int)(C1>1)), r2i = 2*(int)(H2*S2*C2)*((int)(R2>1)-(int)(C2>1));
        for (size_t r=0; r<R; r++, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0; c<C; c++, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0; s<S; s++, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0; h<H; h++, X1+=h1i, X2+=h2i)
                    {
                        *Y++ = *X1 * *X2;
                        *Y++ = *(X1+1) * *(X2+1);
                    }
                }
            }
        }
    }

    return 0;
}


int times_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;

    if (N1==1)
    {
        for (size_t n2=0; n2<2*N; n2+=2) { Y[n2] = X1[0]+X2[n2]; Y[n2+1] = X1[1]+X2[n2+1]; }
    }    
    else if (N2==1)
    {
        for (size_t n2=0; n2<2*N; n2+=2) { Y[n2] = X1[n2]+X2[0]; Y[n2+1] = X1[n2+1]+X2[0]; }
    }
    else if (N1==N2)
    {
        for (size_t n=0; n<2*N; n++) { Y[n] = X1[n] * X2[n]; }
    }
    else if (iscolmajor)
    {
        const int r1i = 2*(int)(R1>1), r2i = 2*(int)(R2>1);
        const int c1i = 2*(int)R1*((int)(C1>1)-(int)(R1>1)), c2i = 2*(int)R2*((int)(C2>1)-(int)(R2>1));
        const int s1i = 2*(int)(R1*C1)*((int)(S1>1)-(int)(C1>1)), s2i = 2*(int)(R2*C2)*((int)(S2>1)-(int)(C2>1));
        const int h1i = 2*(int)(R1*C1*S1)*((int)(H1>1)-(int)(S1>1)), h2i = 2*(int)(R2*C2*S2)*((int)(H2>1)-(int)(S2>1));
        for (size_t h=0; h<H; h++, X1+=h1i, X2+=h2i)
        {
            for (size_t s=0; s<S; s++, X1+=s1i, X2+=s2i)
            {
                for (size_t c=0; c<C; c++, X1+=c1i, X2+=c2i)
                {
                    for (size_t r=0; r<R; r++, X1+=r1i, X2+=r2i)
                    {
                        *Y++ = *X1 * *X2;
                        *Y++ = *(X1+1) * *(X2+1);
                    }
                }
            }
        }
    }
    else
    {
        const int h1i = 2*(int)(H1>1), h2i = 2*(int)(H2>1);
        const int s1i = 2*(int)H1*((int)(S1>1)-(int)(H1>1)), s2i = 2*(int)H2*((int)(S2>1)-(int)(H2>1));
        const int c1i = 2*(int)(H1*S1)*((int)(C1>1)-(int)(S1>1)), c2i = 2*(int)(H2*S2)*((int)(C2>1)-(int)(S2>1));
        const int r1i = 2*(int)(H1*S1*C1)*((int)(R1>1)-(int)(C1>1)), r2i = 2*(int)(H2*S2*C2)*((int)(R2>1)-(int)(C2>1));
        for (size_t r=0; r<R; r++, X1+=r1i, X2+=r2i)
        {
            for (size_t c=0; c<C; c++, X1+=c1i, X2+=c2i)
            {
                for (size_t s=0; s<S; s++, X1+=s1i, X2+=s2i)
                {
                    for (size_t h=0; h<H; h++, X1+=h1i, X2+=h2i)
                    {
                        *Y++ = *X1 * *X2;
                        *Y++ = *(X1+1) * *(X2+1);
                    }
                }
            }
        }
    }

    return 0;
}


int times_inplace_s (float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    if (N1!=N) { fprintf(stderr,"error in times_inplace_s: first input (X1) cannot be broadcast for inplace version\n"); return 1; }
    
    if (N2==1)
    {
        for (size_t n=0; n<N; n++) { X1[n] *= X2[0]; }
    }
    else if (N==N2)
    {
        for (size_t n=0; n<N; n++) { X1[n] *= X2[n]; }
    }
    else if (iscolmajor)
    {
        const int r2i = (int)(R2>1);
        const int c2i = (int)R2*((int)(C2>1)-(int)(R2>1));
        const int s2i = (int)(R2*C2)*((int)(S2>1)-(int)(C2>1));
        const int h2i = (int)(R2*C2*S2)*((int)(H2>1)-(int)(S2>1));
        for (size_t h=0; h<H; h++, X2+=h2i)
        {
            for (size_t s=0; s<S; s++, X2+=s2i)
            {
                for (size_t c=0; c<C; c++, X2+=c2i)
                {
                    for (size_t r=0; r<R; r++, X2+=r2i)
                    {
                        *X1++ *= *X2;
                    }
                }
            }
        }
    }
    else
    {
        const int h2i = (int)(H2>1);
        const int s2i = (int)H2*((int)(S2>1)-(int)(H2>1));
        const int c2i = (int)(H2*S2)*((int)(C2>1)-(int)(S2>1));
        const int r2i = (int)(H2*S2*C2)*((int)(R2>1)-(int)(C2>1));
        for (size_t r=0; r<R; r++, X2+=r2i)
        {
            for (size_t c=0; c<C; c++, X2+=c2i)
            {
                for (size_t s=0; s<S; s++, X2+=s2i)
                {
                    for (size_t h=0; h<H; h++, X2+=h2i)
                    {
                        *X1++ *= *X2;
                    }
                }
            }
        }
    }

    return 0;
}


int times_inplace_d (double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    if (N1!=N) { fprintf(stderr,"error in times_inplace_d: first input (X1) cannot be broadcast for inplace version\n"); return 1; }

    if (N2==1)
    {
        for (size_t n=0; n<N; n++) { X1[n] *= X2[0]; }
    }
    else if (N==N2)
    {
        for (size_t n=0; n<N; n++) { X1[n] *= X2[n]; }
    }
    else if (iscolmajor)
    {
        const int r2i = (int)(R2>1);
        const int c2i = (int)R2*((int)(C2>1)-(int)(R2>1));
        const int s2i = (int)(R2*C2)*((int)(S2>1)-(int)(C2>1));
        const int h2i = (int)(R2*C2*S2)*((int)(H2>1)-(int)(S2>1));
        for (size_t h=0; h<H; h++, X2+=h2i)
        {
            for (size_t s=0; s<S; s++, X2+=s2i)
            {
                for (size_t c=0; c<C; c++, X2+=c2i)
                {
                    for (size_t r=0; r<R; r++, X2+=r2i)
                    {
                        *X1++ *= *X2;
                    }
                }
            }
        }
    }
    else
    {
        const int h2i = (int)(H2>1);
        const int s2i = (int)H2*((int)(S2>1)-(int)(H2>1));
        const int c2i = (int)(H2*S2)*((int)(C2>1)-(int)(S2>1));
        const int r2i = (int)(H2*S2*C2)*((int)(R2>1)-(int)(C2>1));
        for (size_t r=0; r<R; r++, X2+=r2i)
        {
            for (size_t c=0; c<C; c++, X2+=c2i)
            {
                for (size_t s=0; s<S; s++, X2+=s2i)
                {
                    for (size_t h=0; h<H; h++, X2+=h2i)
                    {
                        *X1++ *= *X2;
                    }
                }
            }
        }
    }

    return 0;
}


int times_inplace_c (float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    if (N1!=N) { fprintf(stderr,"error in times_inplace_c: first input (X1) cannot be broadcast for inplace version\n"); return 1; }
    
    if (N2==1)
    {
        for (size_t n2=0; n2<2*N; n2+=2) { X1[n2] *= X2[0]; X1[n2+1] *= X2[1]; }
    }
    else if (N==N2)
    {
        for (size_t n=0; n<2*N; n++) { X1[n] *= X2[n]; }
    }
    else if (iscolmajor)
    {
        const int r2i = 2*(int)(R2>1);
        const int c2i = 2*(int)R2*((int)(C2>1)-(int)(R2>1));
        const int s2i = 2*(int)(R2*C2)*((int)(S2>1)-(int)(C2>1));
        const int h2i = 2*(int)(R2*C2*S2)*((int)(H2>1)-(int)(S2>1));
        for (size_t h=0; h<H; h++, X2+=h2i)
        {
            for (size_t s=0; s<S; s++, X2+=s2i)
            {
                for (size_t c=0; c<C; c++, X2+=c2i)
                {
                    for (size_t r=0; r<R; r++, X2+=r2i)
                    {
                        *X1++ *= *X2;
                        *X1++ *= *(X2+1);
                    }
                }
            }
        }
    }
    else
    {
        const int h2i = 2*(int)(H2>1);
        const int s2i = 2*(int)H2*((int)(S2>1)-(int)(H2>1));
        const int c2i = 2*(int)(H2*S2)*((int)(C2>1)-(int)(S2>1));
        const int r2i = 2*(int)(H2*S2*C2)*((int)(R2>1)-(int)(C2>1));
        for (size_t r=0; r<R; r++, X2+=r2i)
        {
            for (size_t c=0; c<C; c++, X2+=c2i)
            {
                for (size_t s=0; s<S; s++, X2+=s2i)
                {
                    for (size_t h=0; h<H; h++, X2+=h2i)
                    {
                        *X1++ *= *X2;
                        *X1++ *= *(X2+1);
                    }
                }
            }
        }
    }

    return 0;
}


int times_inplace_z (double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    if (N1!=N) { fprintf(stderr,"error in times_inplace_z: first input (X1) cannot be broadcast for inplace version\n"); return 1; }
    
    if (N2==1)
    {
        for (size_t n2=0; n2<2*N; n2+=2) { X1[n2] *= X2[0]; X1[n2+1] *= X2[1]; }
    }
    else if (N==N2)
    {
        for (size_t n=0; n<2*N; n++) { X1[n] *= X2[n]; }
    }
    else if (iscolmajor)
    {
        const int r2i = 2*(int)(R2>1);
        const int c2i = 2*(int)R2*((int)(C2>1)-(int)(R2>1));
        const int s2i = 2*(int)(R2*C2)*((int)(S2>1)-(int)(C2>1));
        const int h2i = 2*(int)(R2*C2*S2)*((int)(H2>1)-(int)(S2>1));
        for (size_t h=0; h<H; h++, X2+=h2i)
        {
            for (size_t s=0; s<S; s++, X2+=s2i)
            {
                for (size_t c=0; c<C; c++, X2+=c2i)
                {
                    for (size_t r=0; r<R; r++, X2+=r2i)
                    {
                        *X1++ *= *X2;
                        *X1++ *= *(X2+1);
                    }
                }
            }
        }
    }
    else
    {
        const int h2i = 2*(int)(H2>1);
        const int s2i = 2*(int)H2*((int)(S2>1)-(int)(H2>1));
        const int c2i = 2*(int)(H2*S2)*((int)(C2>1)-(int)(S2>1));
        const int r2i = 2*(int)(H2*S2*C2)*((int)(R2>1)-(int)(C2>1));
        for (size_t r=0; r<R; r++, X2+=r2i)
        {
            for (size_t c=0; c<C; c++, X2+=c2i)
            {
                for (size_t s=0; s<S; s++, X2+=s2i)
                {
                    for (size_t h=0; h<H; h++, X2+=h2i)
                    {
                        *X1++ *= *X2;
                        *X1++ *= *(X2+1);
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
