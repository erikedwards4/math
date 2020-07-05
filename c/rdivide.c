//2-input elementwise function.
//Does elementwise division: Y = X1./X2
//This has in-place and not-in-place versions.

#include <stdio.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int rdivide_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int rdivide_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int rdivide_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int rdivide_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);

int rdivide_inplace_s (float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int rdivide_inplace_d (double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int rdivide_inplace_c (float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);
int rdivide_inplace_z (double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor);


int rdivide_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);
    if (N1==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X1[0] / X2[n]; }
    }    
    else if (N2==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X1[n] / X2[0]; }
    }
    else if (N==N1 && N==N2)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X1[n] / X2[n]; }
    }
    else if (iscolmajor)
    {
        const size_t RCSH1 = R1*C1*S1*(H1>1), RCSH2 = R2*C2*S2*(H2>1);
        const size_t RCS1 = R1*C1*(S1>1), RCS2 = R2*C2*(S2>1);
        const size_t RC1 = R1*(C1>1), RC2 = R2*(C2>1);
        const size_t r1i = (R1>1), r2i = (R2>1);
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                for (size_t c=0; c<C; c++)
                {
                    size_t n1 = h*RCSH1 + s*RCS1 + c*RC1;
                    size_t n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (size_t r=0; r<R; r++, n++, n1+=r1i, n2+=r2i)
                    {
                        Y[n] = X1[n1] / X2[n2];
                    }
                }
            }
        }
    }
    else
    {
        const size_t HSCR1 = H1*S1*C1*(R1>1), HSCR2 = H2*S2*C2*(R2>1);
        const size_t HSC1 = H1*S1*(C1>1), HSC2 = H2*S2*(C2>1);
        const size_t HS1 = H1*(S1>1), HS2 = H2*(S2>1);
        const size_t h1i = (H1>1), h2i = (H2>1);
        for (size_t r=0, n=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++)
            {
                for (size_t s=0; s<S; s++)
                {
                    size_t n1 = r*HSCR1 + c*HSC1 + s*HS1;
                    size_t n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (size_t h=0; h<H; h++, n++, n1+=h1i, n2+=h2i)
                    {
                        Y[n] = X1[n1] / X2[n2];
                    }
                }
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int rdivide_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    
    if (N1==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X1[0] / X2[n]; }
    }    
    else if (N2==1)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X1[n] / X2[0]; }
    }
    else if (N==N1 && N==N2)
    {
        for (size_t n=0; n<N; n++) { Y[n] = X1[n] / X2[n]; }
    }
    else if (iscolmajor)
    {
        const size_t RCSH1 = R1*C1*S1*(H1>1), RCSH2 = R2*C2*S2*(H2>1);
        const size_t RCS1 = R1*C1*(S1>1), RCS2 = R2*C2*(S2>1);
        const size_t RC1 = R1*(C1>1), RC2 = R2*(C2>1);
        const size_t r1i = (R1>1), r2i = (R2>1);
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                for (size_t c=0; c<C; c++)
                {
                    size_t n1 = h*RCSH1 + s*RCS1 + c*RC1;
                    size_t n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (size_t r=0; r<R; r++, n++, n1+=r1i, n2+=r2i)
                    {
                        Y[n] = X1[n1] / X2[n2];
                    }
                }
            }
        }
    }
    else
    {
        const size_t HSCR1 = H1*S1*C1*(R1>1), HSCR2 = H2*S2*C2*(R2>1);
        const size_t HSC1 = H1*S1*(C1>1), HSC2 = H2*S2*(C2>1);
        const size_t HS1 = H1*(S1>1), HS2 = H2*(S2>1);
        const size_t h1i = (H1>1), h2i = (H2>1);
        for (size_t r=0, n=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++)
            {
                for (size_t s=0; s<S; s++)
                {
                    size_t n1 = r*HSCR1 + c*HSC1 + s*HS1;
                    size_t n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (size_t h=0; h<H; h++, n++, n1+=h1i, n2+=h2i)
                    {
                        Y[n] = X1[n1] / X2[n2];
                    }
                }
            }
        }
    }

    return 0;
}


int rdivide_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
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
    else if (N==N1 && N==N2)
    {
        for (size_t n=0; n<2*N; n++) { Y[n] = X1[n] / X2[n]; }
    }
    else if (iscolmajor)
    {
        const size_t RCSH1 = 2*R1*C1*S1*(H1>1), RCSH2 = 2*R2*C2*S2*(H2>1);
        const size_t RCS1 = 2*R1*C1*(S1>1), RCS2 = 2*R2*C2*(S2>1);
        const size_t RC1 = 2*R1*(C1>1), RC2 = 2*R2*(C2>1);
        const size_t r1i = 2*(R1>1), r2i = 2*(R2>1);
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                for (size_t c=0; c<C; c++)
                {
                    size_t n1 = h*RCSH1 + s*RCS1 + c*RC1;
                    size_t n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (size_t r=0; r<R; r++, n+=2, n1+=r1i, n2+=r2i)
                    {
                        Y[n] = X1[n1] / X2[n2]; Y[n+1] = X1[n1+1] / X2[n2+1];
                    }
                }
            }
        }
    }
    else
    {
        const size_t HSCR1 = 2*H1*S1*C1*(R1>1), HSCR2 = 2*H2*S2*C2*(R2>1);
        const size_t HSC1 = 2*H1*S1*(C1>1), HSC2 = 2*H2*S2*(C2>1);
        const size_t HS1 = 2*H1*(S1>1), HS2 = 2*H2*(S2>1);
        const size_t h1i = 2*(H1>1), h2i = 2*(H2>1);
        for (size_t r=0, n=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++)
            {
                for (size_t s=0; s<S; s++)
                {
                    size_t n1 = r*HSCR1 + c*HSC1 + s*HS1;
                    size_t n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (size_t h=0; h<H; h++, n+=2, n1+=h1i, n2+=h2i)
                    {
                        Y[n] = X1[n1] / X2[n2]; Y[n+1] = X1[n1+1] / X2[n2+1];
                    }
                }
            }
        }
    }

    return 0;
}


int rdivide_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
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
    else if (N==N1 && N==N2)
    {
        for (size_t n=0; n<2*N; n++) { Y[n] = X1[n] / X2[n]; }
    }
    else if (iscolmajor)
    {
        const size_t RCSH1 = 2*R1*C1*S1*(H1>1), RCSH2 = 2*R2*C2*S2*(H2>1);
        const size_t RCS1 = 2*R1*C1*(S1>1), RCS2 = 2*R2*C2*(S2>1);
        const size_t RC1 = 2*R1*(C1>1), RC2 = 2*R2*(C2>1);
        const size_t r1i = 2*(R1>1), r2i = 2*(R2>1);
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                for (size_t c=0; c<C; c++)
                {
                    size_t n1 = h*RCSH1 + s*RCS1 + c*RC1;
                    size_t n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (size_t r=0; r<R; r++, n+=2, n1+=r1i, n2+=r2i)
                    {
                        Y[n] = X1[n1] / X2[n2]; Y[n+1] = X1[n1+1] / X2[n2+1];
                    }
                }
            }
        }
    }
    else
    {
        const size_t HSCR1 = 2*H1*S1*C1*(R1>1), HSCR2 = 2*H2*S2*C2*(R2>1);
        const size_t HSC1 = 2*H1*S1*(C1>1), HSC2 = 2*H2*S2*(C2>1);
        const size_t HS1 = 2*H1*(S1>1), HS2 = 2*H2*(S2>1);
        const size_t h1i = 2*(H1>1), h2i = 2*(H2>1);
        for (size_t r=0, n=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++)
            {
                for (size_t s=0; s<S; s++)
                {
                    size_t n1 = r*HSCR1 + c*HSC1 + s*HS1;
                    size_t n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (size_t h=0; h<H; h++, n+=2, n1+=h1i, n2+=h2i)
                    {
                        Y[n] = X1[n1] / X2[n2]; Y[n+1] = X1[n1+1] / X2[n2+1];
                    }
                }
            }
        }
    }

    return 0;
}


int rdivide_inplace_s (float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    if (N1!=N) { fprintf(stderr,"error in rdivide_inplace_s: first input (X1) cannot be broadcast for inplace version\n"); return 1; }
    
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);
    if (N2==1)
    {
        for (size_t n=0; n<N; n++) { X1[n] /= X2[0]; }
    }
    else if (N==N2)
    {
        for (size_t n=0; n<N; n++) { X1[n] /= X2[n]; }
    }
    else if (iscolmajor)
    {
        const size_t RCSH2 = R2*C2*S2*(H2>1), RCS2 = R2*C2*(S2>1), RC2 = R2*(C2>1), r2i = (R2>1);
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                for (size_t c=0; c<C; c++)
                {
                    size_t n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (size_t r=0; r<R; r++, n++, n2+=r2i)
                    {
                        X1[n] /= X2[n2];
                    }
                }
            }
        }
    }
    else
    {
        const size_t HSCR2 = H2*S2*C2*(R2>1), HSC2 = H2*S2*(C2>1), HS2 = H2*(S2>1), h2i = (H2>1);
        for (size_t r=0, n=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++)
            {
                for (size_t s=0; s<S; s++)
                {
                    size_t n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (size_t h=0; h<H; h++, n++, n2+=h2i)
                    {
                        X1[n] /= X2[n2];
                    }
                }
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int rdivide_inplace_d (double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    if (N1!=N) { fprintf(stderr,"error in rdivide_inplace_d: first input (X1) cannot be broadcast for inplace version\n"); return 1; }
    
    if (N2==1)
    {
        for (size_t n=0; n<N; n++) { X1[n] /= X2[0]; }
    }
    else if (N==N2)
    {
        for (size_t n=0; n<N; n++) { X1[n] /= X2[n]; }
    }
    else if (iscolmajor)
    {
        const size_t RCSH2 = R2*C2*S2*(H2>1), RCS2 = R2*C2*(S2>1), RC2 = R2*(C2>1), r2i = (R2>1);
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                for (size_t c=0; c<C; c++)
                {
                    size_t n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (size_t r=0; r<R; r++, n++, n2+=r2i)
                    {
                        X1[n] /= X2[n2];
                    }
                }
            }
        }
    }
    else
    {
        const size_t HSCR2 = H2*S2*C2*(R2>1), HSC2 = H2*S2*(C2>1), HS2 = H2*(S2>1), h2i = (H2>1);
        for (size_t r=0, n=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++)
            {
                for (size_t s=0; s<S; s++)
                {
                    size_t n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (size_t h=0; h<H; h++, n++, n2+=h2i)
                    {
                        X1[n] /= X2[n2];
                    }
                }
            }
        }
    }

    return 0;
}


int rdivide_inplace_c (float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    if (N1!=N) { fprintf(stderr,"error in rdivide_inplace_c: first input (X1) cannot be broadcast for inplace version\n"); return 1; }
    
    if (N2==1)
    {
        for (size_t n2=0; n2<2*N; n2+=2) { X1[n2] /= X2[0]; X1[n2+1] /= X2[1]; }
    }
    else if (N==N1 && N==N2)
    {
        for (size_t n=0; n<2*N; n++) { X1[n] /= X2[n]; }
    }
    else if (iscolmajor)
    {
        const size_t RCSH2 = 2*R2*C2*S2*(H2>1), RCS2 = 2*R2*C2*(S2>1), RC2 = 2*R2*(C2>1), r2i = 2*(R2>1);
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                for (size_t c=0; c<C; c++)
                {
                    size_t n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (size_t r=0; r<R; r++, n+=2, n2+=r2i)
                    {
                        X1[n] /= X2[n2]; X1[n+1] /= X2[n2+1];
                    }
                }
            }
        }
    }
    else
    {
        const size_t HSCR2 = 2*H2*S2*C2*(R2>1), HSC2 = 2*H2*S2*(C2>1), HS2 = 2*H2*(S2>1), h2i = 2*(H2>1);
        for (size_t r=0, n=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++)
            {
                for (size_t s=0; s<S; s++)
                {
                    size_t n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (size_t h=0; h<H; h+=2, n++, n2+=h2i)
                    {
                        X1[n] /= X2[n2]; X1[n+1] /= X2[n2+1];
                    }
                }
            }
        }
    }

    return 0;
}


int rdivide_inplace_z (double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor)
{
    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    if (N1!=N) { fprintf(stderr,"error in rdivide_inplace_z: first input (X1) cannot be broadcast for inplace version\n"); return 1; }
    
    if (N2==1)
    {
        for (size_t n2=0; n2<2*N; n2+=2) { X1[n2] /= X2[0]; X1[n2+1] /= X2[1]; }
    }
    else if (N==N1 && N==N2)
    {
        for (size_t n=0; n<2*N; n++) { X1[n] /= X2[n]; }
    }
    else if (iscolmajor)
    {
        const size_t RCSH2 = 2*R2*C2*S2*(H2>1), RCS2 = 2*R2*C2*(S2>1), RC2 = 2*R2*(C2>1), r2i = 2*(R2>1);
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                for (size_t c=0; c<C; c++)
                {
                    size_t n2 = h*RCSH2 + s*RCS2 + c*RC2;
                    for (size_t r=0; r<R; r++, n+=2, n2+=r2i)
                    {
                        X1[n] /= X2[n2]; X1[n+1] /= X2[n2+1];
                    }
                }
            }
        }
    }
    else
    {
        const size_t HSCR2 = 2*H2*S2*C2*(R2>1), HSC2 = 2*H2*S2*(C2>1), HS2 = 2*H2*(S2>1), h2i = 2*(H2>1);
        for (size_t r=0, n=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++)
            {
                for (size_t s=0; s<S; s++)
                {
                    size_t n2 = r*HSCR2 + c*HSC2 + s*HS2;
                    for (size_t h=0; h<H; h++, n+=2, n2+=h2i)
                    {
                        X1[n] /= X2[n2]; X1[n+1] /= X2[n2+1];
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
