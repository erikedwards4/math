//Generates vector Y with elements linearly spaced from a to b.
//For complex Y, imag part is set to 0.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int linspace_s (float *Y, const size_t N, const float a, const float b);
int linspace_d (double *Y, const size_t N, const double a, const double b);
int linspace_c (float *Y, const size_t N, const float a, const float b);
int linspace_z (double *Y, const size_t N, const double a, const double b);


int linspace_s (float *Y, const size_t N, const float a, const float b)
{
    const float stp = (b-a)/(N-1);
    float y = a;

    for (size_t n=0u; n<N-1; ++n, ++Y, y+=stp) { *Y = y; }
    *Y = b;

    return 0;
}


int linspace_d (double *Y, const size_t N, const double a, const double b)
{
    const double stp = (b-a)/(N-1);
    double y = a;

    for (size_t n=0u; n<N-1; ++n, ++Y, y+=stp) { *Y = y; }
    *Y = b;

    return 0;
}


int linspace_c (float *Y, const size_t N, const float a, const float b)
{
    const float stp = (b-a)/(N-1);
    float yr = a;

    for (size_t n=0u; n<N-1; ++n, ++Y, yr+=stp) { *Y = yr; *++Y = 0.0f; }
    *Y = b; *++Y = 0.0f;

    return 0;
}


int linspace_z (double *Y, const size_t N, const double a, const double b)
{
    const double stp = (b-a)/(N-1);
    double yr = a;

    for (size_t n=0u; n<N-1; ++n, ++Y, yr+=stp) { *Y = yr; *++Y = 0.0; }
    *Y = b; *++Y = 0.0;

    return 0;
}


#ifdef __cplusplus
}
}
#endif
