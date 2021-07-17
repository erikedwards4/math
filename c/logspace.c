//Generates vector Y with elements logarithmically spaced from 10^a to 10^b.
//For complex Y, imag part is set to 0.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int logspace_s (float *Y, const size_t N, const float a, const float b);
int logspace_d (double *Y, const size_t N, const double a, const double b);
int logspace_c (float *Y, const size_t N, const float a, const float b);
int logspace_z (double *Y, const size_t N, const double a, const double b);


int logspace_s (float *Y, const size_t N, const float a, const float b)
{
    const float stp = (b-a)/(N-1u);
    float y = a;

    for (size_t n=0u; n<N-1u; ++n, ++Y, y+=stp) { *Y = powf(10.0f,y); }
    *Y = powf(10.0f,b);

    return 0;
}


int logspace_d (double *Y, const size_t N, const double a, const double b)
{
    const double stp = (b-a)/(N-1u);
    double y = a;

    for (size_t n=0u; n<N-1u; ++n, ++Y, y+=stp) { *Y = pow(10.0,y); }
    *Y = pow(10.0,b);

    return 0;
}


int logspace_c (float *Y, const size_t N, const float a, const float b)
{
    const float stp = (b-a)/(N-1u);
    float yr = a;

    for (size_t n=0u; n<N-1u; ++n, ++Y, yr+=stp) { *Y = powf(10.0f,yr); *++Y = 0.0f; }
    *Y = powf(10.0f,b); *++Y = 0.0f;

    return 0;
}


int logspace_z (double *Y, const size_t N, const double a, const double b)
{
    const double stp = (b-a)/(N-1u);
    double yr = a;

    for (size_t n=0u; n<N-1u; ++n, ++Y, yr+=stp) { *Y = pow(10.0,yr); *++Y = 0.0; }
    *Y = pow(10.0,b); *++Y = 0.0;

    return 0;
}


#ifdef __cplusplus
}
}
#endif
