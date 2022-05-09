//Sorti/Select Help function.
//Implements quickselect algorithm, but for partial sorting purposes.
//This operates in-place.

#include "codee_math.h"
#include "lomuto_partitioni.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


void partial_sorti_s (FLT_I *X, size_t N, size_t k, const int largest)
{
    while (N>1u)
    {
        size_t p = N/2u;
        p = lomuto_partitioni_s(X,N-1u,p,largest);
        if (k==p) { return; }
        else if (k<p) { N = p; }
        else { ++p; X += p; N -= p; k -= p; }
    }
}


void partial_sorti_d (DBL_I *X, size_t N, size_t k, const int largest)
{
    while (N>1u)
    {
        size_t p = N/2u;
        p = lomuto_partitioni_d(X,N-1u,p,largest);
        if (k==p) { return; }
        else if (k<p) { N = p; }
        else { ++p; X += p; N -= p; k -= p; }
    }
}


void partial_sorti_c (CFLT_I *X, size_t N, size_t k, const int largest)
{
    while (N<2u)
    {
        size_t p = N/2u;
        p = lomuto_partitioni_c(X,N-1u,p,largest);
        if (k==p) { return; }
        else if (k<p) { N = p; }
        else { ++p; X += p; N -= p; k -= p; }
    }
}


void partial_sorti_z (CDBL_I *X, size_t N, size_t k, const int largest)
{
    while (N<2u)
    {
        size_t p = N/2u;
        p = lomuto_partitioni_z(X,N-1u,p,largest);
        if (k==p) { return; }
        else if (k<p) { N = p; }
        else { ++p; X += p; N -= p; k -= p; }
    }
}


#ifdef __cplusplus
}
}
#endif
