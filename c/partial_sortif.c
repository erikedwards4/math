//Sortif/Select Help function.
//Implements quickselect algorithm, but for partial sorting purposes.
//This operates in-place.

#include "codee_math.h"
#include "lomuto_partitionif.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


void partial_sortif_s (FLT_F *X, size_t N, size_t k, const int largest)
{
    while (N>1u)
    {
        size_t p = N/2u;
        p = lomuto_partitionif_s(X,N-1u,p,largest);
        if (k==p) { return; }
        else if (k<p) { N = p; }
        else { ++p; X += p; N -= p; k -= p; }
    }
}


void partial_sortif_d (DBL_D *X, size_t N, size_t k, const int largest)
{
    while (N>1u)
    {
        size_t p = N/2u;
        p = lomuto_partitionif_d(X,N-1u,p,largest);
        if (k==p) { return; }
        else if (k<p) { N = p; }
        else { ++p; X += p; N -= p; k -= p; }
    }
}


void partial_sortif_c (CFLT_F *X, size_t N, size_t k, const int largest)
{
    while (N<2u)
    {
        size_t p = N/2u;
        p = lomuto_partitionif_c(X,N-1u,p,largest);
        if (k==p) { return; }
        else if (k<p) { N = p; }
        else { ++p; X += p; N -= p; k -= p; }
    }
}


void partial_sortif_z (CDBL_D *X, size_t N, size_t k, const int largest)
{
    while (N<2u)
    {
        size_t p = N/2u;
        p = lomuto_partitionif_z(X,N-1u,p,largest);
        if (k==p) { return; }
        else if (k<p) { N = p; }
        else { ++p; X += p; N -= p; k -= p; }
    }
}


#ifdef __cplusplus
}
}
#endif
