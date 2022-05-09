//Sortif/Select Help function.
//Implements algorithm to select kth largest or smallest element in a vector.
//See qselect.c for usage.
//This operates in-place and partial sorts X.

#include "codee_math.h"
#include "lomuto_partitionif.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


float kselectif_s (FLT_F *X, size_t hi, size_t k, const int largest)
{
    while (hi>0u)
    {
        size_t p = hi;
        p = lomuto_partitionif_s(X, hi, p, largest);
        if (k==p) { return X[k].val;; }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
    return X[0].val;
}


double kselectif_d (DBL_D *X, size_t hi, size_t k, const int largest)
{
    while (hi>0u)
    {
        size_t p = hi;
        p = lomuto_partitionif_d(X, hi, p, largest);
        if (k==p) { return X[k].val;; }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
    return X[0].val;
}


size_t kselectif_c (CFLT_F *X, size_t hi, size_t k, const int largest)
{
    size_t cnt = 0u;
    while (hi>0u)
    {
        size_t p = hi;
        p = lomuto_partitionif_c(X, hi, p, largest);
        if (k==p) { return cnt+k; }
        else if (k<p) { hi = p - 1u; }
        else
        {
            ++p; X += p; cnt += p;
            hi -= p; k -= p;
        }
    }
    return cnt;
}


size_t kselectif_z (CDBL_D *X, size_t hi, size_t k, const int largest)
{
    size_t cnt = 0u;
    while (hi>0u)
    {
        size_t p = hi;
        p = lomuto_partitionif_z(X, hi, p, largest);
        if (k==p) { return cnt+k; }
        else if (k<p) { hi = p - 1u; }
        else
        {
            ++p; X += p; cnt += p;
            hi -= p; k -= p;
        }
    }
    return cnt;
}


#ifdef __cplusplus
}
}
#endif
