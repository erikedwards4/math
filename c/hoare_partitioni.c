//Sorti_Help function.
//Paritioning part of quicksorti algorithm.
//See quicksorti.c, etc. for usage.
//This operates in-place.

#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


size_t hoare_partitioni_s (FLT_I *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    FLT_I x = X[0];
    float pivot = x.val;
    if (ascend)
    {
        while (1)
        {
            while (X[i].val<pivot) { ++i; }
            while (X[j].val>pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
    else
    {
        while (1)
        {
            while (X[i].val>pivot) { ++i; }
            while (X[j].val<pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
}


size_t hoare_partitioni_d (DBL_I *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    DBL_I x = X[0];
    double pivot = x.val;
    if (ascend)
    {
        while (1)
        {
            while (X[i].val<pivot) { ++i; }
            while (X[j].val>pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
    else
    {
        while (1)
        {
            while (X[i].val>pivot) { ++i; }
            while (X[j].val<pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
}


size_t hoare_partitioni_c (CFLT_I *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    CFLT_I x = X[0];
    float pivot = x.val;
    if (ascend)
    {
        while (1)
        {
            while (X[i].val<pivot) { ++i; }
            while (X[j].val>pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
    else
    {
        while (1)
        {
            while (X[i].val>pivot) { ++i; }
            while (X[j].val<pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
}


size_t hoare_partitioni_z (CDBL_I *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    CDBL_I x = X[0];
    double pivot = x.val;
    if (ascend)
    {
        while (1)
        {
            while (X[i].val<pivot) { ++i; }
            while (X[j].val>pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
    else
    {
        while (1)
        {
            while (X[i].val>pivot) { ++i; }
            while (X[j].val<pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
}


#ifdef __cplusplus
}
}
#endif
