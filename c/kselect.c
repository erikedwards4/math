//Sort/Select Help function.
//Implements algorithm to select kth largest or smallest element in a vector.
//See qselect.c for usage.

#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

static size_t lomuto_partition_s (float *X, const size_t hi, size_t p, int largest);
static size_t lomuto_partition_d (double *X, const size_t hi, size_t p, int largest);
static size_t lomuto_partition_c (float *X, const size_t hi, size_t p, int largest);
static size_t lomuto_partition_z (double *X, const size_t hi, size_t p, int largest);


static size_t lomuto_partition_s (float *X, const size_t hi, size_t p, int largest)
{
    float x = X[p], pivot = x;
    X[p] = X[hi]; X[hi] = x;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i]>=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i]<=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    x = X[p]; X[p] = X[hi]; X[hi] = x;
    return p;
}


static size_t lomuto_partition_d (double *X, const size_t hi, size_t p, int largest)
{
    double x = X[p], pivot = x;
    X[p] = X[hi]; X[hi] = x;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i]>=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i]<=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    x = X[p]; X[p] = X[hi]; X[hi] = x;
    return p;
}


static size_t lomuto_partition_c (float *X, const size_t hi, size_t p, int largest)
{
    float xr = X[2u*p], xi = X[2u*p+1u], pivot = xr*xr + xi*xi;
    X[2u*p] = X[2u*hi]; X[2u*hi] = xr;
    X[2u*p+1u] = X[2u*hi+1u]; X[2u*hi+1u] = xi;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[2u*i]*X[2u*i]+X[2u*i+1u]*X[2u*i+1u]>=pivot+1e-5f)
            {
                xr = X[2u*i]; X[2u*i] = X[2u*p]; X[2u*p] = xr;
                xi = X[2u*i+1u]; X[2u*i+1u] = X[2u*p+1u]; X[2u*p+1u] = xi;
                ++p;
            }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[2u*i]*X[2u*i]+X[2u*i+1u]*X[2u*i+1u]<=pivot-1e-5f)
            {
                xr = X[2u*i]; X[2u*i] = X[2u*p]; X[2u*p] = xr;
                xi = X[2u*i+1u]; X[2u*i+1u] = X[2u*p+1u]; X[2u*p+1u] = xi;
                ++p;
            }
        }
    }
    xr = X[2u*p]; X[2u*p] = X[2u*hi]; X[2u*hi] = xr;
    xi = X[2u*p+1u]; X[2u*p+1u] = X[2u*hi+1u]; X[2u*hi+1u] = xi;
    return p;
}


static size_t lomuto_partition_z (double *X, const size_t hi, size_t p, int largest)
{
    double xr = X[2u*p], xi = X[2u*p+1u], pivot = xr*xr + xi*xi;
    X[2u*p] = X[2u*hi]; X[2u*hi] = xr;
    X[2u*p+1u] = X[2u*hi+1u]; X[2u*hi+1u] = xi;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[2u*i]*X[2u*i]+X[2u*i+1u]*X[2u*i+1u]>=pivot)
            {
                xr = X[2u*i]; X[2u*i] = X[2u*p]; X[2u*p] = xr;
                xi = X[2u*i+1u]; X[2u*i+1u] = X[2u*p+1u]; X[2u*p+1u] = xi;
                ++p;
            }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[2u*i]*X[2u*i]+X[2u*i+1u]*X[2u*i+1u]<=pivot)
            {
                xr = X[2u*i]; X[2u*i] = X[2u*p]; X[2u*p] = xr;
                xi = X[2u*i+1u]; X[2u*i+1u] = X[2u*p+1u]; X[2u*p+1u] = xi;
                ++p;
            }
        }
    }
    xr = X[2u*p]; X[2u*p] = X[2u*hi]; X[2u*hi] = xr;
    xi = X[2u*p+1u]; X[2u*p+1u] = X[2u*hi+1u]; X[2u*hi+1u] = xi;
    return p;
}


float kselect_s (float *X, size_t hi, size_t k, int largest)
{
    while (1)
    {
        if (hi==0u) { return *X; }
        size_t p = hi/2u;
        p = lomuto_partition_s(X,hi,p,largest);
        if (k==p) { return *(X+k); }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
}


double kselect_d (double *X, size_t hi, size_t k, int largest)
{
    while (1)
    {
        if (hi==0u) { return *X; }
        size_t p = hi/2u;
        p = lomuto_partition_d(X,hi,p,largest);
        if (k==p) { return *(X+k); }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
}


size_t kselect_c (float *X, size_t hi, size_t k, int largest)
{
    size_t cnt = 0u;
    while (1)
    {
        if (hi==0u) { return cnt; }
        size_t p = hi/2u;
        p = lomuto_partition_c(X,hi,p,largest);
        if (k==p) { return cnt+2u*k; }
        else if (k<p) { hi = p - 1u; }
        else
        {
            ++p; X += 2u*p; cnt += 2u*p;
            hi -= p; k -= p;
        }
    }
}


size_t kselect_z (double *X, size_t hi, size_t k, int largest)
{
    size_t cnt = 0u;
    while (1)
    {
        if (hi==0u) { return cnt; }
        size_t p = hi/2u;
        p = lomuto_partition_z(X,hi,p,largest);
        if (k==p) { return cnt+2u*k; }
        else if (k<p) { hi = p - 1u; }
        else
        {
            ++p; X += 2u*p; cnt += 2u*p;
            hi -= p; k -= p;
        }
    }
}


#ifdef __cplusplus
}
}
#endif
