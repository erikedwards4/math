//Sort_Help functions.
//Give comparison (cmp) functions for use by C stdlib qsort.

#include <math.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int cmp_descend_s (const void *a, const void *b)
{
	const float x1 = *(const float*)a, x2 = *(const float*)b;
	if (x1!=x1) { return -1; }
    else if (x2!=x2) { return 1; }
    else if (x1<x2) { return 1; }
    else if (x2<x1) { return -1; }
    else { return 0; }
}


int cmp_descend_d (const void *a, const void *b)
{
	const double x1 = *(const double*)a, x2 = *(const double*)b;
	if (x1!=x1) { return -1; }
    else if (x2!=x2) { return 1; }
    else if (x1<x2) { return 1; }
    else if (x2<x1) { return -1; }
    else { return 0; }
}


int cmp_descend_c (const void *a, const void *b)
{
	const float abs1 = sqrtf(*(const float *)a**(const float *)a + *((const float *)(a)+1)**((const float *)(a)+1));
    const float abs2 = sqrtf(*(const float *)b**(const float *)b + *((const float *)(b)+1)**((const float *)(b)+1));
	if (abs1!=abs1) { return -1; }
    else if (abs2!=abs2) { return 1; }
    else if (abs1<abs2) { return 1; }
    else if (abs2<abs1) { return -1; }
	else if (atan2f(*((const float *)(a)+1),*(const float *)a) < atan2f(*((const float *)(b)+1),*(const float *)b)) { return 1; }
    else if (atan2f(*((const float *)(b)+1),*(const float *)b) < atan2f(*((const float *)(a)+1),*(const float *)a)) { return -1; }
    else { return 0; }
}


int cmp_descend_z (const void *a, const void *b)
{
	const double abs1 = sqrt(*(const double *)a**(const double *)a + *((const double *)(a)+1)**((const double *)(a)+1));
    const double abs2 = sqrt(*(const double *)b**(const double *)b + *((const double *)(b)+1)**((const double *)(b)+1));
	if (abs1!=abs1) { return -1; }
    else if (abs2!=abs2) { return 1; }
    else if (abs1<abs2) { return 1; }
    else if (abs2<abs1) { return -1; }
	else if (atan2(*((const double *)(a)+1),*(const double *)a) < atan2(*((const double *)(b)+1),*(const double *)b)) { return 1; }
    else if (atan2(*((const double *)(b)+1),*(const double *)b) < atan2(*((const double *)(a)+1),*(const double *)a)) { return -1; }
    else { return 0; }
}


#ifdef __cplusplus
}
}
#endif
