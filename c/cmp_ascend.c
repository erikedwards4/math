//Sort_Help functions.
//Give comparison (cmp) functions for use by C stdlib qsort.

//#include <string.h>
#include <math.h>
//#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int cmp_ascend_s (const void *a, const void *b)
{
	const float x1 = *(const float*)a, x2 = *(const float*)b;
	if (x1!=x1) { return 1; }
    else if (x2!=x2) { return -1; }
    else if (x1>x2) { return 1; }
    else if (x2>x1) { return -1; }
    else { return 0; }
}


int cmp_ascend_d (const void *a, const void *b)
{
	const double x1 = *(const double*)a, x2 = *(const double*)b;
	if (x1!=x1) { return 1; }
    else if (x2!=x2) { return -1; }
    else if (x1>x2) { return 1; }
    else if (x2>x1) { return -1; }
    else { return 0; }
}


// static int cmp_ascend_c (const void *a, const void *b)
// {
//     float x1[2], x2[2], abs1, abs2;
//     memcpy(x1,a,2*sizeof(float)); memcpy(x2,b,2*sizeof(float));
//     abs1 = sqrtf(x1[0]*x1[0] + x1[1]*x1[1]);
//     abs2 = sqrtf(x2[0]*x2[0] + x2[1]*x2[1]);
// 	if (abs1!=abs1) { return 1; }
//     else if (abs2!=abs2) { return -1; }
//     else if (abs1>abs2) { return 1; }
//     else if (abs2>abs1) { return -1; }
// 	else if (atan2f(x1[1],x1[0])>atan2f(x2[1],x2[0])) { return 1; }
//     else if (atan2f(x2[1],x2[0])>atan2f(x1[1],x1[0])) { return -1; }
//     else { return 0; }
// }


int cmp_ascend_c (const void *a, const void *b)
{
	const float abs1 = sqrtf(*(const float *)a**(const float *)a + *((const float *)(a)+1)**((const float *)(a)+1));
    const float abs2 = sqrtf(*(const float *)b**(const float *)b + *((const float *)(b)+1)**((const float *)(b)+1));
	if (abs1!=abs1) { return 1; }
    else if (abs2!=abs2) { return -1; }
    else if (abs1>abs2) { return 1; }
    else if (abs2>abs1) { return -1; }
	else if (atan2f(*((const float *)(a)+1),*(const float *)a) > atan2f(*((const float *)(b)+1),*(const float *)b)) { return 1; }
    else if (atan2f(*((const float *)(b)+1),*(const float *)b) > atan2f(*((const float *)(a)+1),*(const float *)a)) { return -1; }
    else { return 0; }
}


int cmp_ascend_z (const void *a, const void *b)
{
	const double abs1 = sqrt(*(const double *)a**(const double *)a + *((const double *)(a)+1)**((const double *)(a)+1));
    const double abs2 = sqrt(*(const double *)b**(const double *)b + *((const double *)(b)+1)**((const double *)(b)+1));
	if (abs1!=abs1) { return 1; }
    else if (abs2!=abs2) { return -1; }
    else if (abs1>abs2) { return 1; }
    else if (abs2>abs1) { return -1; }
	else if (atan2(*((const double *)(a)+1),*(const double *)a) > atan2(*((const double *)(b)+1),*(const double *)b)) { return 1; }
    else if (atan2(*((const double *)(b)+1),*(const double *)b) > atan2(*((const double *)(a)+1),*(const double *)a)) { return -1; }
    else { return 0; }
}


#ifdef __cplusplus
}
}
#endif
