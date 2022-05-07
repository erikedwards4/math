//Sort_Help functions.
//Give comparison (cmp) functions using structs that include the index.
//See sorti.c and ranks.c for usage.

#include <math.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

//It is faster to keep the inds as the same data type as the vals for sorti,
//but ranks requires the size_t ind.
// typedef struct { float val; float ind; } FLT;
// typedef struct { double val; double ind; } DBL;
// typedef struct { float r; float i; float ind; } CFLT;
// typedef struct { double r; double i; double ind; } ZDBL;

typedef struct { float val; size_t ind; } FLT;
typedef struct { double val; size_t ind; } DBL;
typedef struct { float r; float i; size_t ind; } CFLT;
typedef struct { double r; double i; size_t ind; } ZDBL;


int cmpi_ascend_s (const void *a, const void *b)
{
    const FLT x1 = *(const FLT *)a;
    const FLT x2 = *(const FLT *)b;
    if (x1.val>x2.val) { return 1; }
    else if (x2.val>x1.val) { return -1; }
    else { return 0; }
}


int cmpi_ascend_d (const void *a, const void *b)
{
	const DBL x1 = *(const DBL *)a;
    const DBL x2 = *(const DBL *)b;
    if (x1.val>x2.val) { return 1; }
    else if (x2.val>x1.val) { return -1; }
    else { return 0; }
}


int cmpi_ascend_c (const void *a, const void *b)
{
    const CFLT x1 = *(const CFLT *)a;
    const CFLT x2 = *(const CFLT *)b;
	const float sq1 = x1.r*x1.r + x1.i*x1.i;
    const float sq2 = x2.r*x2.r + x2.i*x2.i;
	if (sq1>sq2) { return 1; }
    else if (sq2>sq1) { return -1; }
	else if (atan2f(x1.i,x1.r) > atan2f(x2.i,x2.r)) { return 1; }
    else if (atan2f(x2.i,x2.r) > atan2f(x1.i,x1.r)) { return -1; }
    else { return 0; }
}


int cmpi_ascend_z (const void *a, const void *b)
{
	const ZDBL x1 = *(const ZDBL *)a;
    const ZDBL x2 = *(const ZDBL *)b;
	const double sq1 = x1.r*x1.r + x1.i*x1.i;
    const double sq2 = x2.r*x2.r + x2.i*x2.i;
	if (sq1>sq2) { return 1; }
    else if (sq2>sq1) { return -1; }
	else if (atan2(x1.i,x1.r) > atan2(x2.i,x2.r)) { return 1; }
    else if (atan2(x2.i,x2.r) > atan2(x1.i,x1.r)) { return -1; }
    else { return 0; }
}


#ifdef __cplusplus
}
}
#endif
