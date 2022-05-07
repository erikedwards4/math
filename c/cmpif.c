//Sort_Help functions.
//Give comparison (cmp) functions using structs that include the index (cmpi),
//where the index is represented as a float (cmpif) for faster speed in sorti.
//See sorti.c for usage.

#include <math.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


typedef struct { float val; float ind; } FLT_F;
typedef struct { double val; double ind; } DBL_D;
typedef struct { float r; float i; float ind; } CFLT_F;
typedef struct { double r; double i; double ind; } CDBL_D;


int cmpif_ascend_s (const void *a, const void *b)
{
    const FLT_F x1 = *(const FLT_F *)a;
    const FLT_F x2 = *(const FLT_F *)b;
    if (x1.val>x2.val) { return 1; }
    else if (x2.val>x1.val) { return -1; }
    else { return 0; }
}


int cmpif_ascend_d (const void *a, const void *b)
{
	const DBL_D x1 = *(const DBL_D *)a;
    const DBL_D x2 = *(const DBL_D *)b;
    if (x1.val>x2.val) { return 1; }
    else if (x2.val>x1.val) { return -1; }
    else { return 0; }
}


int cmpif_ascend_c (const void *a, const void *b)
{
    const CFLT_F x1 = *(const CFLT_F *)a;
    const CFLT_F x2 = *(const CFLT_F *)b;
	const float sq1 = x1.r*x1.r + x1.i*x1.i;
    const float sq2 = x2.r*x2.r + x2.i*x2.i;
	if (sq1>sq2) { return 1; }
    else if (sq2>sq1) { return -1; }
	else if (atan2f(x1.i,x1.r) > atan2f(x2.i,x2.r)) { return 1; }
    else if (atan2f(x2.i,x2.r) > atan2f(x1.i,x1.r)) { return -1; }
    else { return 0; }
}


int cmpif_ascend_z (const void *a, const void *b)
{
	const CDBL_D x1 = *(const CDBL_D *)a;
    const CDBL_D x2 = *(const CDBL_D *)b;
	const double sq1 = x1.r*x1.r + x1.i*x1.i;
    const double sq2 = x2.r*x2.r + x2.i*x2.i;
	if (sq1>sq2) { return 1; }
    else if (sq2>sq1) { return -1; }
	else if (atan2(x1.i,x1.r) > atan2(x2.i,x2.r)) { return 1; }
    else if (atan2(x2.i,x2.r) > atan2(x1.i,x1.r)) { return -1; }
    else { return 0; }
}


int cmpif_descend_s (const void *a, const void *b)
{
    const FLT_F x1 = *(const FLT_F *)a;
    const FLT_F x2 = *(const FLT_F *)b;
    if (x1.val>x2.val) { return -1; }
    else if (x2.val>x1.val) { return 1; }
    else { return 0; }
}


int cmpif_descend_d (const void *a, const void *b)
{
	const DBL_D x1 = *(const DBL_D *)a;
    const DBL_D x2 = *(const DBL_D *)b;
    if (x1.val>x2.val) { return -1; }
    else if (x2.val>x1.val) { return 1; }
    else { return 0; }
}


int cmpif_descend_c (const void *a, const void *b)
{
    const CFLT_F x1 = *(const CFLT_F *)a;
    const CFLT_F x2 = *(const CFLT_F *)b;
	const float sq1 = x1.r*x1.r + x1.i*x1.i;
    const float sq2 = x2.r*x2.r + x2.i*x2.i;
	if (sq1>sq2) { return -1; }
    else if (sq2>sq1) { return 1; }
	else if (atan2f(x1.i,x1.r) > atan2f(x2.i,x2.r)) { return -1; }
    else if (atan2f(x2.i,x2.r) > atan2f(x1.i,x1.r)) { return 1; }
    else { return 0; }
}


int cmpif_descend_z (const void *a, const void *b)
{
	const CDBL_D x1 = *(const CDBL_D *)a;
    const CDBL_D x2 = *(const CDBL_D *)b;
	const double sq1 = x1.r*x1.r + x1.i*x1.i;
    const double sq2 = x2.r*x2.r + x2.i*x2.i;
	if (sq1>sq2) { return -1; }
    else if (sq2>sq1) { return 1; }
	else if (atan2(x1.i,x1.r) > atan2(x2.i,x2.r)) { return -1; }
    else if (atan2(x2.i,x2.r) > atan2(x1.i,x1.r)) { return 1; }
    else { return 0; }
}


#ifdef __cplusplus
}
}
#endif
