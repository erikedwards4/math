//Sorti_Help functions.
//Give comparison (cmp) functions using structs that include the index (cmpi),
//where the index is represented as a size_t for faster speed in ranks.
//See ranks.c for usage.

#pragma once

#include <math.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int cmpi_ascend_s (const void *a, const void *b)
{
    const FLT_I x1 = *(const FLT_I *)a;
    const FLT_I x2 = *(const FLT_I *)b;
    if (x1.val>x2.val) { return 1; }
    else if (x2.val>x1.val) { return -1; }
    else { return 0; }
}


int cmpi_ascend_d (const void *a, const void *b)
{
	const DBL_I x1 = *(const DBL_I *)a;
    const DBL_I x2 = *(const DBL_I *)b;
    if (x1.val>x2.val) { return 1; }
    else if (x2.val>x1.val) { return -1; }
    else { return 0; }
}


int cmpi_ascend_c (const void *a, const void *b)
{
    const CFLT_I x1 = *(const CFLT_I *)a;
    const CFLT_I x2 = *(const CFLT_I *)b;
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
	const CDBL_I x1 = *(const CDBL_I *)a;
    const CDBL_I x2 = *(const CDBL_I *)b;
	const double sq1 = x1.r*x1.r + x1.i*x1.i;
    const double sq2 = x2.r*x2.r + x2.i*x2.i;
	if (sq1>sq2) { return 1; }
    else if (sq2>sq1) { return -1; }
	else if (atan2(x1.i,x1.r) > atan2(x2.i,x2.r)) { return 1; }
    else if (atan2(x2.i,x2.r) > atan2(x1.i,x1.r)) { return -1; }
    else { return 0; }
}


int cmpi_descend_s (const void *a, const void *b)
{
    const FLT_I x1 = *(const FLT_I *)a;
    const FLT_I x2 = *(const FLT_I *)b;
    if (x1.val>x2.val) { return -1; }
    else if (x2.val>x1.val) { return 1; }
    else { return 0; }
}


int cmpi_descend_d (const void *a, const void *b)
{
	const DBL_I x1 = *(const DBL_I *)a;
    const DBL_I x2 = *(const DBL_I *)b;
    if (x1.val>x2.val) { return -1; }
    else if (x2.val>x1.val) { return 1; }
    else { return 0; }
}


int cmpi_descend_c (const void *a, const void *b)
{
    const CFLT_I x1 = *(const CFLT_I *)a;
    const CFLT_I x2 = *(const CFLT_I *)b;
	const float sq1 = x1.r*x1.r + x1.i*x1.i;
    const float sq2 = x2.r*x2.r + x2.i*x2.i;
	if (sq1>sq2) { return -1; }
    else if (sq2>sq1) { return 1; }
	else if (atan2f(x1.i,x1.r) > atan2f(x2.i,x2.r)) { return -1; }
    else if (atan2f(x2.i,x2.r) > atan2f(x1.i,x1.r)) { return 1; }
    else { return 0; }
}


int cmpi_descend_z (const void *a, const void *b)
{
	const CDBL_I x1 = *(const CDBL_I *)a;
    const CDBL_I x2 = *(const CDBL_I *)b;
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
