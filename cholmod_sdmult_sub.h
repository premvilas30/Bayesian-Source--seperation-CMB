#ifndef __CHOLMOD_SDMULT_SUB_H__
#define __CHOLMOD_SDMULT_SUB_H__

#include "defs.h"
#include "required_libs.h"

int cholmod_sdmult_sub( cholmod_sparse *A, double *a, double *b, cholmod_dense *x, cholmod_dense *y);

#endif
