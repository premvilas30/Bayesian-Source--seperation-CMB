
#ifndef __NELDER_MEAD_SIMPLEX__
#define __NELDER_MEAD_SIMPLEX__

#include "defs.h"
#include "required_libs.h"
#include "model_gsl.h"


void optimizer_nmmin(int n, double *Bvec, double *X, double *Fmin,
	   int *fail, double abstol, double intol, 
	   double alpha, double bet, double gamm, int trace,
	   int *fncount, int maxit, void *params );
	   

#endif
