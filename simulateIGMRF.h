
#ifndef __SIMULATE_IGMRF_H__
#define __SIMULATE_IGMRF_H__

#include "required_libs.h"
#include "defs.h"
#include "graphs.h"

void simulate_IGMRF( int n_side, double precision, double *x,  int iter ) ;

gsl_rng *simulate_IGMRF_setup_RNG();


#endif
