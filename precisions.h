//precisions.h

#ifndef __PRECISIONS_H__
#define __PRECISIONS_H__

#include "required_libs.h"
#include "defs.h"

//expected ranges of sources in milli Kelvin

#define RANGE_CMB .6//.01 //.3 //1E2
#define RANGE_SYNC .2//.04//1. //1.2E3
#define RANGE_DUST .6//.004//.05 //6E2
#define RANGE_FREE_FREE .3//.02//.5 //4E3

void precision_prior_gamma_get_rates_shapes( double *rates, double *shapes, int nout_map );

void precision_prior_gamma_get_rates_shapes_1( double *rates, double *shapes, int nout_map );

void precision_prior_PC_get_lambda( double *lambda );

void precision_initialize_log_precisions( double *theta );

#endif
