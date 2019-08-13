
#ifndef __MODEL_GSL_H__
#define __MODEL_GSL_H__

#include "defs.h"
#include "required_libs.h"
#include "super.h"
#include "model.h"

struct model_gsl
{
	struct super *s;
	int thread_num;
};

double gsl_objective_function( const gsl_vector *x, void *params ) ;

void gsl_objective_function_gradient( const gsl_vector *x, void *params, gsl_vector *grad ) ;

void gsl_objective_function_and_gradient( const gsl_vector *x, void *params, double *f, gsl_vector *grad ) ;

double gsl_objective_function2( const gsl_vector *x, void *params ) ;

void gsl_objective_function_gradient2( const gsl_vector *x, void *params, gsl_vector *grad ) ;

void gsl_objective_function_and_gradient2( const gsl_vector *x, void *params, double *f, gsl_vector *grad ) ;

#endif 
