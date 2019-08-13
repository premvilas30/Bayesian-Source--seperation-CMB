//optimizer.h

#ifndef __OPTIMIZER_H__
#define __OPTIMIZER_H__

#include "defs.h"
#include "required_libs.h"
#include "block.h"
#include "patch.h"
#include "hyperpar.h"
#include "super.h"
#include "model_gsl.h"
#include "nelder_mead_simplex.h"
//#include "bfgs3.h"

struct optpar
{
	int npar;
	double *opt_theta;
	double f_opt_theta;
	double *grad_opt_theta;
	gsl_matrix *hessian;
};

struct optpar *optimizer_optpar_create( int nhyper );

void optimizer_optpar_destroy( struct optpar *opt );

struct optpar *optimizer( struct super *s, int max_iter, int thread_num, int restarts  );

struct optpar *optimizer_zzz( struct super *s, int max_iter, int thread_num, int restarts  );

struct optpar *optimizer2( struct super *s, int max_iter, int thread_num  );

gsl_matrix *optimizer_get_hessian( void *par, struct optpar *opt, double f );

#endif
