//hyperpar.h

#ifndef __HYPERPAR_H__
#define __HYPERPAR_H__

#include "required_libs.h"
#include "defs.h"
#include "precisions.h"

struct hyperpar
{
	int npar; //no. of unknown hyperparameters
	int nin_map; //no. of maps to be imput
	int nout_map; //no. of maps to be inferred
	double *nu; // frequency in units of GHz
	double *theta; //hyperparameters transformed
	double *prior_mean_spec; //prior mean on normal for spectral  params
	double *prior_sd_spec; //prior sd on normal for spectral params
	double *obs_precision; //known detector precision 
	//parameters for prior on precision
	double *rates; //gamma prior
	double *shapes;
	double *lambda; //PC prior
	double *prior_sd_monopole; //monopole sd prior
	double *offset; //offsets
	gsl_matrix *hessian; //stored hessian at the mode
	gsl_vector *e_vals; //stored eigenvalues (for computing marginal likelihood approximation)
};

struct hyperpar *hyperpar_create( int nin_map, int nout_map, int n_parinfer );

void hyperpar_destroy( struct hyperpar *h );

void hyperpar_initialize( struct hyperpar *h, double *obs_precision, double *spec, double *obs_freq, double *offset, int explo, int *templated );

#endif
