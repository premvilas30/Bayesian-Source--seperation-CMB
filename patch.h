//patch.h

#ifndef __PATCH_H__
#define __PATCH_H__

#include "defs.h"
#include "required_libs.h"
#include "routines.h"
#include "graphs.h"
#include "block.h"

struct patch
{
	//structure to hold information for 
	// a patch
	int *model;
	
	int masked_pix_include;
	
	int nin_map;

	int nout_map;
	
	int nin_template;
	
	int nin_mu_prior;
	
	int ninferred_map;
	
	int n_specidx;
	
	int covariance; //store the covariance  matrix?
	
	int residual; //compute residuals for the fit?
	
	int *specidx;
	
	int *spectype;
	
	int n_parinfer;
	
	int *parinfer; //indicator whether to infer parameters (TRUE) or fix (FALSE)
	
	struct graph *graph; //conditional independence graph
	
	struct graph **graph_models; //holds 2 graphs (independent pixels and nn dependencies)
	
	int *template_source; //indicator vector indicating whether the source is templated
	
	double *conv_ant_to_therm ; // conversion factors from antenna to thermodynamic temperature
	
	double **A; //loadings matrix for sources without templates
	
	double **F; //loadings matrix for sources with templates
	
	cholmod_sparse *B; //kronecker product of loadings with identity
	
	cholmod_sparse *I1; //kroncker product of column vector of ones with identity
	
	cholmod_sparse *G; //kronecker product of loadings (templates) with identity
	
	cholmod_sparse *C; //known observation level noise
	
	double log_det_C; //determinant of C computed from noise and hitrate map
	
	cholmod_sparse *T; //pattern for GMRF precision
	
	cholmod_sparse *Q; //numerical GMRF precision
	
	cholmod_sparse *Q__; //numerical posterior precision
	
	double log_det_Q__; //log of the determinant of Q
	
	cholmod_factor  *L; //symbolic Cholesky
	
	cholmod_dense *y; //observed intensity
	
	cholmod_dense *z; //for calculation of the modal value
	
	cholmod_dense *z1; //for further calculations
	
	cholmod_dense *mu; //mean of GMRF given y
	
	cholmod_dense *synch_index; //spectral index for synchrotron
	
	cholmod_dense *dust_index; //spectral index for dust
	
	cholmod_dense *mu_template; //the 'actual' GMRF for any templated sources
	
	cholmod_dense *mu_prior; //the prior mean for the GMRF
	
	double *mu_prior_ref_freq; //the reference frequency for the prior template
	
	double *template_ref_freq; //the reference frequency for the template  
	
	int *prior_mu_source; //indicator to say whether prior mean included or not
	
	//cholmod_dense *prior_mu; //the prior mean of the field
	
	cholmod_dense *Sig; //covariance matrix of GMRF given y
	
	cholmod_dense *sdev; //sqrt diagonal of covariance matrix of GMRF given y
	
	cholmod_dense *resid; //values of the residual for the fitted model
	
	cholmod_common *chol_comm; //cholmod common workspace :: necessary for multi threading
	
};

struct patch * patch_create( int *model, int nrow, int ncol, int nin_map, int nin_template, int nout_map, int *templated, int *mask, int all, int masked_pix_include,  int covariance, int residual, int n_specidx, int *specidx, int *parinfer, int *prior_mu, double *prior_mu_ref_freq, int *spectype );

void patch_destroy( struct patch *p, int all );

cholmod_sparse * patch_set_up_B( struct graph *graph , int  nin_map, int nout_map, cholmod_common *comm_ptr );

struct patch **patch_create_from_block( int *model, struct block *b, int nin_map, int nin_template, int nout_map, int *templated, int masked_pix_include, int covariance, int residual, int n_specidx, int *specidx, int *parinfer, int *prior_mu, double *prior_mu_ref_freq, int *spectype  );

void patch_destroy_from_block( struct patch **p, struct block *b, int individual );

void patch_combine( struct patch **p, struct block *b );

#endif

