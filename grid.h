// grid.h

#ifndef __GRID_H__
#define __GRID_H__

#include "defs.h"
#include "required_libs.h"
#include "super.h"
#include "model_gsl.h"

#define MAX_GRID_PTS 10000

struct grid
{
	int npar; //col dimension of pts
	int max_pts; //maximum number of points in the grid
	int npts ; //number of points in the grid 0,...,npts-1
	double **pts ; //indexing in to each of the hyperparameter pts
	double *log_function_values ; //holds the values of the objective function
	double max_log_function_value; //hols the maximum value of the objective over the grid
	double *weights ; //normalized exp( log_function_values )
	double delta_z ; 
	double log_marginal_likelihood; //approximate log marginal likelihood
};

struct grid_axes
{
	int npar;
	gsl_vector *e_vals;
	gsl_matrix *e_vecs;
	gsl_vector *sqrt_e_vals;
};


struct grid *grid_create( int npar, int max_pts_par ) ;

void grid_destroy( struct grid *g ) ;

struct grid_axes *grid_create_grid_axes( int npar ) ;

void grid_destroy_grid_axes( struct grid_axes *ga ) ;

struct grid_axes *grid_construct_axes( struct hyperpar *h , struct patch *p ) ;

void grid_z_to_x( double *x, int npar, double *modal_x, double *z, struct grid_axes *ga );

struct grid *grid_construct( struct super *s, int thread_num ) ;

void grid_compute_posterior_functionals( struct grid *g, struct patch *p, struct hyperpar *h, int *mask, int compute_covariance, int compute_residual );

void grid_print( struct grid *g, FILE *fp ) ;

struct grid *grid_construct2( struct super *s, int thread_num );

#endif
