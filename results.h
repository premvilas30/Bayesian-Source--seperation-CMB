
#ifndef __RESULTS_H__
#define __RESULTS_H__

#include "fitsio.h"
#include "super.h"

struct results
{
	int nout_map;
	int nin_map; 
	int covariance;
	int residual;
	double **mu;
	double **sdev;
	double **resid;
	double **precision;
	double *sync_ind;
	double *dust_ind;
};

struct results *results_create( int n_pixel, int nin_map, int nout_map, int covariance, int residual );

void results_destroy( struct results *r );

void results_create_fits_results_files( struct block *b, int n_file, char **file );

void results_write_block_to_input_vector( int k_, int block_id, int residual, double *x, double *resid, struct patch **p, struct block *b, struct hyperpar *h );

void results_write_block_to_results(int block_id, struct results *r, struct patch **p, struct block *b, struct hyperpar *h, int simulated );

int results_write_patch_result_to_file( int patch_id, struct results *r, struct block *b, char **file );

int results_write_patch_result_to_file_0( int patch_id, struct results *r, struct block *b, char **file );

int results_write_patch_result_to_file_2( int patch_id, struct results *r, struct block *b, char **file );

void results_simulated_data_copy_hitrate( int npix, int nin_map, char **infile, char **outfile );

#endif 
