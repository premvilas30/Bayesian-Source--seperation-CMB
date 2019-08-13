
#ifndef __SUPER_H__
#define __SUPER_H__

#include "block.h"
#include "patch.h"
#include "hyperpar.h"
#include "data.h"
#include "update.h"

struct super
{
	int id;
	int explo; // which space exploration: planck == 1, else wmap == 0
	int individual; //analyze blocks of observations individually
	int *individual_id;
	int *thread_num;
	struct patch **p;
	struct block *b;
	struct hyperpar **h;
};

struct super *super_create( int *model, int id, int nrow, int ncol, int nrow_subblock, int ncol_subblock, int nin_map, int nin_template, int nout_map, int *templated, int *mask, int individual, int masked_pix_include, int nthread, int covariance, int residual, int n_specidx, int *specidx, int *parinfer, int *prior_mu, double *prior_mu_ref_freq, int *spectype );

void super_destroy( struct super *s );

void super_initialize(	struct super *s, char **files, char **template_files, double *obs_freq, double *obs_precision, double *offset,
								double *spec , int in_type, double *obs_intensity, double *hit_rate, int explo );

#endif
