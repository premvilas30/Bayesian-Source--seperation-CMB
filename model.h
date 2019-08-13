//model.h

#ifndef __MODEL_H__
#define __MODEL_H__

#include "defs.h"
#include "required_libs.h"
#include "super.h"

double model_log_igmrf_prior( cholmod_dense *x, struct patch *p , struct hyperpar *h) ;

double model_log_likelihood( cholmod_dense *x, struct patch *p, struct hyperpar *h, int masked_pixel_include, int *mask );

double model_log_conditional_igmrf_at_mode( struct patch *p ) ;

double model_log_joint_marginal_hyperpar( struct patch *p0, struct patch **p , struct hyperpar *h, int npatches, int masked_pixel_include, int *mask, int individual_id );

double model_log_prior_hyperpar( struct patch *p, struct hyperpar *h) ;

double model_check_vector_norm( int n, cholmod_dense *v );

double model_ugaussian_log_density( double x, double mu, double sd ) ;

double model_ulaplace_log_density( double x, double mu, double scale ) ;

double model_jeff_spec_sync_log_prior( double sync, struct patch *p, struct hyperpar *h );

double model_jeff_spec_dust_log_prior( double dust, struct patch *p, struct hyperpar *h );

#endif
