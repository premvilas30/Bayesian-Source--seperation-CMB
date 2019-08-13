//update.h

#ifndef __UPDATE_H__
#define __UPDATE_H__

#include "required_libs.h"
#include "defs.h"
#include "hyperpar.h"
#include "block.h"
#include "patch.h"


double update_sync_link( double x );

double update_gdust_link( double x );

double update_precision_link( double x );

void update_temperature_conversion_factors( struct patch *p, struct hyperpar *h );

double update_cmb_tdymc( double nu, double nu_ref, double beta );

double update_pow_law( double nu, double nu_ref, double beta );

void update_A_x( struct patch *p, struct hyperpar *h );

void update_A( struct patch *p, struct hyperpar *h );

void update_A_simulated( struct patch *p, struct hyperpar *h );

void update_B_pixel_specific( struct patch *p, struct hyperpar *h );

void update_F( struct patch *p, struct hyperpar *h );

void update_F_simulated( struct patch *p, struct hyperpar *h );

void update_B( struct patch *p, struct hyperpar *h );

void update_Q( struct patch *p, struct hyperpar *h );

void update_Q_using_mask( struct patch *p, struct hyperpar *h, int *mask );

void update_Q__( struct patch *p, struct hyperpar *h );

void update_L_initial( struct patch *p );

void update_mu( struct patch *p );

cholmod_dense *update_linear_predictor( struct patch *p, struct hyperpar *h );

void update_covariance( struct patch *p );

void update_patch( struct patch *p, struct hyperpar *h, int all, int *mask );

void update_initial( struct block *b, struct patch **p, struct hyperpar **h );

void update_simulate_data( struct patch *p, struct hyperpar *h, int *mask, int seed );

#endif
