//hyperpar.c

#include "hyperpar.h"

struct hyperpar *hyperpar_create( int nin_map, int nout_map, int n_parinfer )
{
	struct hyperpar *h = ( struct hyperpar *)malloc( sizeof( struct hyperpar ) ) ;
	
	h->npar = nout_map + 2 +  (nout_map  - 1 ) + nin_map ; //+ 3;
	h->nin_map = nin_map;
	h->nout_map = nout_map;
	//h->nu = calloc( nin_map, sizeof(double) );
	h->nu = calloc( nin_map, sizeof(double) );
	h->theta = calloc( h->npar , sizeof(double) );
	h->prior_mean_spec = calloc( 2, sizeof(double) );
	h->prior_sd_spec = calloc( 2, sizeof(double) );
	h->obs_precision = calloc( nin_map, sizeof(double) );
	h->rates = calloc( 2 + nout_map, sizeof(double) );
	h->shapes = calloc( 2 + nout_map, sizeof(double) );
	h->lambda =  calloc( 2 + nout_map, sizeof(double) );
	h->prior_sd_monopole = calloc( nin_map, sizeof(double) );
	h->offset = calloc( nin_map, sizeof(double) );
	h->hessian = gsl_matrix_alloc( (size_t) n_parinfer, (size_t) n_parinfer ) ;
	h->e_vals = gsl_vector_alloc( (size_t) n_parinfer ) ;
	return( h );
}

void hyperpar_destroy( struct hyperpar *h )
{
	free( h->theta );
	free( h->prior_mean_spec );
	free( h->prior_sd_spec );
	free( h->obs_precision );
	free( h->nu );
	free( h->rates );
	free( h->shapes );
	free( h->lambda ); 
	free( h->prior_sd_monopole ) ;
	free( h->offset );
	gsl_matrix_free( h->hessian );
	free( h );
	return;
}

void hyperpar_initialize( struct hyperpar *h, double *obs_precision, double *spec, double *obs_freq, double *offset, int explo, int *templated )
{
	//transform the spectral parameters to the real line	
	
	//printf("\n From hyperpar_intialize: %.2f, %.2f \n ", h->theta[0], h->theta[1] );
	
	//need to adjust these for the different types of spectral parameters
	
	h->prior_mean_spec[0] = -2.5;//( SYNC_UPPER - SYNC_LOWER )/2 + SYNC_LOWER ;
	
	h->prior_mean_spec[1] = 1.5;//( DUST_UPPER - DUST_LOWER )/2 + DUST_LOWER ;
	
	h->prior_sd_spec[0] = .15;
	
	h->prior_sd_spec[1] = .15;
	
	h->theta[0] = log( ( spec[0] - SYNC_LOWER ) / ( SYNC_UPPER - spec[0] ) ) ;
	
	h->theta[1] = log( ( spec[1] - DUST_LOWER ) / ( DUST_UPPER - spec[1] ) ) ; 
	
	int k;
	
	for( k=0; k<h->nin_map; k++ ) h->prior_sd_monopole[k] = 0.001;
	
	//put in rates for the gamma prior on the precisions
	//precision_prior_gamma_get_rates_shapes( h->rates, h->shapes , h->nout_map );
	precision_prior_gamma_get_rates_shapes_1( h->rates, h->shapes , h->nout_map );
	
	if( explo ) //make adjustment for Planck
	{
		//for(k=2;k<6;k++) h->rates[k] *= 1E-6;
		for( k=0; k<h->nin_map; k++ ) h->prior_sd_monopole[k] *= 1E-3;
	}
	precision_prior_PC_get_lambda( h->lambda );  
	
	for( k=2; k<6; k++ ) 
	{
		if( !templated[k-2] )
			h->theta[k] = log( h->shapes[k] / h->rates[k]  ); 
		else 
			h->theta[k] = 0.;
	}
	
	
	//free amplitude for sync, dust  and ffem
	for( k=6; k<9; k++ ) h->theta[k] = 0.;
	
	//monopoles
	for( k=9; k<14; k++ ) h->theta[k] = 0.;
	
	if( explo )  
	{
		//printf("\n Initial: \n");
		//for( k=0; k<6; k++ ) printf("\t %.2f", h->theta[k]); 
	}
	//put the precisions of the observations in
	for( k=0; k<h->nin_map; k++ ) h->obs_precision[ k ] = obs_precision[ k ] ; 
	
	for( k=0; k<h->nin_map; k++ ) h->offset[k] = offset[k] ;
	
	for( k=0; k<h->nin_map; k++ ) h->nu[ k ] = obs_freq[ k ];
	
	return;
}

