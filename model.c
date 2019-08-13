//model.c

#include "model.h"

double model_log_igmrf_prior( cholmod_dense *x, struct patch *p, struct hyperpar *h )
{
	//evaluate the igmrf prior at x ( which is of length p->nout_map * p->graph->number_nodes )
	int nl = p->ninferred_map * p->graph->number_nodes, k;
	double 	*a = calloc(2,sizeof(double)),
				*b = calloc(2,sizeof(double));
	a[0] = 1.; 
	
	cholmod_dense *w = cholmod_zeros( nl, 1, CHOLMOD_REAL, p->chol_comm );
	cholmod_dense *v = cholmod_zeros( nl, 1, CHOLMOD_REAL, p->chol_comm );
	
	double 	*vv = (double *)v->x, *xx = (double *)x->x, 
			*mu = (double *)p->mu_prior->x, *ww = (double *)w->x ;
	
	for( k=0; k<nl; k++ ) 
	{
		vv[k] = xx[k] ;//- mu[k] ;
		//if( k<20) printf("\n value at pixel %d is  %.2f", k, mu[k]);	
	}
	
	//compute Q(x-mu)
	cholmod_sdmult( p->Q, 0, a, b, v, w , p->chol_comm );
	
	//cholmod_print_sparse(  p->Q, "Q", p->chol_comm );
	//cholmod_print_dense(  w, "w", p->chol_comm );
	
	//compute  (x-mu)^t [ Q (x-mu) ]
	double qform  = 0.;
	for( k=0; k<nl; k++ ) qform += vv[k] * ww[k];
	
	//printf("\n GMRF: the value of qform is %.10f ",qform);
	
	cholmod_free_dense( &w, p->chol_comm );
	cholmod_free_dense( &v, p->chol_comm );
	
	/*brief change determinant eval to nnode instead of nnodes-1*/
	
	//add the log normalizing constant (for the IGMRF)
	
	int N = 0; 
	for( k=0; k<p->ninferred_map; k++ )
	{
		N += p->graph_models[0]->number_nodes ;
		if( p->model[k] == 1 ) N -= 1 ;
	}
	
	
	//if( p->model == NNPIXEL ) N = p->graph->number_nodes - 1;
	//if( p->model == INDPIXEL ) N = p->graph->number_nodes ;
	
	double log_const = -.5 * N * log_2_pi;
	for( k=0; k<p->nout_map; k++ ) 
	{
		if( !p->template_source[k] )
		{
			log_const += .5 * ( p->graph_models[0]->number_nodes - p->model[k]  ) * h->theta[2+k];
		}
	}
	
	free(a);
	free(b);
	
	return( log_const - .5 * qform ); 
}

double model_log_likelihood( cholmod_dense *x, struct patch *p, struct hyperpar *h, int masked_pixel_include, int *mask )
{
	int nd = p->nin_map * p->graph->number_nodes, k, l;
	int nf = p->ninferred_map * p->graph->number_nodes;
	
	double 	*a = calloc(2,sizeof(double)),
				*b = calloc(2,sizeof(double)),
				qform = 0.;
	a[0] = 1.;
	
	//printf("\n In likelihood and nd = %d, ninmap = %d, numbernodes = %d\n",nd, p->nin_map, p->graph->number_nodes);
	
	//modify the field vector x to allow for offsets
	
	//offsets[0] = 0.; offsets[1] = .25 ;
	//offsets[2] = .05; offsets[3] = 1;
	
	/*double *xx = (double *)x->x ;
	
	for( k=1; k<p->nout_map; k++ )
	{
		for( l=0; l<p->graph->number_nodes; l++ ) xx[ k * p->graph->number_nodes +  l ] += h->theta[k+5] ; 
	}*/
	
	//cholmod_dense *f = cholmod_allocate_dense( nf, 1, nf, CHOLMOD_REAL, p->chol_comm );
	
	int include_field_mean =  TRUE, c = 0;
	
	/*for( k=0; k<p->ninferred_map; k++ )
	{
		
		for( l=0 ; l<p->graph->number_nodes; l++ ) ff[ k * p->graph->number_nodes + l ] = xx[ k * p->graph->number_nodes + l ] ;

	}*/
	
	//printf("\n value of nf = %d", nf);
	
	//cholmod_print_dense( f, "f", p->chol_comm );
	
	
	cholmod_dense *w = cholmod_zeros( nd, 1, CHOLMOD_REAL, p->chol_comm );
	cholmod_dense *z = cholmod_zeros( nd, 1,  CHOLMOD_REAL, p->chol_comm );
	//cholmod_dense *u = cholmod_zeros( nd, 1, CHOLMOD_REAL, p->chol_comm );
	
	double *wx = (double *) w->x;
	//double *ux = (double *) u->x;
	
	/*printf("\n ==== before ===\n");
	
	for( k=0; k<10; k++ ) printf( " \n w[%d] = %lf", k, wx[k] );
	
	printf("\n===== after =====\n");*/
	
	if( p->ninferred_map > 0 )
	{
		//printf("\n  in the B multiply part");
		cholmod_sdmult( p->B, 0, a, b, x, w , p->chol_comm );
	}
	
	/*for( k=0; k<10; k++ ) printf( " \n w[%d] = %lf", k, wx[k] );
	
	printf("\n===== after =====\n");*/

	if( p->nin_template > 0 )
	{
		//add G T to w 
		
		b[0] = 1.;
		
		cholmod_sdmult( p->G, 0, a, b, p->mu_template, w, p->chol_comm );
	}
	
	double *yx = (double *) p->y->x;
	//double *mu_t = (double *) p->mu_template->x;
	
	/*FILE *fo = fopen("../SIM_small_data/testy.txt", "w");

	for( l=0; l<p->graph->number_nodes; l++ )
	{
		for( k=0; k<p->nin_map; k++ )
		{
			fprintf( fo, "%lf\t", wx[ k * p->graph->number_nodes + l ] );
		}
		fprintf(fo,"\n");
	}
	
	fclose(fo);*/
	
	/*FILE *fo = fopen("../SIM_small_data/eta.txt", "w");
	
	for( k=0; k<nd; k++ ) fprintf( fo, "%.20f\n", ((double*)w->x)[k]);
	
	fclose(fo);*/
	
	/*fo = fopen("../SIM_small_data/mu_temp.txt", "w");
	
	for( k=0; k<p->graph->number_nodes * p->nin_template; k++ ) fprintf( fo, "%.20f\n", ((double*)p->mu_template->x)[k]);
	
	fclose(fo);*/
	
	
	//for( k=0; k<100; k++ ) printf("\n %lf \t %lf ", yx[k], mu_t[k]);

	//printf("\n Vector norm 2 : %lf ", model_check_vector_norm( nd, w ) ) ;
	
	//for( k=0; k<10; k++ ) printf( " \n w[%d] = %lf", k, wx[k] );
	
	//incorporate the offsets here
	/*for( k=0; k<p->nin_map; k++ )
	{
		for( l=0; l<p->graph->number_nodes; l++ ) wx[ k * p->graph->number_nodes + l ] += h->theta[k+9] ;
	}*/
	
	wx = (double *) w->x;
	//ux = (double *) u->x;
	
	//cholmod_sparse *I = cholmod_speye( nd , nd, CHOLMOD_REAL, p->chol_comm ); //alloc memory and initialize later
	
	//b[0] = -1.;
	
	//double *yx = (double *) p->y->x;
	
	//cholmod_sdmult( I, 0, a, b, p->y, w, p->chol_comm ); 
	
	for( k=0; k<nd; k++ ) wx[k] = yx[k] - wx[k] /*- ux[k]*/; 
	
	
	//for( k=0; k<nd; k++ ) wx[k] = yx[k] - wx[k];
	
	//printf("\n Vector norm 3 : %lf ", model_check_vector_norm( nd, w ) ) ;
	
	b[0] = 0.;
	
	cholmod_sdmult( p->C, 0 , a, b, w, z, p->chol_comm );
	
	//printf("\n Vector norm 4 : %lf ", model_check_vector_norm( nd, z ) ) ;

	//double *wx = (double*) w->x ;	
	double *zx = (double*)z->x ;
	
	//can insert a mask term here... 
	
	for( k=0; k<nd; k++ ) 
	{ 
		qform += zx[k] * wx[k] ;  //else qform += mask[k] * ( wz[k] * wx[k] ) ;  
	}
	
	//cholmod_free_sparse( &I, p->chol_comm );
	
	/*printf("\n **** LOGLIKE **** ");
	printf("\n nconst =  %.10f ", - .5 * nd * log_2_pi);
	printf("\n logdet =  %.10f ",  .5 * p->log_det_C);
	printf("\n qform =  %.10f ", - .5 * qform);*/
	
	/*for( k=1; k<p->nout_map; k++ )
	{
		for( l=0; l<p->graph->number_nodes; l++ ) xx[ k * p->graph->number_nodes +  l ] -= h->theta[k+5] ; 
	}*/
	
	cholmod_free_dense( &w, p->chol_comm );
	cholmod_free_dense( &z, p->chol_comm );
	//cholmod_free_dense( &u, p->chol_comm );
	//cholmod_free_dense( &f, p->chol_comm );
	//if( p->nin_template > 0 ) cholmod_free_dense( &u, p->chol_comm );
	
	free(a);
	free(b);
	
	//printf("\n LOG LIKELIHOOD: %.5f",  - .5 * nd * log_2_pi + .5 * p->log_det_C - .5 * qform ) ;
	
	//printf("\n Model log joint likelihood =  %lf", - .5 * nd * log_2_pi + .5 * p->log_det_C - .5 * qform );
	
		//printf("\n Leaving likelihood ");
	
	return( - .5 * nd * log_2_pi + .5 * p->log_det_C - .5 * qform );
	
}

double model_check_vector_norm( int n, cholmod_dense *v )
{
	double norm = 0.;
	int k;
	
	double *vv = (double *) v->x ;
	
	for( k=0; k<n; k++ ) norm += vv[k] * vv[k] ;

	return( sqrt( norm ) );
}



double model_log_conditional_igmrf_at_mode( struct patch *p )
{
	return(  -.5 * p->graph->number_nodes * p->ninferred_map * log_2_pi + .5 * p->log_det_Q__ ) ;
}


double model_log_prior_hyperpar( struct patch *p, struct hyperpar *h)
{
	//this is where the priors are
	
	int k;
	double log_prior = 0., ldens;
	
	//take uniform priors over admissable range on the two spectral parameters and include a jacobian 
	
	if( p->parinfer[0] )
	{
	
		//need to adjust these for the different types of spectral parameters
	
		double sync = update_sync_link( h->theta[0] ) ;
	
		//trial for new Jeffrey's prior
	
		log_prior += model_jeff_spec_sync_log_prior( sync, p, h ) 
				 	+ model_ugaussian_log_density( sync, h->prior_mean_spec[0], h->prior_sd_spec[0] ) 
						 + gsl_sf_log( SYNC_UPPER  - SYNC_LOWER ) + h->theta[0] - 2.*gsl_sf_log( 1. + gsl_sf_exp( h->theta[0] ) )  ; 
		
		//log_prior +=  /*model_ugaussian_log_density( sync, h->prior_mean_spec[0], h->prior_sd_spec[0] ) 
			//	+*/ gsl_sf_log( SYNC_UPPER  - SYNC_LOWER ) + h->theta[0] - 2. * gsl_sf_log( 1. + gsl_sf_exp( h->theta[0] ) ) ;
	
	}
	
	/*log_prior +=  model_ulaplace_log_density( sync, h->prior_mean_spec[0], h->prior_sd_spec[0] ) 
				+ gsl_sf_log( SYNC_UPPER  - SYNC_LOWER ) + h->theta[0] - 2.*gsl_sf_log( 1. + gsl_sf_exp( h->theta[0] ) ) ;*/
	
	if( p->parinfer[1] )
	{
	
	double dust = update_gdust_link( h->theta[1] );
	
	//trial for Jeffrey's  priormodel
	
	log_prior += model_jeff_spec_dust_log_prior( dust, p , h ) +
				 + model_ugaussian_log_density( dust, h->prior_mean_spec[1], h->prior_sd_spec[1] )
					 + gsl_sf_log( DUST_UPPER  - DUST_LOWER ) + h->theta[1] - 2.*gsl_sf_log( 1. + gsl_sf_exp( h->theta[1] ) )  ; 
	
	//log_prior +=  /*model_ugaussian_log_density( dust, h->prior_mean_spec[1], h->prior_sd_spec[1] )
		//	+*/ gsl_sf_log( DUST_UPPER - DUST_LOWER ) + h->theta[1] - 2. * gsl_sf_log( 1. + gsl_sf_exp( h->theta[1] ) ) ;
	
	/*log_prior +=  model_ulaplace_log_density( dust, h->prior_mean_spec[1], h->prior_sd_spec[1] )
				+ gsl_sf_log( DUST_UPPER - DUST_LOWER ) + h->theta[1] - 2.*gsl_sf_log( 1. + gsl_sf_exp( h->theta[1] ) ) ;*/
	
	}
	
	//prior for gmrf precision
	
	for( k=0; k<4; k++ ) 
	{
		//try the PC prior
		
		//ldens = - gsl_sf_log( 2. ) + gsl_sf_log( h->lambda[k+2] ) - 1.5 * h->theta[k+2] - h->lambda[k+2] * gsl_sf_exp( -.5 * h->theta[k+2] )  ; 
		
		if( p->parinfer[ k + 2 ] )//!p->template_source[k] ) 
		{
			//ldens = - gsl_sf_log( 2. ) + gsl_sf_log( h->lambda[k+2] ) - 1.5 * h->theta[k+2] - h->lambda[k+2] * gsl_sf_exp( -.5 * h->theta[k+2] )  ; 
		
			ldens =  h->shapes[k+2] * gsl_sf_log( h->rates[k+2] ) -  gsl_sf_lngamma( h->shapes[k+2] ) + ( h->shapes[k+2] - 1. )  *  h->theta[k+2] - h->rates[k+2] * gsl_sf_exp( h->theta[k+2] ) ; 
		
			log_prior += ldens + h->theta[ k+2 ] ;  /*jacobian*/
		
		}
		else
		{
			//assume a uniform prior on the  amplitudes
			
			//ldens = model_ugaussian_log_density( h->theta[k+2] , 0, 1.) ;
		
			//log_prior += ldens ; /*jacobian*/
		}
	}
	
	//priors for the amplitudes
	
	for( k=0; k<3; k++ )
	{
		if( p->parinfer[ k+6 ] ) 
		{
			ldens =  3. * gsl_sf_log( 3. ) -  gsl_sf_lngamma( 3. ) + ( 3. - 1. )  *  h->theta[k+6] - 3. * gsl_sf_exp( h->theta[k+6] ) ; 
			
			log_prior += ldens + h->theta[ k+6 ];
		}
	}
	
	//priors for the monopoles
	for( k=0; k<h->nin_map; k++ )
	{
		if( p->parinfer[ k+9 ] )
		{
			log_prior += model_ugaussian_log_density( h->theta[k+9], 0., h->prior_sd_monopole[k] ) ;
		}
	}
			
	
	//printf("\n Model log joint prior =  %lf", log_prior );		
	
	return( log_prior ) ;
}


double model_log_joint_marginal_hyperpar( struct patch *p0, struct patch **p , struct hyperpar *h, int npatches, int masked_pixel_include, int *mask, int individual_id )
{
	double log_lik = model_log_likelihood( p0->mu , p0, h, masked_pixel_include, mask ),
			 log_prior_hpar = model_log_prior_hyperpar( p0 , h ), log_igmrf_prior ;
			 
	if( p0->ninferred_map > 0 ) log_igmrf_prior = model_log_igmrf_prior( p0->mu, p0, h ); else log_igmrf_prior = 0.;
	
	double log_joint = log_lik + log_igmrf_prior + log_prior_hpar ;
	
	double log_conditional = 0.;
	int k;

	if( p0->ninferred_map > 0 )
	{

		if( npatches == 1 )
		{
			log_conditional = model_log_conditional_igmrf_at_mode( p[ individual_id ] );
		}
		else
		{
			for( k=1; k< npatches+1; k++) log_conditional += model_log_conditional_igmrf_at_mode( p[k] ) ;
		}
	
	}
	//printf("\n log like = %.20f ", log_lik);
	/*printf("\n log igmrf pr = %.20f ", log_igmrf_prior );
	printf( "\n log prior hpar = %.20f", log_prior_hpar );
	
	printf("\n Model log joint hyperpar =  %lf", log_joint - log_conditional );*/
	
	return( log_joint - log_conditional );
}

double model_ugaussian_log_density( double x, double mu, double sd )
{
	double z = ( x - mu ) / sd ;
	
	return( - gsl_sf_log( sd ) - log_2_pi - .5 * z * z ) ;
}

double model_ulaplace_log_density( double x, double mu, double scale )
{
	double z = fabs( x - mu ) / scale ;
	
	return( - log( 2. * scale )  - z ) ;
}

double model_jeff_spec_sync_log_prior( double sync, struct patch *p, struct hyperpar *h )
{
	int k;
	double s = 0.;
	
	//there a problem here- instead of nu[0] on the bottom it should be the reference frequency!
	
	for( k=0; k<h->nin_map; k++ ) s += gsl_pow_2( p->conv_ant_to_therm[k] ) * h->obs_precision[k] * pow( h->nu[k] / p->mu_prior_ref_freq[1] , 2.*sync ) * gsl_pow_2( log( h->nu[k] / p->mu_prior_ref_freq[1] ) ) ;
	
	return( .5 * log(s) );
}

double model_jeff_spec_dust_log_prior( double dust, struct patch *p, struct hyperpar *h )
{
	int k;
	double s = 0.;
	
	//there a problem here- instead of nu[0] on the bottom it should be the reference frequency!
	
	for( k=0; k<h->nin_map; k++ ) s += gsl_pow_2( p->conv_ant_to_therm[k] ) * h->obs_precision[k] * pow( h->nu[k] / p->mu_prior_ref_freq[2] , 2.*dust ) * gsl_pow_2( log( h->nu[k] / p->mu_prior_ref_freq[2] ) ) ;
	
	return( .5 * log(s) );
}



