//simulate IGMRF

#include "simulateIGMRF.h"

void simulate_IGMRF( int n_side, double precision, double *x,  int iter )
{

	gsl_rng *r = simulate_IGMRF_setup_RNG() ;

	int k, l, t, n_pixel, nn, *mask;
	
	double sd = 1./sqrt( precision ) , mu; 
	
	n_pixel = (n_side) * (n_side) ;
	
	mask = calloc( n_pixel, sizeof(int) );
	
	for( k=0; k<n_pixel; k++ ) mask[k] = 1;
	
	struct graph  *gr ;  
	
	gr = graph_construct_exact_graph_ignoring_masked_pixels( n_side, n_side, mask ) ;
	
	//initialize 
	
	for( k=0; k<n_pixel; k++ ) x[k] = gsl_ran_gaussian( r, sd ) ;
	
	//sampling
	
	for( t=0; t<iter; t++ )
	{
	
		for( k=0 ; k<gr->number_nodes; k++ )
		{
			mu = 0.;
			
			nn = gr->number_neighbour_nodes[k] ; 
			
			//compute mean
			for( l=0; l<nn; l++ ) mu += x[ gr->neighbour_nodes[k][l] ]  ; 
			
			mu /= nn ;
			
			x[k] = mu + gsl_ran_gaussian( r, sd/sqrt(nn) );	
		}
	
	}
	
	graph_destroy( gr );

	free( mask );
	
	gsl_rng_free( r );
	
	return;

}


gsl_rng *simulate_IGMRF_setup_RNG()
{
	const gsl_rng_type *T;
	gsl_rng *r;
	
  	gsl_rng_env_setup();
  	T = gsl_rng_default;
  	r = gsl_rng_alloc(T);
  	
  	gsl_rng_set(r,(int)time(NULL));
  	
  	return(r);
 }
