//simulate a cmb style dataset

#include "simulate.h"

struct simulation *simulation_create( int n_side, int n_freq, int n_src )
{
	int k, n_pixel = n_side * n_side;
	
	struct simulation *sim = (struct simulation *)malloc(sizeof(struct simulation));
	
	sim->n_side = n_side;
	sim->n_freq = n_freq;
	sim->n_src = n_src;
	
	sim->src = calloc( n_src, sizeof( double *) );
	for( k=0; k<n_src; k++ )  sim->src[k] = calloc( n_pixel, sizeof(double) );
	
	sim->obs = calloc( n_freq, sizeof( double *) );
	for( k=0; k<n_freq; k++ ) sim->obs[k] = calloc( n_pixel, sizeof(double) );
	
	return( sim );
}

void simulation_destroy( struct simulation *sim )
{
	int k;
	
	for( k=0; k<sim->n_src; k++ ) free( sim->src[k] );
	free( sim->src );
	
	for( k=0; k<sim->n_freq; k++ ) free( sim->obs[k] );
	free( sim->obs );
	
	free( sim );
}

void simulation_simulate_data( int n_side, int n_freq, int n_src, double *nu, double *specpar, double *obs_precisions, double *src_precisions, char **obs_files, char **src_files, char *info_file  ) 
{

	int k, l, n_pixel = n_side * n_side ;
	
	int *mask = calloc( n_pixel , sizeof(int) );
	for( k=0; k<n_pixel; k++ ) mask[k] = 1; 
	
	struct hyperpar *h = hyperpar_create( n_freq, n_src );
	hyperpar_initialize( h, obs_precisions, specpar, nu, 0 );
	
	struct simulation *sim = simulation_create( n_side, n_freq, n_src );
	
	//simulate the IGMRF components
	
	for( k=0; k<n_src; k++ ) simulate_IGMRF( n_side, src_precisions[k], sim->src[k], 500 ) ;
	
	struct super *s = super_create( 0, n_side, n_side, n_side, n_side, n_freq, n_src, mask, 0, 1, 1, 0, 0 ) ;
	
	cholmod_free_sparse( &s->p[1]->Q__ , s->p[1]->chol_comm );
	patch_destroy( s->p[1], 0 );
	
	//put the src's into super
	for( k=0; k<n_src; k++ )
	{
		for( l=0; l<n_pixel; l++ ) ((double *)s->p[0]->mu)[ k*n_pixel + l ] = sim->src[k][l] ; 
	}
	
	update_patch( s->p[0], h, 0, mask );
	
	//multiply s->p->B by s->p->mu to get the mean of y
	
	double 	*a = calloc(2,sizeof(double)), *b = calloc(2,sizeof(double)) ; 
	a[0] = 1.;
	
	cholmod_dense *w = cholmod_allocate_dense( n_pixel * n_freq, 1, n_pixel * n_freq, CHOLMOD_REAL, s->p[0]->chol_comm );
	
	cholmod_sdmult( s->p[0]->B, 0, a, b, s->p[0]->mu, w , s->p[0]->chol_comm );
	
	free( a );
	free( b );	
	
	//add noise to w and store in sim->obs
	
	gsl_rng *r = simulate_IGMRF_setup_RNG();
	
	for( k=0; k<n_freq; k++ )
	{
		for( l=0; l<n_pixel; l++ )
		{
			sim->obs[k][l] = ((double *)w->x)[ n_pixel * k + l ] + gsl_ran_gaussian( r , 1./sqrt( obs_precisions[k]) ) ;
		}
	}
	
	gsl_rng_free( r );
	cholmod_free_dense( &w, s->p[0]->chol_comm );
	
	//now write the observations and generated sources to file
	
	FILE *fout ;
	
	for( k=0; k<n_freq; k++ )
	{
		fout = fopen( obs_files[k], "w" );
		
		for( l=0; l<n_pixel; l++ ) fprintf( fout, "%lf\n", sim->obs[k][l] );
		
		fclose( fout );
	}
	
	for( k=0; k<n_src; k++ )
	{
		fout = fopen( src_files[k], "w" );
		
		for( l=0; l<n_pixel; l++ ) fprintf( fout, "%lf\n", sim->src[k][l] );
		
		fclose( fout );
	}
	
	fout = fopen( info_file, "w" );
	
	fprintf( fout, "n_side \t %d\nn_freq \t %d\nn_src \t %d", n_side, n_freq, n_src );
	fprintf( fout, "\nsync \t %lf\ndust \t %lf", specpar[0], specpar[1] );
	for( k=0; k<n_freq; k++ ) fprintf( fout, "\nobsprec[%d] \t %lf", k , obs_precisions[k] );
	for( k=0; k<n_src; k++ ) fprintf( fout, "\nsrcprec[%d] \t %lf", k, src_precisions[k]);
	
	fclose( fout );
	
	
	hyperpar_destroy( h );
	
	free( s->individual_id );
	
	free( s->thread_num);
	
	patch_destroy( s->p[0], 0 );
	
	for( k=0; k< 2; k++) hyperpar_destroy( s->h[k] );
	
	block_destroy( s->b ) ;
	
	simulation_destroy( sim );
	
	free( mask );
	
	return;
}






