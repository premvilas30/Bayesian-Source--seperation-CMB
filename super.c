//super.c : functions for the super struct holding all information

#include "super.h"

struct super *super_create( int *model, int id, int nrow, int ncol, int nrow_subblock, int ncol_subblock, int nin_map, int nin_template, int nout_map, int *templated, int *mask, int individual, int masked_pix_include, int nthread, int covariance, int residual, int n_specidx, int *specidx, int *parinfer, int *prior_mu, double *prior_mu_ref_freq, int *spectype )
{
	//id gives the base resolution id for the patch being analyzed going 0,...,11
	int k, Npar = nout_map +2 + (nout_map-1) + nin_map ;
	
	struct super *s = ( struct super *)malloc( sizeof( struct super ) ) ;
	
	s->id = id ;
	
	s->individual = individual ;
	
	s->individual_id = calloc( nthread, sizeof(int) );
	
	s->thread_num =  calloc( nthread, sizeof(int) ); 
	
	s->b = block_create( nrow, ncol, nrow_subblock, ncol_subblock, nin_map, nin_template, nout_map, mask, masked_pix_include ) ;
	
	s->p = patch_create_from_block( model, s->b, nin_map, nin_template, nout_map, templated, masked_pix_include, covariance, residual, n_specidx, specidx, parinfer, prior_mu, prior_mu_ref_freq, spectype ) ;
	
	s->h = malloc( ( s->b->n_block + 1 )  * sizeof( struct hyperpar *) ) ;
	
	int n_parinfer = 0;
	for( k=0; k< Npar ; k++ ) n_parinfer += parinfer[k] ; 
	
	for( k=0; k<s->b->n_block+1 ; k++) s->h[ k ] = hyperpar_create( nin_map, nout_map, n_parinfer ) ;
	
	return( s ) ;
}

void super_destroy( struct super *s )
{
	int k, n_block = s->b->n_block;
	
	free( s->individual_id );
	
	free( s->thread_num);
	
	if( !s->individual ) patch_destroy_from_block( s->p, s->b, s->individual ) ;
	
	block_destroy( s->b ) ;
	
	for( k=0; k< n_block+1; k++) hyperpar_destroy( s->h[k] );
	free( s->h );
	
	free( s );
	
	return;
}

void super_initialize(	struct super *s, char **files, char **template_files, double *obs_freq, double *obs_precision, double *offset,
								double *spec , int in_type, double *obs_intensity, double *hit_rate, int explo )
{
	int k;
	
	for( k=0; k< s->b->n_block+1; k++ )
	{
		hyperpar_initialize( s->h[ k ], obs_precision, spec, obs_freq, offset, explo, s->p[k]->template_source );
	}
	
	data_read_maps( s->id, s->b, s->p, s->h[0], files, in_type, obs_intensity, hit_rate, explo );
	
	//printf("\n Leaving data_read_maps");
	//printf("\n Leaving data_read_maps");
	
	//if( s->p[0]->nin_template > 0 )
	//{
	//	data_read_templates( s->id, s->b, s->p, template_files, explo );
	//}
	
	if( s->p[0]->ninferred_map > 0 )
	{
		data_read_templates( s->id, s->b, s->p, template_files, explo, TRUE );
	}
	
	if( s->p[0]->nin_template > 0 )
	{
		data_read_templates( s->id, s->b, s->p, template_files, explo, FALSE );
	}
	
	
	/*FILE *out = fopen( "yout.txt", "w");
	int k;
	for( k=0; k<s->b->n_unmasked[2] * s->p[0]->nin_map; k++ ) fprintf( out, "%.10f\n", ((double *) s->p[2]->y->x)[k] );
	fclose( out );*/

	//initialize the matrices and the Cholesky factors and the conversion factors
	update_initial( s->b, s->p, s->h ) ;
	
	if( s->individual ) patch_destroy( s->p[0], FALSE ) ; //leave the zero'th place free to hold arbitrary patch

	return;
}
