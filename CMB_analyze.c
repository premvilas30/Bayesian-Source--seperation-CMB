
#include "super.h"

void cmb_check( 	/*dimensions of the problem*/
						int *nside, int *nrow_subblock, int *ncol_subblock, int *nin_map, int *nout_map,
						/*observation maps*/
						int *in_type, int *mask, double *obs_intensity, double *hit_rate, char **file,
						/*model parameters*/
						double *obs_freq, double *obs_precision, double *src_precision, double *spec,
						/*output*/
						double *nz_B, double *nz_C, double *nz_Q__, double *nz_mu  )
{

	struct super *s = super_create( 0, *nside, *nside, *nrow_subblock, *nrow_subblock, *nin_map, *nout_map, mask) ;
	
	super_initialize( s, file, obs_freq, obs_precision, src_precision, spec, *in_type, obs_intensity, hit_rate ); 
	
	update_patch( s->p[0], s->h, FALSE ) ;
	update_patch( s->p[1], s->h, TRUE ) ;
	//patch_combine( s->p, s->b );
	
	double *x;
	int n, k;
	
	//non zero values in B
	x = (double *) s->p[1]->B->x ;
	n = s->p[1]->B->nzmax, k ;
	for( k=0; k<n ; k++ ) nz_B[ k ] = x[ k ] ;
	
	//non zero values in C
	x = (double *) s->p[1]->C->x ;
	n = s->p[1]->C->nzmax, k ;
	for( k=0; k<n ; k++ ) nz_C[ k ] = x[ k ] ;
	
	//non zero values in Q__
	x = (double *) s->p[1]->Q__->x ;
	n = s->p[1]->C->nzmax, k ;
	for( k=0; k<n ; k++ ) nz_C[ k ] = x[ k ] ;
	
	//non zero values in mu
	x = (double *) s->p[1]->mu->x ;
	n = s->p[1]->mu->nrow ;
	for( k=0; k<n ; k++ ) nz_mu[ k ] = x[ k ] ;
	
	super_destroy( s ) ;
	
	return;	
}
