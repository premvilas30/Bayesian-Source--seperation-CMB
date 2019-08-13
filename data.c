 //data.c : I/O for fits with healpix orderings

//functions for input and output of data / results

#include "data.h"

int *data_get_mask_ones( int n_pixel )
{
	int k, *mask = calloc( n_pixel, sizeof(int) );
	
	for( k=0; k<n_pixel; k++ ) mask[k] = 1;
	
	return( mask );
}


int *data_get_mask( int id, int n_pixel, char *file_name)
{
	double *mask_d = calloc( n_pixel, sizeof(double) ), 
			  nulval = FITS_NULL_VAL;
	int k, hdu_num = 2, hdu_type, any_null, colno = 1, status = 0;
	fitsfile *fp;
	LONGLONG	n = (LONGLONG) n_pixel,
				elem = 1,
				row = id * n + 1 ;
	
	//printf("\n we are definitely in this function : npixel = %d ", n_pixel);
	//printf("\n we are definitely in this function");	
	//printf("\n the file name here is %c",  file_name[1]);
	//printf("\n the file name here is %c",  file_name[1]);

	if( fits_open_file( &fp, file_name, READONLY, &status ) ) printf( "\n Mask file: fits error: unable to open file " ) ;

	if( fits_movabs_hdu( fp, hdu_num, &hdu_type, &status ) ) printf( "\n Mask file: fits error: unable to move within fits file ") ;

	fits_read_col( fp, TDOUBLE, colno, row, elem, n, &nulval, mask_d, &any_null, &status) ;
	
	fits_close_file( fp, &status );
	
	int *mask = calloc( n_pixel, sizeof(int) );
	
	//FILE  *mask_chk = fopen( "../PLANCK_data/mask_chk.txt", "w");

	for( k=0; k<n_pixel; k++ ) 
	{
		mask[k] = (int) mask_d[k];
		//fprintf( mask_chk, "%lf\n", mask[k]);
		//printf("\n%d",mask[k]);
	}
	
	//fclose( mask_chk );
	
	free( mask_d );
	
	return( mask );
}

void data_read_maps( int id, struct block *b, struct patch **p, struct hyperpar *h, char **file, int in_type, double *obs_intensity, double *hit_rate, int explo )
{
	//different types of input (in_type) 0: fits 1: ascii txt 2: provided by args obs_intensity and hit_rate
	
	fitsfile *fp;
	
	int 	k, l, status,
			nin_map = h->nin_map, 
			hdu_num = 2,
			hdu_type,
			any_null;
			
	LONGLONG	n = (LONGLONG)b->n_pixel,
				elem = 1,
				row = id * n + 1; 
	
	FILE *fp_a, *fp_a_hit;
			
	double *y = calloc( b->n_pixel, sizeof(double) ),
			 *z = calloc( b->n_pixel , sizeof(double) ), //z plays the role of hitrate OR pixel noise variance
			 nulval = FITS_NULL_VAL ;
	
	for( k=0; k<nin_map; k++ )
	{
		status = 0;
		
		printf("\n reading map %d", k );
		
		switch( in_type )
		{
			
			case 0:
				if( fits_open_file( &fp, file[k], READONLY, &status ) ) printf( "\n fits error: unable to open file " ) ;
				if( fits_movabs_hdu( fp, hdu_num, &hdu_type, &status ) ) printf( "\n fits error: unable to move within fits file ") ;
				fits_read_col( fp, TDOUBLE, 1, row, elem, n, &nulval, y, &any_null, &status) ;
				if( !explo )
				{
					fits_read_col( fp, TDOUBLE, 2, row, elem, n, &nulval, z, &any_null, &status) ; //hitrate if WMAP
				}
				else
				{
					fits_read_col( fp, TDOUBLE, 3, row, elem, n, &nulval, z, &any_null, &status ) ; //pixel noise variance if PLANCK
				}
				fits_close_file( fp, &status );
			break;
			
			case 1:
				fp_a = fopen( file[k], "r" ) ;
				fp_a_hit = fopen( file[ nin_map ], "r") ;
				for( l=0; l<b->n_pixel; l++ )
				{
					fscanf( fp_a, "%lf", &y[l] ) ;
					fscanf( fp_a_hit, "%lf", &z[l] ) ;
				}
				fclose( fp_a );
				fclose( fp_a_hit );
			break;
			
			case 2:
				for( l=0; l<b->n_pixel; l++ )
				{
					y[ l ] = obs_intensity[ k * b->n_pixel + l ] ;
					z[ l ] = z[ k * b->n_pixel + l ] ;
				}
			break;
			
		}
		
		if( explo )
		{
			//apply a conversion to units of mK for Planck data
			for( l=0; l<b->n_pixel; l++ )
			{
				y[l] *= 1E3 ;
				z[l] *= 1E6 ;
			}
		}
		
		//add an offset to each element of y
		
		for( l=0 ; l<b->n_pixel; l++ )
		{
			y[l] += h->offset[k] ; 
		}
		
		
		data_healpix_reorder( b, y, z, in_type, 0 );
		data_put_values( b, p, h, k, y, z, explo ) ;
		
		//put bit in here to read the maps of spectral parameters
		
		if( p[0]->n_specidx > 0 )
		{
			if( p[0]->specidx[0] )
			{
				if( fits_open_file( &fp, file[k+1], READONLY, &status ) ) printf( "\n fits error: unable to open file " ) ;
				if( fits_movabs_hdu( fp, hdu_num, &hdu_type, &status ) ) printf( "\n fits error: unable to move within fits file ") ;
				fits_read_col( fp, TDOUBLE, 1, row, elem, n, &nulval, y, &any_null, &status) ;
				fits_close_file( fp, &status );
				
				data_healpix_reorder( b, y, z, in_type, 0 );
				data_put_spec_values( b, p, y, 1);
				
			}		
		
		}
		
	}
	
	data_compute_det_C(  b, p ); 
	
	free( y );
	free( z );
	
	return;
}


void data_read_templates( int id, struct block *b, struct patch **p, char **template_file, int explo, int prior )
{
	fitsfile *fp;
	
	int 	k, l, c =  0, status, *g,
			nin_map = p[0]->nin_template, 
			hdu_num = 2,
			hdu_type,
			any_null;
			
	LONGLONG	n = (LONGLONG)b->n_pixel,
				elem = 1,
				row = id * n + 1; 
	
	FILE *fp_a, *fp_a_hit;
			
	double *t = calloc( b->n_pixel, sizeof(double) ),
		   *z = calloc( b->n_pixel, sizeof(double) ),
			 nulval = FITS_NULL_VAL, a_nu, arg ;
	
	if( prior ) g = p[0]->prior_mu_source ; else g = p[0]->template_source ;
	
	for( k=0; k<p[0]->nout_map; k++ )
	{ 
	
		if( g[k] )
		{
		
			//printf("\n template name : %s ", template_file[k] );
		
			status = 0;
		
			if( fits_open_file( &fp, template_file[k], READONLY, &status ) ) printf( "\n fits error: unable to open file : status = %d",status ) ;
			
			if( fits_movabs_hdu( fp, hdu_num, &hdu_type, &status ) ) printf( "\n fits error: unable to move within fits file ") ;
			fits_read_col( fp, TDOUBLE, 1, row, elem, n, &nulval, t, &any_null, &status) ;
			fits_close_file( fp, &status );
			
			//conversion for the Haslam map from WMAP explainer (Bennett et al 2003) pg 15 to mK Antenna temperature
			if( k == 1 && !UPDATE_SIMULATED )
			{
				for( l=0; l<b->n_pixel; l++ ) t[l] = 1000. * ( t[l] - 5.9) * pow( 23/.408, -2.9 ) ; 
			}
			
			//From LAMBDA
			//Dust emission "Model 8" from Finkbeiner et.al. (FDS 1999) was evaluated at 
			//94 GHz and used as the basis for the dust intensity templates. This 94GHz dust map, smoothed to 1 degree, 
			//was used as the dust emission template for Stokes I/temperature analysis
			

			if( k==2 )
			{
				//don't do anything with the dust template 
			}
			
			
			/*
			For free-free emission, we estimate the prior, P ff (p), using the extinction-– 15 –
corrected Hα map (Finkbeiner 2003). This is converted to a free-free signal using a conversion
factor of 11.4 μK R −1 (units of antenna temperature at K-band). From 'Three-Year Wilkinson Microwave Anisotropy Probe (WMAP 1 )
Observations:
Temperature Analysis' Bennett et al
			*/ 
			if( k == 3 && !UPDATE_SIMULATED )
			{
				//arg = pow( 10., log10( p[0]->mu_prior_ref_freq[k] ) + 9. + log10PLANCK - log10BOLTZMANN - log10T0 );
				//a_nu = exp( -  ( 2.*log( arg ) + arg - 2.* log( gsl_expm1(arg) ) ) ) ; 
				for( l=0; l<b->n_pixel; l++ )
				{
					t[l] *= .0114; // * a_nu ;
					//printf("\n t[%d] = %lf", l, t[l]);
				}
			}
			
			//for( l=0; l<10; l++ ) printf("\n value t[%d] = %lf ", l, t[l] );
			
			data_healpix_reorder( b, t, z, 0, 0 );
			
			//for( l=0; l<10; l++ ) printf("\n value t[%d] = %lf ", l, t[l] );
			
			data_put_template( b, p, c, t, prior ) ; 
			
			c++;
		
		}
	}	
	
	free( t );
	free( z );

}




void data_healpix_reorder( struct block *b, double *y, double *n_hit, int in_type, int ordering )
{
	//ordering gives the type of ordering 0 for nested, 1 for ring

	if( in_type > 0 ) return; //don't do anything if in ascii format
	
	//make copy of y and n_hit to re-order
	double 	*y_sort = calloc( b->n_pixel, sizeof(double) ),
				*n_hit_sort = calloc( b->n_pixel, sizeof(double) );
	
	int k, l;
	
	for( k=0; k<b->n_pixel; k++ ) 
	{
		if( ordering == 0 ) 
			l = b->healpix_2d_lattice_to_idx_nested[ k ] ; 
		else 
			l = b->healpix_2d_lattice_to_idx_ring[ k ] ;
			
		y_sort[ k ] = y[ l ] ;
		n_hit_sort[ k ] = n_hit[ l ] ;
	}
	
	FILE *rey,  *ren;
	//rey = fopen("../WMAP_data/test_of_remap_y.txt", "w");
	//ren = fopen("../WMAP_data/test_of_remap_n.txt", "w");
	
	for( k=0; k<b->n_pixel; k++ )
	{
		y[ k ] = y_sort[ k ] ;
		n_hit[ k ] = n_hit_sort[ k ] ;
		
		//if( k < 10) printf("\n entry %d = %lf ", k, y[k]);
		
		//fprintf( rey, "%lf\n", y[k] );
		//fprintf( ren, "%lf\n", (double) n_hit[k] );
	}
	
	//fclose(rey);
	//fclose(ren);
	
	free( y_sort );
	free( n_hit_sort );
	
	return;
}

void data_put_values( struct block *b, struct patch **p, struct hyperpar *h, int map, double *y, double *z, int explo )
{
	//y is a vector of length b->n_pixel-- direct mapping to patches
	//	can be done for patch
	
	int k, l, counter, masked_pix_inlcude = p[0]->masked_pix_include, *mask;
	
	FILE  *C_file;
	
	double *y_p, *C_diag; //C_diag is the precision of 
	for( k=0; k<b->n_block+1; k++ )
	{
		mask = b->mask[k] ;
		
		y_p = (double *) p[k]->y->x ;
		C_diag = (double *)p[k]->C->x;
		
		//C_file = fopen("../PLANCK_data/C_check.txt", "w");
		
		for( l=0; l<b->n_unmasked[k]; l++ )
		{
			y_p[ map * b->n_unmasked[k] + l ] = y[ b->idx[k][l] ] ;
			if( mask[l] ) 
			{
				if( !explo )
				{
					C_diag[ map * b->n_unmasked[k] + l ] = z[ b->idx[k][l] ] * h->obs_precision[ map ];
					//printf("\n val is  %lf", z[ b->idx[k][l] ] * h->obs_precision[ map ]);
				}
				else
				{
					C_diag[ map * b->n_unmasked[k] + l ] = z[ b->idx[k][l] ];
				//	fprintf( C_file, "%.5f\n", z[ b->idx[k][l] ] );
				}
			}
			else
			{
				C_diag[ map * b->n_unmasked[k] + l ] = 1E-6;
			}
		}
		
		//fclose( C_file  );
		
	}
	
	return;
}

void data_put_spec_values(  struct block *b, struct patch **p, double *v, int synch)
{

	int k, l, counter, masked_pix_inlcude = p[0]->masked_pix_include, *mask;
	
	double *v_p ; //C_diag is the precision of 
	for( k=0; k<b->n_block+1; k++ )
	{
		mask = b->mask[k] ;
		
		if( synch ) v_p = (double *) p[k]->synch_index->x ; else v_p = (double *) p[k]->dust_index->x ;
		
		for( l=0; l<b->n_unmasked[k]; l++ ) v_p[ l ] = v[ b->idx[k][l] ] ;
	}
	
	return;
}

void data_put_template( struct block *b, struct patch **p, int map, double *t, int prior)
{
	//y is a vector of length b->n_pixel-- direct mapping to patches
	//	can be done for patch
	
	int k, l, counter, masked_pix_inlcude = p[0]->masked_pix_include, *mask;
	
	double *t_p, *t_p1; 
	
	for( k=0; k<b->n_block+1; k++ )
	{
		mask = b->mask[k] ;
		
		if( prior) 
			t_p = (double *) p[k]->mu_prior->x ; 
		else 
			t_p = (double *) p[k]->mu_template->x ; //p[k]->mu_prior->x ;
		
		//t_p1 = (double *) p[k]->mu->x ;	
		
		for( l=0; l<b->n_unmasked[k]; l++ )
		{
			t_p[ map * b->n_unmasked[k] + l ] = t[ b->idx[k][l] ] ;
			//printf("\n  value[%d] = %lf ",map * b->n_unmasked[k] + l, t_p[ map * b->n_unmasked[k] + l ] );
			//t_p1[ map * b->n_unmasked[k] + l ] = t[ b->idx[k][l] ] ;
		}
		
		//if( k==0 ) cholmod_print_dense(  p[k]->mu_prior, "mu", p[k]->chol_comm );
	}
	
	return;
}


void data_compute_det_C( struct block *b, struct patch **p )
{
	int k, l, m, masked_pix_include = p[0]->masked_pix_include, *mask;
	double log_det_C;
	
	for( k=0; k<b->n_block; k++ )
	{
		log_det_C = 0.;
		mask = b->mask[k] ;
		for( m=0; m<p[k]->nin_map; m++ )
		{
			for( l=0; l<b->n_unmasked[k] ; l++ ) 
			{
				if( TRUE )
					log_det_C += log( ( (double *)p[k]->C->x )[ l + m * b->n_unmasked[k] ] );
				else
					log_det_C += log( ( (double *)p[k]->C->x )[ l + m * b->n_unmasked[k] ] ) * mask[ l ] ;
			}
		}
		p[k]->log_det_C = log_det_C;
		//printf("\n Block %d: n_unmasked = %d", k, b->n_unmasked[k] ) ; 
		//printf("\n Block %d: log|C| = %.2f", k , log_det_C);
	}
	
}


void data_write_out_patch_mean( int patch_idx, int nout_map, struct patch *p, struct block *b, char *file_name )
{
	int k, m, m_, masked_pix_include = p->masked_pix_include, n_pixel = b->nrow_subblock * b->ncol_subblock,
	n_unmasked = b->n_unmasked[ patch_idx ], n_masked = n_pixel - n_unmasked ;
	
	double *mu = calloc( nout_map * n_pixel, sizeof( double ));
	
	int *idx = calloc( n_unmasked , sizeof( int ));
	//int *idx_ = calloc( n_masked, sizeof(int) );
	
	m = 0;
	m_ = 0 ;
	for( k=0; k<n_pixel; k++ )
	{
		if( b->mask[ patch_idx ][k] )
		{
			idx[ m  ] = k;
			m++;
		}else{
			//idx_[ m_ ] = k;
			//m_++;
		}
	}
	
	for( m=0; m<nout_map; m++ )
	{
		for( k=0; k<n_unmasked; k++ ) 
		{
			mu[ m * n_pixel + idx[ k ]  ] = ((double *)p->mu->x)[ m * n_unmasked + k ];
			//if( m==0 ) printf("\n Put %d in slot %d",m * n_unmasked + k, m * b->n_pixel + b->idx[ patch_idx ][ k ] );
		}		
		// if there is masked pixels replace these with the prior...
		
		
	}
	
	FILE *fp;
	fp = fopen( file_name, "w" );
	
	for( k=0; k<nout_map*n_pixel; k++ ) fprintf( fp, "%lf\n", mu[k] );
	
	fclose( fp );
	
	free( mu );
	
	free( idx );
	//free( idx_ );
	
	return;

}



