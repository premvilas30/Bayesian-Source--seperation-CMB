//results.c : write out results

#include "results.h"

struct results *results_create( int n_pixel, int nin_map, int nout_map, int covariance, int residual )
{
	int k ;
	
	struct results *r = (struct results *)malloc( sizeof( struct results ));
	
	r->nout_map = nout_map;
	r->nin_map = nin_map ;
	r->covariance = covariance ;
	r->residual = residual ;
	r->mu = calloc( nout_map, sizeof( double *) ) ; 
	for( k=0; k<nout_map; k++ ) r->mu[k] = calloc( n_pixel, sizeof(double) );
	if( covariance )
	{
		r->sdev = calloc( nout_map, sizeof(double *) );
		for( k=0; k<nout_map; k++ ) r->sdev[k] = calloc( n_pixel, sizeof(double) );
	}
	if( residual )
	{
		r->resid = calloc( nin_map, sizeof(double *) );
		for( k=0; k<nin_map; k++ ) r->resid[k] = calloc( n_pixel, sizeof(double) );
	}
	r->precision = calloc( nout_map, sizeof(double *) );
	for( k=0; k<nout_map; k++ ) r->precision[k] = calloc( n_pixel, sizeof(double) );
	
	r->sync_ind = calloc( n_pixel, sizeof(double) );
	r->dust_ind = calloc( n_pixel, sizeof(double) );
	
	return( r );
}

void results_destroy( struct results *r )
{
	int k;
	for( k=0; k<r->nout_map; k++ ) 
	{
		free( r->mu[k] );
		free( r->precision[k] );
		if( r->covariance ) free( r->sdev[k] ) ;
	}
	for( k=0; k<r->nin_map; k++ )
	{
		if( r->residual ) free( r->resid[k] ) ;
	}
	free( r->mu );
	if( r->covariance ) free( r->sdev ) ;
	if( r->residual ) free( r->resid );
	free( r->precision );
	free( r->sync_ind );
	free( r->dust_ind );
	return; 
}


void results_create_fits_results_files( struct block *b, int n_file, char **file )
{
	fitsfile *fptr;
	int status=0, anynul, hdunum, hdutype;
	double nulval;
	int naxis2 = 12 * b->n_pixel  ;

	char **ttype;
	char **tform;
	char **tunit;
	char *extname = "result";

	int tfields = 2;

	ttype = malloc( tfields * sizeof(char *) );
	tform = malloc( tfields * sizeof(char *) );
	tunit = malloc( tfields * sizeof(char *) );
	
	ttype[0] = "SOURCE";
	ttype[1] = "MASK";

	tform[0] = "E";
	tform[1] = "E";

	tunit[0] = "omitted";
	tunit[1] = "omitted";

	char *pixname = "PIXTYPE";
	char *orderingname = "ORDERING";
	char *pixtype = "HEALPIX";
	char *ordering = "NESTED";
	char *comment = NULL;
	char *nsidename = "NSIDE";
	int *nside;
	int *firstpix;
	int *lastpix;
	int fpix = 0;
	int lpix = naxis2-1;
	int nsi = b->nrow;
	
	printf("\n Lastpix  is  %d",lpix);
	
	nside = &nsi;
	firstpix = &fpix;
	lastpix = &lpix;
	
	int k ;
	for( k=0; k<n_file; k++ )
	{
		if( fits_create_file( &fptr, file[k], &status) )
		{
			printf("\n results_create_fits_results_files:: The result files specified already exist, these will be overwritten");
			//continue;
		}
		
		if( fits_create_tbl( fptr, BINARY_TBL, naxis2, tfields, ttype, tform, tunit, extname, &status ) )
		{
			//printf("\n results_create_fits_results_files:: Binary table creation error");
			//exit(-1);
		}
			
		fits_write_key( fptr, TSTRING, pixname, pixtype, comment, &status );
		fits_write_key( fptr, TSTRING, orderingname, ordering, comment, &status );
		fits_write_key( fptr, TINT, nsidename, nside, comment, &status );
		fits_write_key( fptr, TINT, "FIRSTPIX", firstpix, comment, &status );
		fits_write_key( fptr, TINT, "LASTPIX", lastpix, comment, &status );
		
		fits_close_file( fptr, &status );
	}
	
	
	free( ttype );
	free( tform );
	free( tunit );
}

void results_write_block_to_results(int block_id, struct results *r, struct patch **p, struct block *b, struct hyperpar *h, int simulated)
{
	int k, m, m_ = 0, nout_map = p[ block_id ]->nout_map, n_pixel = b->nrow_subblock * b->ncol_subblock ;
	
	if( simulated ) nout_map +=  1; // change flag to 1 for simulation
	
	for( m=0; m<nout_map; m++ )
	{	
		if( !p[block_id]->template_source[m] )
		{
	
			for( k=0; k<b->n_unmasked[ block_id ]; k++ )
			{
				//if( simulated ) //change flag to 1 to simulate data
				//	r->mu[ m ][ b->idx[ block_id ][ k ] ] = ( (double *) p[ block_id ]->y->x )[ m_ * b->n_unmasked[ block_id ] + k ] ;
				//else
				r->mu[ m ][ b->idx[ block_id ][ k ] ] = ( (double *) p[ block_id ]->mu->x )[ m_ * b->n_unmasked[ block_id ] + k ] ;
				
				//printf("\n values= %.5f", r->mu[ m ][ b->idx[ block_id ][ k ] ] );
			
				if( r->covariance ) r->sdev[ m ][ b->idx[ block_id ][ k ] ] = ( (double *) p[ block_id ]->sdev->x )[ m_ * b->n_unmasked[ block_id ] + k ] ;
			
				r->precision[ m ][ b->idx[ block_id ][ k ] ] = update_precision_link( h->theta[ m + 2] ) ;
			}
		
			m_++;
		}
	}	
	
	if( r->residual )
	{
	
	for( m=0; m < r->nin_map ; m++ )
	{
		for( k=0; k<b->n_unmasked[ block_id ]; k++ )
		{
			r->resid[ m ][ b->idx[ block_id ][ k ] ] = ( (double *) p[ block_id ]->resid->x )[ m * b->n_unmasked[ block_id ] + k ] ;
		}
	}
	}
	
	for( k=0; k<b->n_unmasked[ block_id ]; k++ )
	{
		r->sync_ind[ b->idx[ block_id ][ k ] ] = update_sync_link( h->theta[0] );
		r->dust_ind[ b->idx[ block_id ][ k ] ] = update_gdust_link( h->theta[1] );
	}
	
	return;	
}

void results_write_block_to_input_vector( int k_, int block_id, int residual, double *x, double *resid, struct patch **p, struct block *b, struct hyperpar *h ) 
{
	int k, m, m_ = 0, nin_map = p[ block_id ]->nin_map, nout_map = p[ block_id ]->nout_map, n_pixel = b->nrow_subblock * b->ncol_subblock ;
	
	for( m=0; m<nout_map; m++ )
	{	

			for( k=0; k<b->n_unmasked[ block_id ]; k++ )
			{
				x[ nout_map * n_pixel * k_ + m * n_pixel + b->idx[ block_id ][ k ] ] = ( (double *) p[ block_id ]->mu->x )[ m_ * b->n_unmasked[ block_id ] + k ] ;
				//r->mu[ m ][ b->idx[ block_id ][ k ] ] = ( (double *) p[ block_id ]->mu->x )[ m_ * b->n_unmasked[ block_id ] + k ] ;
				//if( r->covariance ) r->sdev[ m ][ b->idx[ block_id ][ k ] ] = ( (double *) p[ block_id ]->sdev->x )[ m_ * b->n_unmasked[ block_id ] + k ] ;
			}
		
			m_++;
	}	

	if( residual )
	{
		for( m=0; m < nin_map ; m++ )
		{
			for( k=0; k<b->n_unmasked[ block_id ]; k++ )
			{
				resid[ nin_map * n_pixel * k_ + m * n_pixel + b->idx[ block_id ][ k ] ] =  ( (double *) p[ block_id ]->resid->x )[ m * b->n_unmasked[ block_id ] + k ] ;
				//r->resid[ m ][ b->idx[ block_id ][ k ] ] = ( (double *) p[ block_id ]->resid->x )[ m * b->n_unmasked[ block_id ] + k ] ;
			}
		}
	}	

	return;

}

int results_write_patch_result_to_file( int patch_id, struct results *r, struct block *b, char **file )
{
	
	fitsfile *fptr;
	int a, k, l, status = 0, anynul, hdunum=2, hdutype, count,
	    frow, felem=1, nelem = b->n_pixel, longnull=0;
	double nullval, *x, *y;
	
	x = calloc( nelem, sizeof(double) );
	
	if( r->covariance ) y = calloc( nelem, sizeof(double) );
	
	//first row to be read for patch
	frow = patch_id * nelem + 1 ;
	
	//first pass through a does the field and the second pass does the precision
	for( a=0; a<2; a++ )
	{
		for( k=0; k < r->nout_map ; k++ )
		{
			for( l=0; l<nelem; l++ )
			{
				if( !a )
				{
					x[ b->healpix_2d_lattice_to_idx_nested[l]  ] = r->mu[ k ][ l ];
					//Engineering style fix-- revisit later
					//if( k > 0 && r->mu[ k ][ l ] < 0 ) x[ b->healpix_2d_lattice_to_idx_nested[l]  ] = 0.;
					//printf("\n value %d = %.3f", l, r->mu[k][l] );
					if( r->covariance ) y[ b->healpix_2d_lattice_to_idx_nested[l]  ] = r->sdev[ k ][ l ];
				}
				else
				{
					x[ b->healpix_2d_lattice_to_idx_nested[l]  ] = r->precision[ k ][ l ];
				}
			}
			
			if( fits_open_file( &fptr, file[ k + a * 10 ] , READWRITE, &status ) ) return( TRUE ) ;
		
			if( fits_movabs_hdu( fptr, hdunum, &hdutype, &status ) ) return( TRUE );
		
			fits_write_col( fptr, TDOUBLE, 1, frow, felem, nelem, (void  *)x, &status ) ;
	
			fits_close_file( fptr, &status );
			
			if( r->covariance && !a )
			{
				if( fits_open_file( &fptr, file[ k + 4 ] , READWRITE, &status ) ) return( TRUE ) ;
				
				if( fits_movabs_hdu( fptr, hdunum, &hdutype, &status ) ) return( TRUE );
			
				fits_write_col( fptr, TDOUBLE, 1, frow, felem, nelem, (void  *)y, &status ) ;
				
				fits_close_file( fptr, &status );
			}
	
		}
		
		for( l=0; l<nelem; l++)
		{
			if( !a )
				x[ b->healpix_2d_lattice_to_idx_nested[l]  ] = r->sync_ind[ l ];
			else
				x[ b->healpix_2d_lattice_to_idx_nested[l]  ] = r->dust_ind[ l ];
		}
		
			
		if( fits_open_file( &fptr, file[ 8 + a ] , READWRITE, &status ) ) return(TRUE);
		
		if( fits_movabs_hdu( fptr, hdunum, &hdutype, &status ) ) return(TRUE);
		
		fits_write_col( fptr, TDOUBLE, 1, frow, felem, nelem, (void *)x, &status ) ;
	
		fits_close_file( fptr, &status );
		
	}
	
	free( x );
	
	if( r->covariance ) free( y );
	
	//part for the residuals
	

	if( r->residual )
	{
		
		x = calloc( nelem, sizeof( double ) );
		
		for( k=0; k<r->nin_map; k++ )
		{
		
			for( l=0; l<nelem; l++ ) x[ b->healpix_2d_lattice_to_idx_nested[l]  ] = r->resid[ k ][ l ];
				
			if( fits_open_file( &fptr, file[ k + 14 ], READWRITE, &status ) ) return( TRUE );
			
			if( fits_movabs_hdu( fptr, hdunum, &hdutype, &status ) ) return( TRUE );
				
			fits_write_col( fptr, TDOUBLE, 1, frow, felem, nelem, (void  *)x, &status ) ;
		
			fits_close_file( fptr, &status );
		
		}
		
		free( x );
			
	}
	
	return( FALSE  ) ;
	
}


int results_write_patch_result_to_file_0( int patch_id, struct results *r, struct block *b, char **file )
{
	
	FILE *fptr;
	int a, k, l, t, count, felem=1, frow, nelem = b->n_pixel, longnull=0;
	double nullval, *x, *y;
	
	x = calloc( nelem, sizeof(double) );
	
	if( r->covariance ) y = calloc( nelem, sizeof(double) );
	
	//first row to be read for patch
	frow = patch_id * nelem + 1 ;
	
	//first pass through a does the field and the second pass does the precision
	for( a=0; a<2; a++ )
	{
		for( k=0; k < r->nout_map ; k++ )
		{
			for( l=0; l<nelem; l++ )
			{
				if( !a )
				{
					x[ l  ] = r->mu[ k ][ l ];
					if( r->covariance ) y[ l ] = r->sdev[ k ][ l ];
				}
				else
				{
					x[ l ] = r->precision[ k ][ l ];
				}
			}
			
			fptr = fopen( file[ k + a*10 ], "w" ) ;
		
			for( l=0; l<nelem; l++ ) fprintf( fptr, "%lf\n", x[l] ) ;
		
			fclose( fptr );
			
			if( r->covariance && !a )
			{
				
				fptr = fopen( file[ k + 4 ], "w" ) ;
				
				for( l=0; l<nelem; l++ ) fprintf( fptr, "%lf\n", y[l] );
			
				fclose( fptr );
			}
	
		}
		
		for( l=0; l<nelem; l++)
		{
			if( !a )
				x[ l ] = r->sync_ind[ l ];
			else
				x[ l ] = r->dust_ind[ l ];
		}
		
			
		fptr = fopen( file[ 8 + a], "w" );	
			
		for( l=0; l<nelem; l++ ) fprintf( fptr, "%lf\n", x[l] );
		
		fclose( fptr );
		
	}
	
	free( x );
	
	if( r->covariance ) free( y );
	
	//part for the residuals
	

	if( r->residual )
	{
		
		x = calloc( nelem, sizeof( double ) );
		
		for( k=0; k<r->nin_map; k++ )
		{
		
			for( l=0; l<nelem; l++ ) x[ l ] = r->resid[ k ][ l ];
				
			fptr = fopen( file[ k + 14 ], "w" );
			
			for( l=0; l<nelem; l++ ) fprintf( fptr, "%lf\n", x[l] ) ;
			
			fclose( fptr );
		
		}
		
		free( x );
			
	}
	
	return( FALSE  ) ;
	
}



int results_write_patch_result_to_file_2( int patch_id, struct results *r, struct block *b, char **file )
{
	
	fitsfile *fptr;
	int a, k, l, status = 0, anynul, hdunum=2, hdutype, count,
	    frow, felem=1, nelem = b->n_pixel, longnull=0;
	double nullval, *x, *y;
	
	x = calloc( nelem, sizeof(double) );
	
	//first row to be read for patch
	frow = patch_id * nelem + 1 ;
	
	//first pass through a does the field and the second pass does the precision

	for( k=0; k < r->nout_map ; k++ )
	{
		for( l=0; l<nelem; l++ ) x[ b->healpix_2d_lattice_to_idx_nested[l]  ] = r->mu[ k ][ l ];
			
		if( fits_open_file( &fptr, file[ k ] , READWRITE, &status ) ) return( TRUE ) ;
		
		if( fits_movabs_hdu( fptr, hdunum, &hdutype, &status ) ) return( TRUE );
		
		fits_write_col( fptr, TDOUBLE, 1, frow, felem, nelem, (void  *)x, &status ) ;
	
		fits_close_file( fptr, &status );
	
	}	
	
	free( x );

	return( FALSE  ) ;
	
}



void results_simulated_data_copy_hitrate( int npix, int nin_map, char **infile, char **outfile )
{
	//different types of input (in_type) 0: fits 1: ascii txt 2: provided by args obs_intensity and hit_rate
	
	fitsfile *fp;
	
	int 	k, l, status,
			hdu_num = 2,
			hdu_type,
			any_null,
			in_type = 0;
			
	LONGLONG	n = (LONGLONG)npix,
				elem = 1,
				row = 1; 
	
	FILE *fp_a, *fp_a_hit;
			
	double *y = calloc( npix, sizeof(double) ),
			 //*z = calloc( b->n_pixel , sizeof(double) ), //z plays the role of hitrate OR pixel noise variance
			 nulval = FITS_NULL_VAL ;
	
	for( k=0; k<nin_map; k++ )
	{
		status = 0;
		
		printf("\n reading map %d", k );
		
		switch( in_type )
		{
			
			case 0:
				
				if( fits_open_file( &fp, infile[k], READONLY, &status ) ) printf( "\n fits error: unable to open file " ) ;
				if( fits_movabs_hdu( fp, hdu_num, &hdu_type, &status ) ) printf( "\n fits error: unable to move within fits file ") ;
				fits_read_col( fp, TDOUBLE, 2, row, elem, n, &nulval, y, &any_null, &status) ;
				fits_close_file( fp, &status );
				
				if( fits_open_file( &fp, outfile[k], READWRITE, &status ) ) printf( "\n fits error: unable to open file " ) ;
				if( fits_movabs_hdu( fp, hdu_num, &hdu_type, &status ) ) printf( "\n fits error: unable to move within fits file ") ;
				fits_write_col( fp, TDOUBLE, 2, row, elem, n, (void *)y, &status ) ;
				fits_close_file( fp, &status );
				
			break;
			
		}
		
	}
	
	
	free( y );
	
	return;
}







