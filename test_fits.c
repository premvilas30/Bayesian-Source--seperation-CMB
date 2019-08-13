#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include "fitsio.h"
#include "planck.h"

#define FITS_NULL_VAL 3800000.00

int main()
{

	//try to open the fits file...
	
	fitsfile *fp;
	int status = 0, hdu_num = 2, hdu_type, n_side = 64, n_pixel = n_side * n_side, any_null, id = 0 ;
	double nulval = FITS_NULL_VAL; 
	double *nsamp = calloc( n_pixel, sizeof(double) );
	double *var = calloc( n_pixel, sizeof(double) );
	
	LONGLONG	n = (LONGLONG) n_pixel,
				elem = 1,
				row = id * n + 1 ;
	
	int  k, i;
	
	char **file_name = planck_get_data_file_names();
	
	file_name[0] = "~/astro/Simulation/wmap_simulated_K.fits";
	
	double **storevals  = calloc( 7, sizeof(double*));
	for( k=0; k<7; k++ ) storevals[k] = calloc( n_pixel, sizeof(double));
	
	for( k=0; k<n_pixel; k++ ) storevals[0][k] = 1000.;
	for( k=0; k<n_pixel; k++ ) storevals[1][k] = 2000.;
	
	for( k=0; k<1 ; k++ )
	{
	
		if( fits_open_file( &fp, file_name[k], READWRITE, &status ) ) printf( "\n Mask file: fits error: unable to open file " ) ;

		if( fits_movabs_hdu( fp, hdu_num, &hdu_type, &status ) ) printf( "\n Mask file: fits error: unable to move within fits file ") ;

		fits_write_col( fp, TDOUBLE, 1, row, elem, n, storevals[0] , &status) ;
		printf("\n value of status %d", status );
		fits_write_col( fp, TDOUBLE, 2, row, elem, n, storevals[1] , &status) ;
		printf("\n value of status %d", status );
		//fits_read_col( fp, TDOUBLE, 3, row, elem, n, &nulval, var, &any_null, &status ) ;
	
		fits_close_file( fp, &status );
		
		//for( i=0; i<n_pixel; i++ ) storevals[k][i] = nsamp[i] * var[i] ;
		
	}
	
	FILE *fileout = fopen( "../TEMPLATES/Outvar.txt", "w" );
	
	
	for( k=0; k<1; k++ )
	//	for( i=0; i<n_pixel; i++ ) printf( "%lf\n", storevals[k][i] ); 
	
	fclose( fileout );
	
	for( k=0; k<7; k++ ) free( storevals[k] );
	free(storevals);
	free( nsamp );
	free( var );
	
	return(1);
	
}
