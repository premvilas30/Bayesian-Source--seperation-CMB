// code to convert a stacked (vector column) fits file to 
// a single column file

//for this the number of fields in the vector column TFIELDS must equal to 1

void fitsroutine_convert_vector_column_to_scalar_column( char *file_name, char *new_file_name, colno, char *new_ordering )
{

	int tfields = 1 ; //only does one column at a time!
	
	fitsfile *fp;
	int 	status = 0, hdu_num = 2, 
			hdu_type, 
			n_pixel = n_side * n_side,
			n_axis2 = 12 * n_pixel, 
			any_null, id = 0 ;
			
	double nulval = FITS_NULL_VAL; 
	
	double *x = calloc( n_pixel, sizeof(double) );
	
	LONGLONG	n = (LONGLONG) n_pixel,
				elem = 1,
				row =  1 ;

	if( fits_open_file( &fp, file_name, READONLY, &status ) ) return(1) ;
	
	//find the nside
	
	//get the essential keywords from the file
	char *pixtype, *ordering, *indxschm,
		*ttype, *tform, *tunit, *comment = NULL, *extname = " ",;
	
	fits_read_keyword( fp, "PIXTYPE", pixtype, comment, &status );
	fits_read_keyword( fp, "ORDERING", ordering, comment, &status );
	fits_read_keyword( fp, "INDXSCHM", indxschm, comment, &status );
	
	fits_read_keyword( fp, "TTYPE1", ttype, comment, &status );
	fits_read_keyword( fp, "TFORM1", tform, comment, &status );
	fits_read_keyword( fp, "TUNIT1", tunit, comment, &status );

	if( fits_movabs_hdu( fp, hdu_num, &hdu_type, &status ) ) ;

	fits_read_col( fp, TDOUBLE, colno, row, elem, n, &nulval, x, &any_null, &status) ;
		
	fits_close_file( fp, &status );
	
	//now write x back out to a single column fits file
	
	if( fits_create_file( &fp, new_file_name, &status ) ) return(1);
	
	if( fits_create_table( fp, BINARY_TBL, n_axis2, tfields, ttype,  tform, tunit, extname, &status ) ) return(1);
	
	fits_write_key( fp, TSTRING, )

}
