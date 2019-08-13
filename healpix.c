/*
	healpix.c
	
	Functions to do operations on Healpix orderings-- to allow for index referencing
	
	Written by: Jason Wyse,
				School of Computer Science and Statistics,
				Trinity College Dublin,
				Dublin 2, 
				Ireland
			
	Last modified: Thu 21 Jul 2016 16:38:26 IST 
	
*/

#include "healpix.h"

void healpix_convert_nested_to_2dlattice_ordering( int nside, int *lattice_pix_number )
{
	//this gives the number of the pixel in the 2d lattice corresponding to pixel numbers 0:(nside^2-1)
	// in the nested healpix ordering for a base resolution patch-- any set of indexes are these values
	// modulo nside^2
	
	int k, kk, npixel = nside * nside, *pix = calloc( npixel, sizeof(int) ), nbit = 32;
	
	//store the index numbers as 32 bit integers in binary
	
	int **binary_form = calloc( npixel, sizeof(int*) ); 
	for( k=0; k<npixel; k++ ) binary_form[k] = calloc( nbit, sizeof(int) );
	
	int *largest_one = calloc( npixel, sizeof(int) );

	int c, n, done;
	
	for( k=0; k<npixel; k++ )
	{
		done = 0;
		kk =  k;
		
		for( c=nbit-1 ; c > -1; c-- )
		{
			n = kk >> c ;
			if( n & 1 ) binary_form[k][c] = 1; 
		}
		
	}
	
	//printf("\n The binary form of 200  is \n");
	//for( c = nbit-1; c>-1; c-- ) printf( " %d, ", binary_form[200][c]);
	
	//now convert the binary form into x and y coordinates on the lattice
	
	int *xpos = calloc( npixel, sizeof(int) ), *ypos = calloc( npixel, sizeof(int) ), l ;
	
	for( k=0; k<npixel; k++ )
	{
		//compute the x position using the binary representation of the 
		//  even positions
		
		//even positions
		
		l = 0 ;
		
		for( c = 0; c < nbit ; c+=2 ) 
		{
			xpos[k] += binary_form[k][c] * pow( 2 , c );
		}

		//odd positions
		
		for( c = 1; c < nbit ; c+=2 ) 
		{
			ypos[k] += binary_form[k][c] * pow( 2 , c );
		}
		
	}
	
	//now compute the 2d lattice pixel number
	
	for( k=0; k<npixel; k++ ) lattice_pix_number[ k ] = ypos[k] * nside + xpos[k] ;
	
	
	for( k=0; k<npixel; k++ ) free( binary_form[k] );
	free( binary_form );
	free( largest_one );
	
	free( xpos );
	free( ypos );
	
	free( pix );
	
	return;
	
}
