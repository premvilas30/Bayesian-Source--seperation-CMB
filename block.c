//block.c

#include "block.h"

struct block *block_create( int nrow, int ncol, int nrow_subblock, int ncol_subblock, int nin_map, int nin_template, int nout_map, int *mask, int masked_pix_include )
{
	//construct the back referencing
	
	if( nrow % nrow_subblock != 0 || ncol % ncol_subblock != 0 ) printf("\n Error-- cannot create block");
	
	int 	n_block = nrow * ncol / ( nrow_subblock * ncol_subblock ), 
			n_pixel = nrow * ncol,
			n_pixel_subblock = nrow_subblock * ncol_subblock,
			n_r = nrow/nrow_subblock,
			n_c = ncol/ncol_subblock ; 
	
	struct block *b = (struct block *)malloc( sizeof( struct block ) );
	b->idx = calloc( n_block+1, sizeof(int *) );
	b->reverse_idx = calloc( n_pixel, sizeof(int) );
	b->n_unmasked = calloc( n_block+1, sizeof(int) );
	b->mask = calloc( n_block+1, sizeof(int *) );
	b->n_pixel = nrow * ncol ;
	b->nrow = nrow;
	b->ncol = ncol;
	b->nrow_subblock = nrow_subblock;
	b->ncol_subblock = ncol_subblock;
	b->n_r = n_r;
	b->n_c = n_c;
	b->n_block = n_block;
	b->x = calloc( n_pixel * (nout_map - nin_template ) , sizeof(double) );
	b->healpix_idx_to_2d_lattice_nested = calloc( n_pixel, sizeof(int) );
	b->healpix_2d_lattice_to_idx_nested = calloc( n_pixel, sizeof(int) );
	b->healpix_idx_to_2d_lattice_ring = calloc( n_pixel, sizeof(int) );
	b->healpix_2d_lattice_to_idx_ring = calloc( n_pixel, sizeof(int) );
	
	
	block_compute_healpix_mappings( b );
	
	//remember to remap the mask entries here as they haven't been mapped
	
	int unmasked = 0, k, l, r, c, i, a, pix, pix__;
	
	for( k=0; k<n_pixel; k++ ) 
	{
		l = b->healpix_2d_lattice_to_idx_nested[ k ];
		if( !masked_pix_include ) unmasked += mask[l]; else unmasked += 1; //trick to force mask
	}
	
	b->n_unmasked[0] = unmasked;
	
	//do the main block first
	b->idx[0] = calloc( unmasked, sizeof(int) );
	b->mask[0] = calloc( n_pixel, sizeof(int) );

	FILE *mask_chk = fopen("../WMAP_data/mask_check.txt","w");

	i=0;
	for( k=0; k<n_pixel; k++ )
	{
		l = b->healpix_2d_lattice_to_idx_nested[ k ];
		b->mask[0][k] = mask[l];
		
		fprintf(  mask_chk, "%d\n", mask[l]);
		
		if( !masked_pix_include )
		{
			if( mask[ l ] )
			{ 
				b->idx[0][i] = k; 
				b->reverse_idx[k] = i;
				i++; 
			}
		}
		else
		{
			b->idx[0][i] = k; 
			b->reverse_idx[k] = i;
			i++; 
		}
	
	}
	
	fclose(  mask_chk );
	
	
	//create a structure to store the indexes in each block
	for( k=0; k<n_r; k++ )
	{
		for( l=0; l<n_c; l++ )
		{
			unmasked = 0;
			i = 0;
			a = 0;
			
			b->mask[ k * n_c + l + 1 ] = calloc( n_pixel_subblock, sizeof(int) );
			
			while( a<2 )
			{
				
				if( a ) b->idx[ k * n_c + l + 1 ] = calloc( unmasked, sizeof(int) ) ;
				for( r=0; r<nrow_subblock; r++ )
				{
					for( c=0; c<ncol_subblock; c++ )
					{
						pix = ( k * nrow_subblock + r ) * ncol + l * ncol_subblock + c ;
						pix__ = b->healpix_2d_lattice_to_idx_nested[ pix ];
						if( a ) b->mask[  k * n_c + l + 1 ][ r * ncol_subblock + c ] = mask[ pix__ ];
						
						if( !masked_pix_include )
						{
							if( mask[ pix__ ] )
							{  
								if( a ) { b->idx[ k * n_c + l + 1 ][i] = pix;  i++; } else unmasked++;
							}
						}
						else
						{
							if( a ) { b->idx[ k * n_c + l + 1 ][i] = pix;  i++; } else unmasked++;
						}
					}
				}
				if( !a ) b->n_unmasked[ k * n_c + l + 1 ] = unmasked;
				//printf("\n The value of i is  %d",i);
				a++;
				
			}
			
		}
	}
	
	return( b );
	
}

void block_destroy( struct block *b)
{
	int k;
	for( k=0; k<b->n_block+1 ; k++ ) 
	{
		free( b->idx[k] );
		free( b->mask[k] );
	}
	free( b->idx );
	free( b->mask );
	free( b->reverse_idx );
	free( b->n_unmasked );
	free( b->x );
	free( b->healpix_idx_to_2d_lattice_nested );
	free( b->healpix_2d_lattice_to_idx_nested );
	free( b->healpix_idx_to_2d_lattice_ring );
	free( b->healpix_2d_lattice_to_idx_ring );
	free( b );
	return;
}


void block_compute_healpix_mappings( struct block *b )
{
	//this scans a pre-computed file for each nside case
	
	int n_side = (int) sqrt( b->n_pixel ), x, y;
	FILE *fp;
	char *file;
	
	/*if( n_side == 512)
			{ file = "../HEALPIX_orderings/HEALPix_nested_and_GMRF_native_pixel_indexing_Nside_512.txt"; }
	else if( n_side == 1024)
			{ file = "../HEALPIX_orderings/HEALPix_nested_and_GMRF_native_pixel_indexing_Nside_1024.txt"; }
	else
		 { file = "../HEALPIX_orderings/sequential_ordering.txt" ; }*/
	
	int k;
	
	//Try using the chealpix library to compute the mappings : lattice is numbered
	
	for( k=0; k<b->n_pixel; k++ )
	{
		y = k / n_side ;
		x = k - y * n_side ;
		
		b->healpix_2d_lattice_to_idx_nested[k] = xyf2nest2( n_side, x, y, 0) ;
		b->healpix_idx_to_2d_lattice_nested[ b->healpix_2d_lattice_to_idx_nested[k] ] = k ;
		
		//b->healpix_2d_lattice_to_idx_ring[k] = xyf2ring2( n_side, x, y, 0 );
		//b->healpix_idx_to_2d_lattice_ring[ b->healpix_2d_lattice_to_idx_ring[k] ] = k ;
		
	}
	
	
	/*if( FALSE )
	{
	
	fp = fopen( file, "r" );
	
	for( k=0; k<b->n_pixel; k++ ) fscanf( fp, "%d", &b->healpix_2d_lattice_to_idx[k] ); //gives the healpix position corresponding to the lattice index 
	
	for( k=0; k<b->n_pixel; k++ ) b->healpix_idx_to_2d_lattice[ b->healpix_2d_lattice_to_idx[k] ] = k; //gives the lattice position corresponding to the healpix index
	
	fclose( fp );
	
	}*/
	
	//for(k=0; k<b->n_pixel; k++) printf( "\n idx_to_2d: %d \t 2d_to_idx: %d", b->healpix_idx_to_2d_lattice[k], b->healpix_2d_lattice_to_idx[k] ) ;
	
	return;	
}


void block_print( struct block *b , char *file)
{
	int k, l, n_pixel_subblock = b->n_pixel/( b->n_r * b->n_c);
	FILE *fp;
	
	fp = fopen( file, "w" );
	
	for( k=0; k< 1/*b->n_block+1*/; k++ )
	{
		//fprintf( fp, "\n  *** block %d *** \n index in image \n", k);
		//for( l=0; l<b->n_unmasked[k]; l++) fprintf( fp, "%d\t", b->idx[k][l] );
		//fprintf( fp,  "\n mask  \n ");
		for( l=0; l< (k>0)*n_pixel_subblock + (k==0)*b->n_pixel ; l++) fprintf( fp, "%d\t", b->mask[k][l] );
	}
	//fprintf(fp,"\n reverse map :\n");
	//for( l=0; l<b->n_pixel; l++)  fprintf( fp, "%d\t", b->reverse_idx[l] );
	fclose( fp );
	return;
}






