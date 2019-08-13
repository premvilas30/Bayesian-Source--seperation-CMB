//blocks.h

#ifndef __BLOCK_H__
#define __BLOCK_H__

#include "required_libs.h"
#include "graphs.h"

struct block
{
	int n_pixel;
	int n_r;
	int n_c;
	int nrow;
	int ncol;
	int nrow_subblock;
	int ncol_subblock;
	int n_block;
	int *n_unmasked;
	double *x;
	int **idx;
	int **mask;
	int *reverse_idx;
	int *healpix_idx_to_2d_lattice_nested;
	int *healpix_2d_lattice_to_idx_nested;
	int *healpix_idx_to_2d_lattice_ring;
	int *healpix_2d_lattice_to_idx_ring;
	
};

struct block *block_create( int nrow, int col, int nrow_subblock, int ncol_subblock, int nin_map, int nin_template, int nout_map, int *mask, int masked_pix_include );

void block_destroy( struct block *b);

void block_compute_healpix_mappings( struct block *b );

void block_print( struct block *b , char *file);

#endif
