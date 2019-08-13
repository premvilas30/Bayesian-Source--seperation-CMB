//blocks.h

#ifndef __BLOCKS_H__
#define __BLOCKS_H__

#include "required_libs.h"

struct block
{
	int n_r;
	int n_c;
	int n_block;
	double *x;
	int **idx;
	int *reverse_idx;
};

struct block *block_create( int nrow, int col, int nrow_subblock, int ncol_subblock, int *mask );

void destroy_block( struct block *b);

#endif
