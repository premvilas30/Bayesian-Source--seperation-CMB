//test program for blocks

#include "block.h"

int main()
{
	int n = 4, N = n*n, k; 
	int *mask = calloc( N , sizeof(int) );
	for( k=0; k<N; k++ ) mask[k] = TRUE;
	mask[3] = FALSE;
	
	struct block *b = block_create( n, n, 2, 2, mask );
	
	block_print( b , "block_out.txt" );
	
	block_destroy( b );
	
	free(mask);
	
	return(1);
}
