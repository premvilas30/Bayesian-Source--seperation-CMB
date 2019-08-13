//test program

#include "patch.h"
#include <math.h>

int main()
{
	int n = 16, N = n*n, nin_map = 1, nout_map = 1, k; 
	int *mask = calloc( N , sizeof(int) );
	for( k=0; k<N; k++ ) mask[k] = TRUE;
	
	struct patch *p = patch_create( n, n, nin_map, nout_map, mask );
	
	patch_destroy( p );
	
	free(mask);
	
	return(1);
}
