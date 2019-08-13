
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include "chealpix.h"
#include "healpix.h"

int main()
{
	
	int nside = 512;
	
	int *idx = calloc( nside * nside , sizeof( int ) );
	
	healpix_convert_nested_to_2dlattice_ordering( nside, idx );
	
	FILE *fp = fopen("../TEST_data/healpix_nside_512.txt", "w");
	
	int k, l ;
	
	printf("\n value %d,  nested %d", nside , xyf2nest2(nside, nside-1, nside-1, 0) );
	
	for( k=0; k<nside*nside; k++ ) fprintf( fp, "%d\n", idx[k] );
	
	fclose( fp );
	
	return( 1 );
	
}
