#include "cholmod_sdmult_sub.h"

int cholmod_sdmult_sub( cholmod_sparse *A, double *a, double *b, cholmod_dense *x, cholmod_dense *y)
{
	int j, k, r;
	
	double *yv = (double *) y->x;
	double *xv = (double *) x->x;
	
	int *p = (int *) A->p ;
	int *i = (int *) A->i ;
	double *Av =  (double*) A->x;
	
	for( k=0; k<y->nrow; k++ ) yv[k] = b[0] * yv[k] ;
	
	for( j=0; j<A->ncol; j++ )
	{
		for( k=p[j]; k<p[j+1]; k++ )	
		{
			r = i[ k ] ;
			
			yv[ r ] += a[0] * Av[ k ] *  xv[ j ] ; 
		}
		
	}
	

	return(TRUE);
}
