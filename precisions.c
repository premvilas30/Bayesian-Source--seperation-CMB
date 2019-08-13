//precisions.c

#include "precisions.h"

void precision_prior_gamma_get_rates_shapes( double *rates, double *shapes, int nout_map )
{
	//the precisions should have quite different rates
	
	//decreasing alpha puts more weight on larger precisions (hence smaller sd's) and smoother maps
	
	int k;
	double c = 1., d = 100., *sig0sq = calloc(  6, sizeof(double) ) ; 
	
	for( k=2; k<6; k++ ) shapes[k] = d * d ;
	
	/*sig0sq[2] =  RANGE_CMB * RANGE_CMB  ;
	sig0sq[3] =  RANGE_SYNC * RANGE_SYNC ;
	sig0sq[4] =  RANGE_DUST * RANGE_DUST ;
	sig0sq[5] =  RANGE_FREE_FREE * RANGE_FREE_FREE ;*/
	
	sig0sq[2] =  1./14049.97 ; //.001;//RANGE_CMB * RANGE_CMB  ;
	sig0sq[3] =  1./10000.00;//.001;//RANGE_SYNC * RANGE_SYNC ;
	sig0sq[4] =    1./10000.00;//.001;//RANGE_DUST * RANGE_DUST ;
	sig0sq[5] =  1./10000.00; //.001;//RANGE_FREE_FREE * RANGE_FREE_FREE ;
	
	//for( k=2; k<6; k++ ) sig0sq[k] /= 16. ;
	
	for( k=2; k<6; k++ ) rates[k]  =  c * d * d * sig0sq[k] ;
	
	free( sig0sq );
	
	return;
}

void precision_prior_gamma_get_rates_shapes_1( double *rates, double *shapes, int nout_map )
{
	//the precisions should have quite different rates
	
	//decreasing alpha puts more weight on larger precisions (hence smaller sd's) and smoother maps
	
	int k;
	double c = 1., d = 1., *sig0sq = calloc(  6, sizeof(double) ) ; 
	
	for( k=2; k<6; k++ ) shapes[k] = d * d ;
	
	for( k=2; k<6; k++ ) sig0sq[k] = 0.0001 ;
	
	//sig0sq[2] =  0.0000715 ; //.00845 ; 
	
	for( k=2; k<6; k++ ) rates[k]  =  c * d * d * sig0sq[k] ;
	
	free( sig0sq );
	
	return;
}



void precision_prior_PC_get_lambda( double *lambda )
{
	double lnalpha = log( .05 ), f = -lnalpha , c = 1., *z;
	int k;
	
	z = calloc( 6, sizeof(double) );
	
	/*z[2] = RANGE_CMB / 4. ;
	z[3] = RANGE_SYNC / 4. ;
	z[4] = RANGE_DUST / 4. ;
	z[5] = RANGE_FREE_FREE / 4. ;*/
	
	
	
	for( k=2; k<6; k++ ) 
	{
		z[k] = .001;
		lambda[k] = -lnalpha * c / z[k] ; 
	}
	
	//int k;
	//for( k=2; k<6; k++ ) printf( "\n lambda[%d] = %lf ",k , lambda[k]);
	
	free( z ) ;
	
	return;
}

void precision_initialize_log_precisions( double *theta )
{
	double a = log( 10. );
	
	theta[2] = 2. * ( a - gsl_sf_log( RANGE_CMB / 4. ) );//6.14    12.21    16.06    12.75 
	theta[3] = 2. * ( a - gsl_sf_log(  RANGE_SYNC / 4.) ) ;
	theta[4] = 2. * ( a - gsl_sf_log(  RANGE_DUST / 4.) );
	theta[5] = 2. * ( a - gsl_sf_log(  RANGE_FREE_FREE / 4. ) );
	
	int k;
	//for( k=2; k<6; k++) theta[k] = log( 1000. ) ;
	
	return;
	
}
