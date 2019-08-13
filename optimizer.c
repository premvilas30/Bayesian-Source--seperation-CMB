//optimizer.c : optimize the log marginal of the hyperparameters

#include "optimizer.h"

struct optpar *optimizer_optpar_create( int nhyper )
{
	
	struct optpar *opt = (struct optpar *)malloc( sizeof( struct optpar ) ) ;
	opt->npar = nhyper;
	opt->opt_theta = calloc( nhyper, sizeof(double) );
	opt->f_opt_theta = -DBL_MAX ;
	opt->grad_opt_theta = calloc( nhyper, sizeof(double) );
	opt->hessian = gsl_matrix_alloc( nhyper, nhyper ) ;
	
	return( opt );
}

void optimizer_optpar_destroy( struct optpar *opt )
{
	int k;
	gsl_matrix_free( opt->hessian );
	free( opt->opt_theta )  ;
	free( opt->grad_opt_theta ) ;
	free( opt ) ;
}


struct optpar *optimizer( struct super *s, int max_iter, int thread_num, int restarts  )
{
	//set up the gsl optimizer
	
	//max_iter=1; 
	
	size_t iter = 0;
	int status, k, id, c, start;
	if( s->individual ) id = s->individual_id[ thread_num ] ; else id = 0;
	int nhyper = s->p[id]->n_parinfer; //s->h[ id ]->nout_map + 2; // + 3;
	//printf("\n Number of inferred parameters is %d", nhyper);
	double dx, size;
	double epsf = 1E-6;//pow(0.005,1.5);
	
	int verbose =  FALSE;
	
	struct optpar *opt = optimizer_optpar_create( nhyper ) ;
	
	struct model_gsl *par = ( struct model_gsl *)malloc( sizeof( struct model_gsl ) );
	par->s = s;
	par->thread_num = thread_num;

	/*allocate vector to hold iterates and initialize from cmbmod*/
	
	gsl_vector *x;
	x = gsl_vector_alloc( nhyper );
	
	if( verbose )
	{
		printf("\n Values in the h[id] struct  : \n ");
		for(k=0; k<nhyper; k++ ) printf( "\t %.2f ", s->h[ id ]->theta[ k ] );
	}
	
	c = 0;
	for( k=0 ; k<s->h[id]->npar; k++ ) 
	{
		if( s->p[id]->parinfer[k] ) { gsl_vector_set( x, c, s->h[ id ]->theta[ k ] ); c++; }
	}
	
	double *startval = calloc( nhyper, sizeof(double) );
	
	for( k=0; k<nhyper; k++ ) startval[k] = gsl_vector_get( x, k );
	
	double *outval = calloc( nhyper, sizeof(double) );
	double Fmin, alpha = 1.0, beta = 0.5, gamma = 2.0,
	        abstol = -DBL_MAX, eps = DBL_EPSILON, intol = sqrt(eps);
	
	int trace = 0, fail=0, fncount, maxit=5000;
	
	for( start=0; start<restarts; start++ )
	{
	
		optimizer_nmmin( nhyper, startval, outval, &Fmin, &fail, abstol, intol,
				alpha, beta, gamma, trace, &fncount, maxit, par  )	;
				
		if( fail ){ printf("\n Optimisation failed....\n"); break; }
		
		for( k=0; k<nhyper; k++ ) startval[k] = outval[k];
		
		if( restarts > 1 ){ fail=0; fncount=0; }
				
	}
	
	
	
	for( k=0; k<opt->npar; k++ ) 
	{
		opt->opt_theta[ k ] = outval[k] ;
	}
	
	opt->f_opt_theta = Fmin;
	
	/*void nmmin(int n, double *Bvec, double *X, double *Fmin,
	   int *fail, double abstol, double intol, void *ex,
	   double alpha, double bet, double gamm, int trace,
	   int *fncount, int maxit, void *params )*/
	
	opt->hessian = optimizer_get_hessian( (void *)par, opt, -1. ) ;
	
	
	/*
	//get the step size
	double step;
	for( k=0; k<nhyper; k++ )
	{
		if( 0.1 * fabs( gsl_vector_get( x, k ) ) > step ) step = 0.1 * fabs( gsl_vector_get( x, k) );
	}
	if( step == 0. ) step = 0.1;
	
	//initial stepsizes for the Nelder-Mead algorithm
	ss = gsl_vector_alloc( nhyper );
	gsl_vector_set_all( ss, step );

	T = gsl_multimin_fminimizer_nmsimplex;

	sf = gsl_multimin_fminimizer_alloc( T, nhyper );
	
	gsl_multimin_fminimizer_set( sf, &G, x, ss ); 

	gsl_vector *x_prev;
	double f_prev = gsl_multimin_fminimizer_minimum( sf );
		
	xx = gsl_multimin_fminimizer_x( sf );
	
	
	x_prev = gsl_vector_alloc( xx->size );
	gsl_vector_memcpy( x_prev, xx );
	
	int status_size = GSL_CONTINUE;
	int status_x = GSL_CONTINUE;
	int status_f = GSL_CONTINUE;
	
	int restarted = 0;
	
	size = gsl_multimin_fminimizer_size( sf );
	
	printf("\n=== Size at beginning of optimisation... %.20f ===\n",  size);

	do{
		
		//if( iter == max_iter-1 && restarted < 1)//restarts ) 
		//{
			//gsl_multimin_fminimizer_restart( sf ) ; 
			//restarted += 1 ; 
		//}
		
		if( verbose )
		{
			printf("\n ------ parameter values --------");
			for( k=0; k<opt->npar; k++ ) printf(" %.20f   ", gsl_vector_get( xx , k) ) ;
		}
	
		iter++;
		status = gsl_multimin_fminimizer_iterate( sf );
		
		if( status ) break;
		
		size = gsl_multimin_fminimizer_size( sf );
		status_size = gsl_multimin_test_size( size, .001);
		
		printf("\n=== Size at beginning of optimisation... %.20f ===\n",  size);
		
		xx = gsl_multimin_fminimizer_x( sf );
		
		if(x_prev)
		{
			
			size_t i_s;
			
			for( i_s=0, dx=0.; i_s<xx->size; i_s++ )
			{
				dx += gsl_pow_2( gsl_vector_get( xx, i_s ) - gsl_vector_get( x_prev, i_s ) );
			}
			dx = sqrt( dx / xx->size );
			status_x = gsl_multimin_test_size( dx, 0.0001 );
		}
		else
		{
			status_x = GSL_CONTINUE;
		}

		if(!x_prev){
			x_prev = gsl_vector_alloc( xx->size );
		}
		
		gsl_vector_memcpy( x_prev, xx );
	
		double df = fabs( f_prev - gsl_multimin_fminimizer_minimum( sf ) );
		status_f = gsl_multimin_test_size( df, epsf );
		f_prev = gsl_multimin_fminimizer_minimum( sf );

	}while( (status_size == GSL_CONTINUE) && iter < max_iter ); 
	
	
	if( verbose )
		{
			printf("\n ------ parameter values --------");
			for( k=0; k<opt->npar; k++ ) printf(" %.2f   ", gsl_vector_get( xx , k) ) ;
		}
	
	xx = gsl_multimin_fminimizer_x( sf );
	
	c = 0;
	for( k=0; k<opt->npar; k++ ) 
	{
		opt->opt_theta[ k ] = gsl_vector_get( xx, k ) ;
	}
	
	opt->f_opt_theta = gsl_multimin_fminimizer_minimum( sf ) ;
	
	//compute the Hessian at the mode
	
	opt->hessian = optimizer_get_hessian( (void *)par, opt, -1. ) ;
	
	//free vectors used for minimization
	
	gsl_multimin_fminimizer_free( sf );
	
	gsl_vector_free( ss );
	
	if(x_prev){
		gsl_vector_free(x_prev);
	}*/
	
	free( par );
	gsl_vector_free( x );
	free(startval);
	free(outval);
	
	return( opt );
}




struct optpar *optimizer_zzz( struct super *s, int max_iter, int thread_num, int restarts  )
{
	//set up the gsl optimizer
	
	//max_iter=1; 
	
	size_t iter = 0;
	int status, k, id, c;
	if( s->individual ) id = s->individual_id[ thread_num ] ; else id = 0;
	int nhyper = s->p[id]->n_parinfer; //s->h[ id ]->nout_map + 2; // + 3;
	//printf("\n Number of inferred parameters is %d", nhyper);
	double dx;
	double epsf = 1E-6;//pow(0.005,1.5);
	
	int verbose =  TRUE;
	
	struct optpar *opt = optimizer_optpar_create( nhyper ) ;
	
	struct model_gsl *par = ( struct model_gsl *)malloc( sizeof( struct model_gsl ) );
	par->s = s;
	par->thread_num = thread_num;
	
	//printf("\nallocate gsl_multimin");
	
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *sf;
	
	gsl_vector *x;
	gsl_multimin_function_fdf G;
	
	gsl_vector *xx;
	
	G.n = nhyper;
	G.f = gsl_objective_function ;
	G.df = gsl_objective_function_gradient ;
	G.fdf = gsl_objective_function_and_gradient  ; //to write this
	G.params = (void *)par;
	
	/*allocate vector to hold iterates and initialize from cmbmod*/
	x = gsl_vector_alloc( nhyper );
	
	//printf("\n Values in the h[id] struct  : \n ");
	//for(k=0; k<nhyper; k++ ) printf( "\t %.2f ", s->h[ id ]->theta[ k ] );
	
	c = 0;
	for( k=0 ; k<s->h[id]->npar; k++ ) 
	{
		if( s->p[id]->parinfer[k] ) { gsl_vector_set( x, c, s->h[ id ]->theta[ k ] ); c++; }
	}
	
	//printf("\n Before call to gsl_sf :\n");
	//for( k=0; k<nhyper; k++ ) printf("\t %.2f  ", gsl_vector_get( x , k) ) ;
	
		
	//T = gsl_multimin_fdfminimizer_vector_bfgs2;
	T = gsl_multimin_fdfminimizer_vector_bfgs2;

	sf = gsl_multimin_fdfminimizer_alloc( T, nhyper );
	
	gsl_multimin_fdfminimizer_set( sf, &G, x, 1, 0.1 ); 

	gsl_vector *x_prev;
	double f_prev = gsl_multimin_fdfminimizer_minimum( sf );
	//printf("\nset values %lf",-f_prev);
		
	xx = gsl_multimin_fdfminimizer_x( sf );
	
	//printf("\n Values initial gsl :\n");
	//for( k=0 ; k<nhyper; k++ ) printf("\t %.2f " , gsl_vector_get( xx, k ) ); 
	
	x_prev = gsl_vector_alloc( xx->size );
	gsl_vector_memcpy( x_prev, xx );
	
	int status_g = GSL_CONTINUE;
	int status_f = GSL_CONTINUE;
	int status_x = GSL_CONTINUE;
	
	//exit(1);
	
	//printf("\n ------ parameter values --------");
	//for( k=0; k<opt->npar; k++ ) printf(" %.2f  ", gsl_vector_get( sf->gradient , k) ) ;
	
	int restarted = 0;

	do{
		
		if( iter == max_iter-1 && restarted < restarts ) 
		{
			gsl_multimin_fdfminimizer_restart( sf ) ; 
			restarted += 1 ; 
		}
		
		if( verbose )
		{
			printf("\n ------ parameter values --------");
			for( k=0; k<opt->npar; k++ ) printf(" %.2f   ", gsl_vector_get( xx , k) ) ;
			printf("\n ------ gradient values --------");
			for( k=0; k<opt->npar; k++ ) printf(" %.2f   ", gsl_vector_get( sf->gradient , k) ) ;
		}
	
		iter++;
		status = gsl_multimin_fdfminimizer_iterate( sf );
		
		//printf("\n Value of : \n \t f :: %lf", gsl_multimin_fdfminimizer_minimum( sf ));
			 
		status_g = gsl_multimin_test_gradient( sf->gradient, 0.01 ); /*tolerance on gradient may need to be reset*/
		//double gnrm2 = gsl_blas_dnrm2( sf->gradient );
		
		xx = gsl_multimin_fdfminimizer_x( sf );
		
		if(x_prev)
		{
			
			size_t i_s;
			
			for( i_s=0, dx=0.; i_s<xx->size; i_s++ )
			{
				dx += gsl_pow_2( gsl_vector_get( xx, i_s ) - gsl_vector_get( x_prev, i_s ) );
			}
			dx = sqrt( dx / xx->size );
			status_x = gsl_multimin_test_size( dx, 0.0001 );
		
		}
		else
		{
			status_x = GSL_CONTINUE;
		}

		if(!x_prev){
			x_prev = gsl_vector_alloc( xx->size );
		}
		
		gsl_vector_memcpy( x_prev, xx );
	
		double df = fabs( f_prev - gsl_multimin_fdfminimizer_minimum( sf ) );
		status_f = gsl_multimin_test_size( df, epsf );
		f_prev = gsl_multimin_fdfminimizer_minimum( sf );
		
		//printf("\n Value of : \n \t f :: %lf", gsl_multimin_fdfminimizer_minimum( sf ));

	}while( (status_g == GSL_CONTINUE) && iter < max_iter ); /*reset max its in preamble*/
	
	
	if( verbose )
		{
			printf("\n ------ parameter values --------");
			for( k=0; k<opt->npar; k++ ) printf(" %.2f   ", gsl_vector_get( xx , k) ) ;
			printf("\n ------ gradient values --------");
			for( k=0; k<opt->npar; k++ ) printf(" %.2f   ", gsl_vector_get( sf->gradient , k) ) ;
		}
	
	xx = gsl_multimin_fdfminimizer_x( sf );
	
	//if( s->p[id]->parinfer[k] )
	
	c = 0;
	for( k=0; k<opt->npar; k++ ) 
	{
		opt->opt_theta[ k ] = gsl_vector_get( xx, k ) ;
		opt->grad_opt_theta[ k ] = gsl_vector_get( sf->gradient, k );
	}
	
	opt->f_opt_theta = gsl_multimin_fdfminimizer_minimum( sf ) ;
	
	//compute the Hessian at the mode
	
	opt->hessian = optimizer_get_hessian( (void *)par, opt, -1. ) ;
	
	/*free vectors used for minimization*/
	
	gsl_multimin_fdfminimizer_free( sf );
	gsl_vector_free( x );
	
	if(x_prev){
		gsl_vector_free(x_prev);
	}
	
	free( par );
	
	return( opt );
}

gsl_matrix *optimizer_get_hessian( void *par, struct optpar *opt, double f )
{

//	struct super *s = (struct super *) par->s ;
	
	int id;

	int k, l, c, npar = opt->npar ;
	
	//printf("\n The value of npar is %d", npar ); 
	
	double f00, fm1p0, fp1p0, fp1p1, fp1m1, fm1p1, fm1m1, deriv, absx, dx, absy, dy, eps = DBL_EPSILON, relStep = sqrt(eps), *t;
	
	t = opt->opt_theta ;

	gsl_vector *x = gsl_vector_alloc( (size_t) npar );

	//problem here with params inferred


	for(k = 0; k < opt->npar; k++ )  
	{
		gsl_vector_set( x, k, t[ k ] ) ; 
		//printf("\n opt[%d] = %.20f", k , t[k]);
	}
	

	gsl_matrix *Hess = gsl_matrix_alloc( (size_t) npar, (size_t) npar ) ;

	f00 = gsl_objective_function( x, par ) ;
	
	//compute Hessian entries using finite differences
	
	//printf("\n sqrt(double eps) = %.20f", relStep );

	for( k=0; k<npar ; k++ )
	{
		absx = fabs( t[k] );
		
		if( absx < 2E-5 ) dx = 1E-4; else dx = relStep * fabs( t[k] ) ;
	
		dx = 0.01;
	
		//if( k==0 ) dx = 0.1;
	
		//get the non-mixed derivative
		gsl_vector_set( x, k , t[k] - dx );
		fm1p0 = gsl_objective_function( x, par ) ;
		
		gsl_vector_set( x, k, t[k] + dx ) ;
		fp1p0 = gsl_objective_function( x, par ) ; 
		
		deriv =  fabs( f * ( - fm1p0 + 2. * f00 - fp1p0 ) / ( dx * dx ) ) ;
		
		//printf("\n par : %d \n\t f00 = %.20f \n\t fm1p0 = %.20f \n\t fp1p0 = %.20f \n\t deriv = %.20f", k, f00, fm1p0, fp1p0, deriv );
		
		gsl_matrix_set( Hess, k, k, deriv )  ;
		
		for( l=0; l<k; l++ )
		{
			//mixed derivative 
			
			absy = fabs( t[l] );
			if( absy < 2E-5 ) dy = 1E-4; else dy = relStep * fabs( t[l] ) ;
			
			dx = 0.01; dy = 0.01;
			
			gsl_vector_set( x, k, t[k] + dx );
			gsl_vector_set( x, l, t[l] + dy );
			fp1p1 = gsl_objective_function( x, par );
			
			gsl_vector_set( x, k, t[k] + dx );
			gsl_vector_set( x, l, t[l] - dy );
			fp1m1 = gsl_objective_function( x, par );
			
			gsl_vector_set( x, k, t[k] - dx );
			gsl_vector_set( x, l, t[l] + dy );
			fm1p1 = gsl_objective_function( x, par );
			
			gsl_vector_set( x, k, t[k] - dx );
			gsl_vector_set( x, l, t[l] - dy );
			fm1m1 = gsl_objective_function( x, par );
			
			deriv = f * ( fp1p1 - fp1m1 - fm1p1 + fm1m1 ) / ( 4. * dx * dy ) ;
			
			gsl_matrix_set( Hess, k, l, deriv );
			gsl_matrix_set( Hess, l, k, deriv );
		
		}	
	}
	
	printf("\n **** Hessian matrix **** \n");
	for( k=0; k<npar; k++ )
	{
		for( l=0; l<npar; l++ )
		{
			printf(" %.3f , " , gsl_matrix_get( Hess, k, l ) );
		}
		printf("\n");
	}
	
	return( Hess ) ; 
}

/***************************** trial optimizer for lower no. of parameters ************************************/

struct optpar *optimizer2( struct super *s, int max_iter, int thread_num  )
{
	//set up the gsl optimizer 
	
	size_t iter = 0;
	int status, k, id = s->individual_id[ thread_num ] , nhyper = 4;
	double dx;
	double epsf = 1E-6;//pow(0.005,1.5);
	
	int verbose = FALSE;
	
	struct optpar *opt = optimizer_optpar_create( nhyper ) ;
	
	struct model_gsl *par = ( struct model_gsl *)malloc( sizeof( struct model_gsl ) );
	par->s = s;
	par->thread_num = thread_num;
	
	//printf("\nallocate gsl_multimin");
	
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *sf;
	
	gsl_vector *x;
	gsl_multimin_function_fdf G;
	
	gsl_vector *xx;
	
	G.n = nhyper;
	G.f = gsl_objective_function2 ;
	G.df = gsl_objective_function_gradient2 ;
	G.fdf = gsl_objective_function_and_gradient2  ; //to write this
	G.params = (void *)par;
	
	/*allocate vector to hold iterates and initialize from cmbmod*/
	x = gsl_vector_alloc( nhyper );
	
	for( k=0 ; k<nhyper; k++ ) gsl_vector_set( x, k, s->h[ id ]->theta[ k+2 ] );
	
	//for( k=0; k<opt->npar; k++ ) printf(" %.2f  ", gsl_vector_get( x , k) ) ;
	
		
	//T = gsl_multimin_fdfminimizer_vector_bfgs2;
	T = gsl_multimin_fdfminimizer_vector_bfgs2;

	sf = gsl_multimin_fdfminimizer_alloc( T, nhyper );
	
	gsl_multimin_fdfminimizer_set( sf, &G, x, 1, 0.1 ); 

	gsl_vector *x_prev;
	double f_prev = gsl_multimin_fdfminimizer_minimum( sf );
	//printf("\nset values %lf",-f_prev);
		
	xx = gsl_multimin_fdfminimizer_x( sf );
	x_prev = gsl_vector_alloc( xx->size );
	gsl_vector_memcpy( x_prev, xx );
	
	int status_g = GSL_CONTINUE;
	int status_f = GSL_CONTINUE;
	int status_x = GSL_CONTINUE;
	
	//printf("\n ------ parameter values --------");
	//for( k=0; k<opt->npar; k++ ) printf(" %.2f  ", gsl_vector_get( sf->gradient , k) ) ;

	do{
		
		if( verbose )
		{
			printf("\n ------ parameter values --------");
			for( k=0; k<opt->npar; k++ ) printf(" %.2f   ", gsl_vector_get( xx , k) ) ;
			printf("\n ------ gradient values --------");
			for( k=0; k<opt->npar; k++ ) printf(" %.2f   ", gsl_vector_get( sf->gradient , k) ) ;
		}
	
		iter++;
		status = gsl_multimin_fdfminimizer_iterate( sf );
			 
		status_g = gsl_multimin_test_gradient( sf->gradient, 0.05 ); /*tolerance on gradient may need to be reset*/
		//double gnrm2 = gsl_blas_dnrm2( sf->gradient );
		
		xx = gsl_multimin_fdfminimizer_x( sf );
		
		if(x_prev)
		{
			
			size_t i_s;
			
			for( i_s=0, dx=0.; i_s<xx->size; i_s++ )
			{
				dx += gsl_pow_2( gsl_vector_get( xx, i_s ) - gsl_vector_get( x_prev, i_s ) );
			}
			dx = sqrt( dx / xx->size );
			status_x = gsl_multimin_test_size( dx, 0.0001 );
		
		}
		else
		{
			status_x = GSL_CONTINUE;
		}

		if(!x_prev){
			x_prev = gsl_vector_alloc( xx->size );
		}
		
		gsl_vector_memcpy( x_prev, xx );
	
		double df = fabs( f_prev - gsl_multimin_fdfminimizer_minimum( sf ) );
		status_f = gsl_multimin_test_size( df, epsf );
		f_prev = gsl_multimin_fdfminimizer_minimum( sf );
		
		//printf("\n Value of : \n \t f :: %lf", gsl_multimin_fdfminimizer_minimum( sf ));

	}while( (status_g == GSL_CONTINUE) && ( status_f == GSL_CONTINUE )  && iter < max_iter ); /*reset max its in preamble*/
	
	
	xx = gsl_multimin_fdfminimizer_x( sf );
	
	for( k=0; k<opt->npar; k++ ) 
	{
		opt->opt_theta[ k ] = gsl_vector_get( xx, k ) ;
		opt->grad_opt_theta[ k ] = gsl_vector_get( sf->gradient,  k );
	}
	
	opt->f_opt_theta = gsl_multimin_fdfminimizer_minimum( sf ) ;
	
	//compute the Hessian at the mode
	
	//opt->hessian = optimizer_get_hessian( (void *)par, opt, -1. ) ;
	
	/*free vectors used for minimization*/
	
	gsl_multimin_fdfminimizer_free( sf );
	gsl_vector_free( x );
	
	if(x_prev){
		gsl_vector_free(x_prev);
	}
	
	free( par );
	
	return( opt );
}





