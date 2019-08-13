/**** functions for gsl ****/

#include "model_gsl.h"

double gsl_objective_function( const gsl_vector *x, void *params )
{
	struct model_gsl *par = ( struct model_gsl *)params;
	struct super *s = par->s;
	
	int id;
	if( s->individual ) id = s->individual_id[ par->thread_num ]; else id = 0;
	
	struct hyperpar *h = s->h[ id ];
	struct patch **p = s->p;
	struct patch *p0;
	int *mask;
	struct block *b = s->b;
	double *store_theta = calloc( h->npar, sizeof(double) );
	
	//printf("\n gsl_objective_function:  evaluate at  " )
	
	int k, c=0;
	for( k=0; k<h->npar; k++ )
	{
		store_theta[k] = h->theta[k];
		if( s->p[id]->parinfer[k] ) 
		{ 
			h->theta[k] = gsl_vector_get( x, c ); 
			c++; 
		}
	}
	
	
	int nbl; //number of blocks (if analyzing individually)
	if( s->individual ) nbl = 1; else nbl = b->n_block;
	
	if( s->individual )
	{
		k = id ;  
		//printf("\n the value of  k is %d ", id);
		update_patch( p[k], h, TRUE, b->mask[k] ) ;
		p0 = p[k] ;
		mask = b->mask[k];
	}
	else
	{
		update_patch( p[0], h, FALSE, b->mask[0] );
		
		#pragma omp parallel for 
		for( k=1; k< b->n_block + 1 ; k++ )
		{
			if(b->n_unmasked[k] > 0) update_patch( p[k] , h , TRUE, b->mask[k] );
		}
		patch_combine( p , b );
		p0 = p[0] ; 
		mask = b->mask[0];
	}
	
	double log_obj = model_log_joint_marginal_hyperpar( p0, p , h , nbl, p0->masked_pix_include, mask, id ) ;
	
	for( k=0; k<h->npar; k++ ) 
	{
		h->theta[k] = store_theta[k];
	}
	free(store_theta);
	
	//printf("\n");
	//for( k=0; k<h->npar; k++) printf(" %.3f ," , h->theta[k]);
	//printf("\n log_posterior_function: val %lf", log_obj );
	
	return( -log_obj ); //taking minimum here as gsl minimize the function
}

void gsl_objective_function_gradient( const gsl_vector *x, void *params, gsl_vector *grad )
{	

		double eps = DBL_EPSILON, absx, dx;///FINITE_DIFFERENCE_STEPSIZE;
		double eps__ = sqrt(eps), xnew, relStep = sqrt(eps);
		
		//Richardson interpolation
		int RICHARDSON_DERIV = TRUE;
		double tp2h, tph, tmh, tm2h, b, a;
		
		gsl_vector *x_cpy = gsl_vector_alloc( x->size ) ; 
		
		int k, l;
		
		for( k=0; k<x->size; k++ ) gsl_vector_set( x_cpy, k, gsl_vector_get( x, k ) ) ; 
		
		double f0 = gsl_objective_function( x, params ), f1,f2;
		
		//double *eps_vec = calloc( 6, sizeof(double) ) ;
		//double *ieps_vec = calloc( 6,sizeof(double) ) ; 
		
		//cycle through each co-ordinate
		for( k=0; k< x->size; k++)
		{
			
			//dx = .001 ; 
			
			// Try a different dx that scales
			
			if( !RICHARDSON_DERIV )
			{
			
				dx = relStep * fabs( gsl_vector_get(x,k) ) ;
			
				xnew = gsl_vector_get( x, k)  + dx ; 
			
				gsl_vector_set( x_cpy, k, xnew ) ; 
			
				f1 = gsl_objective_function( x_cpy, params ) ;
			
				//printf("\n The value of \n \t dx = %.25f\n \t f0 = %.25f\n \t f1 = %.25f\n \t diff = %.25f", dx, f0, f1, f1-f0);
			
				gsl_vector_set( grad, k, (f1-f0)/dx  );
			
				//printf("\n grad entry %d is %lf",k,gsl_vector_get(grad,k)) ;
			
				gsl_vector_set( x_cpy, k, gsl_vector_get( x, k ) ) ;
			
			}else{
			
				absx = fabs( gsl_vector_get(x,k) );
			
				if( absx < 1.781029E-5 ) dx = 1E-4; else dx = relStep * fabs( gsl_vector_get(x,k) ) ;
				
				// x+h
				
				xnew = gsl_vector_get( x, k)  + dx ;
				
				gsl_vector_set( x_cpy, k, xnew ) ;
				
				tph = gsl_objective_function( x_cpy, params ) ;
				
				// x+2h
				
				xnew = gsl_vector_get( x, k) + 2.*dx ;
				
				gsl_vector_set( x_cpy, k, xnew );
				
				tp2h = gsl_objective_function( x_cpy, params );
				
				//x-h
				
				xnew = gsl_vector_get( x, k) - dx ;
				
				gsl_vector_set( x_cpy, k, xnew );
				
				tmh = gsl_objective_function( x_cpy, params );
			
				//x-2h
				
				xnew = gsl_vector_get( x, k) - 2.*dx ;
				
				gsl_vector_set( x_cpy, k, xnew );
				
				tm2h = gsl_objective_function( x_cpy, params );
				
				b = ( -tp2h + 8.*tph - 8.*tmh + tm2h ) ;
				
				a = 12. * dx ; 
				
				gsl_vector_set( grad, k, b/a  );
			
				gsl_vector_set( x_cpy, k, gsl_vector_get( x, k ) ) ;
			
			}
			
		}
		
		gsl_vector_free( x_cpy ) ; 
		
		//free( eps_vec );
		//free( ieps_vec ) ;
		return;
}

void gsl_objective_function_and_gradient( const gsl_vector *x, void *params, double *f, gsl_vector *grad )
{
	
	*f = gsl_objective_function( x,  params ) ;
	
	gsl_objective_function_gradient( x, params, grad ) ;
	
	return;
}

/*double gsl_matrix *gsl_hessian_at_mode( gsl_vector *x, void *params)
{
	
	struct model_gsl *par = ( struct model_gsl *)params;
	struct super *s = par->s;
	
}*/


double gsl_objective_function2( const gsl_vector *x, void *params )
{
	struct model_gsl *par = ( struct model_gsl *)params;
	struct super *s = par->s;
	
	int id = s->individual_id[ par->thread_num ];
	
	struct hyperpar *h = s->h[ id ];
	struct patch **p = s->p;
	struct patch *p0;
	int *mask;
	struct block *b = s->b;
	double *store_theta = calloc( h->npar, sizeof(double) );
	
	//printf("\n gsl_objective_function:  evaluate at  " );
	
	int k;
	for( k=0; k<4; k++ )
	{
		store_theta[k+2] = h->theta[k+2];
		h->theta[k+2] = gsl_vector_get( x, k );
	}
	
	int nbl; //number of blocks (if analyzing individually)
	if( s->individual ) nbl = 1; else nbl = b->n_block;
	
	if( s->individual )
	{
		k = id ;  
		update_patch( p[k], h, TRUE, b->mask[k] ) ;
		p0 = p[k] ;
		mask = b->mask[k];
	}
	else
	{
		update_patch( p[0], h, FALSE, b->mask[0] );
		for( k=1; k< b->n_block + 1 ; k++ ) update_patch( p[k] , h , TRUE, b->mask[k] );
		patch_combine( p , b );
		p0 = p[0] ; 
		mask = b->mask[0];
	}
	
	double log_obj = model_log_joint_marginal_hyperpar( p0, p , h , nbl, p0->masked_pix_include, mask, id ) ;
	
	for( k=0; k<4; k++ ) h->theta[k+2] = store_theta[k+2];
	free(store_theta);
	
	//printf("\n gsl_objective_function: val %lf", -log_obj );
	
	return( -log_obj ); //taking minimum here as gsl minimize the function
}

void gsl_objective_function_gradient2( const gsl_vector *x, void *params, gsl_vector *grad )
{	

		printf("\n using the  objective 2");

		double eps = DBL_EPSILON,  dx;///FINITE_DIFFERENCE_STEPSIZE;
		double eps__ = sqrt(eps), xnew;
		
		gsl_vector *x_cpy = gsl_vector_alloc( x->size ) ; 
		
		int k, l;
		
		for( k=0; k<x->size; k++ ) gsl_vector_set( x_cpy, k, gsl_vector_get( x, k ) ) ; 
		
		double f0 = gsl_objective_function2( x, params ), f1,f2;
		
		//double *eps_vec = calloc( 6, sizeof(double) ) ;
		//double *ieps_vec = calloc( 6,sizeof(double) ) ; 
		
		//cycle through each co-ordinate
		for( k=0; k< x->size; k++)
		{
			
			dx = .0001 ; 
			
			xnew = gsl_vector_get( x, k)  + dx ; 
			
			gsl_vector_set( x_cpy, k, xnew ) ; 
			
			f1 = gsl_objective_function2( x_cpy, params ) ;
			
			//printf("\n The value of \n \t dx = %.25f\n \t f0 = %.25f\n \t f1 = %.25f\n \t diff = %.25f", dx, f0, f1, f1-f0);
			
			gsl_vector_set( grad, k, (f1-f0)/dx  );
			
			//printf("\n grad entry %d is %lf",k,gsl_vector_get(grad,k)) ;
			
			gsl_vector_set( x_cpy, k, gsl_vector_get( x, k ) ) ;
			
		}
		
		gsl_vector_free( x_cpy ) ; 
		
		//free( eps_vec );
		//free( ieps_vec ) ;
		return;
}

void gsl_objective_function_and_gradient2( const gsl_vector *x, void *params, double *f, gsl_vector *grad )
{
	
	*f = gsl_objective_function2( x,  params ) ;
	
	gsl_objective_function_gradient2( x, params, grad ) ;
	
	return;
}



