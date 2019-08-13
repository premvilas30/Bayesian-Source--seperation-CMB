// grid.c

#include "grid.h"

struct grid *grid_create( int npar, int max_pts_par )
{
	int k, max_pts = max_pts_par ;
	//for( k = 1; k < npar; k++ ) max_pts *= max_pts_par ;

	struct grid *g = ( struct grid * )malloc( sizeof( struct grid ) ) ; 
	
	g->npar = npar ;
	g->max_pts = MAX_GRID_PTS;
	g->npts = 0;
	g->pts = calloc( MAX_GRID_PTS , sizeof( double * ) ) ;
	g->log_function_values = calloc( MAX_GRID_PTS , sizeof(double) ) ;
	g->weights = calloc( MAX_GRID_PTS, sizeof(double) );
	
	return( g );
}

void grid_destroy( struct grid *g )
{
	int k;
	for( k=0; k<g->npts; k++ ) free(  g->pts[ k ] ) ;  
	free( g->pts );
	free( g->log_function_values );
	free( g->weights );
	free( g );
	return;
}

struct grid_axes *grid_create_grid_axes( int npar )
{
	struct grid_axes *ga = ( struct grid_axes * )malloc( sizeof( struct grid_axes ) ) ;
	ga->e_vals = gsl_vector_alloc( npar );
	ga->sqrt_e_vals = gsl_vector_alloc( npar );
	ga->e_vecs = gsl_matrix_alloc( npar, npar ); 
	return( ga );
}

void grid_destroy_grid_axes( struct grid_axes *ga )
{
	gsl_vector_free( ga->e_vals );
	gsl_vector_free( ga->sqrt_e_vals );
	gsl_matrix_free( ga->e_vecs );
	return;
}

struct grid_axes *grid_construct_axes( struct hyperpar *h , struct patch *p )
{
	int k, l, j, npar = p->n_parinfer;
	
	struct grid_axes *ga = grid_create_grid_axes( npar );
	
	gsl_eigen_symmv_workspace *workspace;
	workspace = gsl_eigen_symmv_alloc((size_t) npar);
	
	gsl_vector *e_vals = ga->e_vals;
	gsl_vector *sqrt_e_vals = ga->sqrt_e_vals;
	gsl_matrix *e_vecs = ga->e_vecs;
	
	gsl_eigen_symmv( h->hessian, e_vals, e_vecs, workspace );
	
	gsl_eigen_symmv_free( workspace );
	
	double min_pos_e_val = DBL_MAX, eig_v;
	
	for( k=0; k < npar; k++ ) 
	{
		eig_v = gsl_vector_get( e_vals, k );
		if( eig_v > 0. ) min_pos_e_val = min_pos_e_val < eig_v ? min_pos_e_val : eig_v;
	}
	
	//if all the eigenvalues are negative set the minimum to 1.0
	if ( min_pos_e_val == DBL_MAX ) min_pos_e_val = 1.0; 
	
	int a_change = 0, all_negative = 1;

	for ( k=0; k<npar; k++) 
	{
		eig_v = gsl_vector_get( e_vals , k );

		all_negative = ( all_negative && ( eig_v <= 0.0 || eig_v == 0.0) );
		if (eig_v < 0.0) 
		{
			printf("\n Warning: Eigenvalue %1d of the Hessian is %.6g... setting to %.6g", k, eig_v, min_pos_e_val );
			gsl_vector_set( e_vals, k, min_pos_e_val );
			a_change++;
		}
		
	}
	
	//want to store the e_vals to compute marginal likelihood approximation
	
	for ( k = 0; k < npar; k++ ) 
	{
		gsl_vector_set( sqrt_e_vals, k, sqrt( gsl_vector_get(e_vals, k) ) );
		gsl_vector_set( h->e_vals, k, gsl_vector_get(e_vals, k) ) ;
	}
	
	if( a_change )
	{
		//if any of the eigenvalues have been changed, update the Hessian
		
		if( all_negative )
		{
		  //use the identity for the Hessian
		  for( k=0; k<npar; k++ )
		  {
		  	for( l=0; l<k+1 ; l++ ) 
		  	{
		  		if( l != k ) 
		  		{ 
		  			gsl_matrix_set( h->hessian, k, l, 0. ); 
		  			gsl_matrix_set(h->hessian, l, k, 0. ) ;
		  		} 
		  		else gsl_matrix_set( h->hessian, k, l, 1. );
		  	}
		  }
		}
		else
		{
			
			double sum;
			
			for( k=0; k<npar; k++ )
			{
				for( l=k; l<npar; l++ )
				{
					sum = 0.;
					for( j=0; j<npar; j++ ) sum += gsl_matrix_get( e_vecs, k, j ) * gsl_matrix_get( e_vecs, l, j ) * gsl_vector_get( e_vals, j ) ;
					gsl_matrix_set( h->hessian, k, l, sum );
					gsl_matrix_set( h->hessian, l, k, sum );
				}
			}	
				
		}
		
		/*printf("\n **** modified covariance **** \n");
		for( k=0; k<npar; k++ )
		{
			for( l=0; l<npar; l++ )
			{
				printf(" %.3f , ", gsl_matrix_get( h->hessian, k, l) );
			}
			printf("\n");
		}*/
		
	}
	
	//forgot the last part down here where we have to reupdate the e_vals and e_vecs!!!
	
	
	return( ga ) ;
}

void grid_z_to_x( double *x, int npar, double *modal_x, double *z, struct grid_axes *ga )
{
	int k, l;
	double a = 0., v, *u;
	
	u = calloc( npar, sizeof(double) );
	
	for( k=0; k<npar; k++ ) u[ k ] = z[ k ] / gsl_vector_get( ga->sqrt_e_vals, k ) ;
	
	for( k=0; k<npar; k++ )
	{
		for( l=0; l<npar; l++ )
		{
			v = gsl_matrix_get( ga->e_vecs, k, l ) ;
			a += v * u[ l ] ;
		}
		x[ k ] = modal_x[ k ] + a ;
	}
	
	free( u );
	return;
}


struct grid *grid_construct( struct super *s, int thread_num )
{
	int id;
	if( s->individual ) id = s->individual_id[ thread_num ]; else id = 0;
	
	struct hyperpar *h = s->h[ id ] ;
	int k, d, l, c, npar = s->p[id]->n_parinfer, count;
	double *delta_z = calloc( npar, sizeof(double) ), diff_thresh = 2.5;
	for( k=0; k<npar; k++ ) delta_z[k] = 1.; 
	
	gsl_vector *x = gsl_vector_alloc( npar );
	double *modal_x = calloc( npar, sizeof(double) );//gsl_vector_alloc( npar );
	
	//makes refinement to the Hessian if necessary, gets directions for grid
	struct grid_axes *ga = grid_construct_axes( h, s->p[id] );
	
	struct grid *grid = grid_create( npar, 8 ) ;
	
	grid->delta_z = delta_z[0];
	
	//create this for evaluation of the gsl_model formulation
	struct model_gsl *par = ( struct model_gsl *)malloc( sizeof( struct model_gsl ) );
	par->s = s;
	par->thread_num = thread_num;
	
	int *gmax = calloc( npar, sizeof(int) );
	int *gmaxx = calloc( npar, sizeof(int) );
	int *gmin = calloc( npar, sizeof(int) );
	int *gminn = calloc( npar, sizeof(int) );	
	
	double *z = calloc( npar, sizeof(double) );
	int *iz = calloc( npar, sizeof(int) );
	
	c=0;
	
	for( l=0; l< h->npar; l++ ) 
	{
		if( s->p[id]->parinfer[l] ) 
		{ 
			gsl_vector_set( x , c, s->h[ id ]->theta[l] ) ;
			modal_x[ c ] = s->h[id]->theta[l] ;
			c++; 
		} 
	}
	 
	
	double log_function_mode = -gsl_objective_function( x, par );
	
	for( k=0; k<npar; k++ )
	{
		//change direction in this loop
		for( d=-1; d<2 ; d+=2 )
		{
		
			memset( z, 0, npar * sizeof(double) );
			memset( iz, 0, npar * sizeof(int) ); 
			
			count = 0;
			
			while( grid->npts < MAX_GRID_PTS )
			{
			
				//if( grid->npts % 100 == 0 ) printf("\n \t on grid pt number %d", grid->npts );
			
				//increment the grid counter
			
				grid->pts[ grid->npts ] = calloc( npar, sizeof(double) ) ;
			
				z[k] += d * delta_z[k];
				
				iz[k] += d;
				
				//convert the z value to a theta value // done to here...
				grid_z_to_x( grid->pts[ grid->npts ] , npar, modal_x, z, ga ) ;
				
				for( l=0; l<npar; l++ ) gsl_vector_set( x , l, grid->pts[ grid->npts ][ l ] ) ;
				
				grid->log_function_values[ grid->npts ] = -gsl_objective_function( x, par );
				
				//printf("\n The log function value at mode %lf, at pt %lf, diff %lf", log_function_mode, grid->log_function_values[ grid->npts ] ,log_function_mode - grid->log_function_values[ grid->npts ] );
				
				if( log_function_mode - grid->log_function_values[ grid->npts ] > diff_thresh  )
				{
					 free( grid->pts[ grid->npts ] ) ;
					 break;
				}
				
				gmax[k] = gmax[k] > iz[k] ? gmax[k] : iz[k] ;
				
				gmin[k] = gmin[k] < iz[k] ? gmin[k] : iz[k] ;
				
				grid->npts += 1;
				
				count++;
				
				//printf("\n in the first part %d", grid->npts );
			
			}
		
		}
		
		//if the step size is too large, reduce by half and try again
		/*if( gmax[k] - gmin[k] + 1 < 3 ) 
		{
			delta_z[k] *= .5;
			k -= 1;
		}*/
	
	}
	
	//sfor( l=0; l<npar; l++ ) printf("\n Number of points in dir %d is %d", l, gmax[l] - gmin[l] + 1);
	
	
	//fill in the rest of the grid
	
	int *len = calloc( npar, sizeof(int) ), len_length ;
	
	for( k=0, len_length = 1; k<npar; k++)
	{
		len[k] = 1 + gmax[k] - gmin[k];
		len_length *= len[k];
	}
	
	int *iz_axes = calloc( npar, sizeof(int) );
	int *izz = calloc( npar, sizeof(int) );
	memset( izz, 0, npar * sizeof(int) );
	
	for( k=0; k<len_length; k++ )
	{
		
		//if( grid->npts % 100 == 0 ) printf("\n \t on grid pt number %d",  grid->npts );
		
		if( grid->npts == MAX_GRID_PTS ) break;
		
		for( l=0; l<npar; l++ )
		{
			iz[l] = (izz[l] <= gmax[l] ? izz[l] : gmax[l] - izz[l]);
			z[l] = iz[l] * delta_z[l];
		}
		
		grid->pts[ grid->npts ] = calloc( npar, sizeof(double) ) ;
		
		//convert the z value to a theta value 
		grid_z_to_x( grid->pts[ grid->npts ] , npar, modal_x, z, ga ) ;
		
		
		//printf("\n");
		for( l=0; l<npar; l++ ) 
		{
			gsl_vector_set( x , l, grid->pts[ grid->npts ][ l ] ) ;
			//printf(" %.3f , ", gsl_vector_get( x, l ) );
		}
				
		grid->log_function_values[ grid->npts ] = - gsl_objective_function( x, par ) ;
		
		//check if the difference from the mode is too large- discard point if so
		if( log_function_mode - grid->log_function_values[ grid->npts ] > diff_thresh ) 
		{
			free( grid->pts[ grid->npts ] );
		}
		else
		{
			grid->npts += 1;
			//if( grid->npts%200 == 0 ) printf("\n Up to %d grid pts at iteration %d ",grid->npts, k);
		}

		//compute the next configuration
		for( l=npar-1; l>=0; l-- )
		{
			if( izz[l] == ( izz[l]+1)%len[l] ) break;
			
		}
		
	}
	
	printf("\n No. grid pts is %d", grid->npts ); 
	
	//now compute the grid weights
	
	double max = -DBL_MAX, nconst = 0. ;
	k = 0;
	while( k < grid->npts )
	{
		if(  grid->log_function_values[ k ] > max )  max = grid->log_function_values[ k ] ;
		k++;
	}
	
	grid->max_log_function_value = max ;
	
	//compute the value of the approximate log marginal likelihood
	
	grid->log_marginal_likelihood = max;	
	
	
	for( k=0; k<grid->npts; k++ ) 
	{
		grid->weights[ k ] = exp( grid->log_function_values[ k ] - max ) ; 
		
		grid->log_marginal_likelihood += grid->weights[k] ; 
		
	 	nconst += grid->weights[ k ] ;
	}
	
	for( k=0; k<grid->npts; k++ ) grid->weights[ k ] /= nconst ;
	
	grid->log_marginal_likelihood *= grid->delta_z ; 
	
	//add Jacobian for the log marginal likelihood approx
	
	for( k=0; k<grid->npar; k++ ) 
		grid->log_marginal_likelihood -= 0.5 * log( gsl_vector_get( s->h[id]->e_vals, (unsigned int) k ) );
	
	grid_destroy_grid_axes( ga ) ;
	
	free( z );
	free( iz );
	gsl_vector_free( x );
	
	free( iz_axes );
	free( izz );
	free( len );
	free( gmax );
	free( gmaxx );
	free( gmin );
	free( gminn );
	
	free( modal_x );
	
	free( delta_z );
	
	return( grid );
}

void grid_compute_posterior_functionals( struct grid *g, struct patch *p, struct hyperpar *h, int *mask , int covariance, int residual )
{
	int k, l, i, j, c, nl = p->ninferred_map * p->graph->number_nodes, nd = p->nin_map * p->graph->number_nodes ;
	double *x = calloc( nl, sizeof(double) ) , *r ;
	cholmod_dense *z = cholmod_allocate_dense( nd, 1, nd, CHOLMOD_REAL, p->chol_comm );
	double **C; 
	
	if( covariance )
	{
		C = calloc( nl, sizeof(double *) );
		for( i=0; i<nl; i++ ) C[i] = calloc( nl, sizeof(double) );
	}
	
	//FILE *fp = fopen("../WMAP_result/GRID_out.txt", "a") ;
	
	if( residual ) 
	{
		r = calloc( nd, sizeof(double) ) ;
		
		//initialize the residual to the data value
		//for( i=0; i<nd; i++ ) r[i] = ( (double *) p->y->x )[i] ;
	}
	
	for( k=0; k<g->npts; k++ )
	{
		//if( k % 100 == 0 ) printf("\n Processed %d grid points", k);
	
		//put in the grid points from grid
		c = 0 ;
		for( l=0; l<h->npar; l++ )
		{
			if( p->parinfer[l] ) { h->theta[ l ] = g->pts[k][c] ; c++ ; }
		}
		
		//compute the value of the optimized field
		update_patch( p, h, TRUE, mask ) ;
		
		//add this to the overall mean
		if( p->ninferred_map > 0 )
		{
			for( i=0; i < nl; i++ ) x[ i ] +=  g->weights[ k ] * (   ( (double *) p->mu->x )[ i ]     );
		}	
		//fprintf(	fp, "%lf\n", g->weights[k] ) ;
		
		if( covariance && p->ninferred_map > 0 )
		{
			//invert the posterior precision
			update_covariance( p ) ;
			
			for( i=0; i<nl; i++ )
			{
				for( j=0; j<nl; j++ ) C[i][j] += g->weights[k] * ( ( (double *)p->Sig->x )[ i + j*nl ] 
																+  ( (double *)p->mu->x )[ i ] * ( (double *)p->mu->x )[ j ] ) ;
			}
		}
		
		if( residual )
		{
			//multiply B by the posterior MAP ordinate
			cholmod_free_dense( &z, p->chol_comm );
			
			z = update_linear_predictor( p, h ) ;
			
			//printf("\n the weight is weights[%d] = %.20f", k , g->weights[k] );
			
			//subtract the weighted predictor from the data value
			for( i=0; i<nd; i++ ) r[i] += g->weights[k] * ( (double *) z->x )[i] ;
		
		}
		
	}
	
	//fprintf( fp , "\n **************\n" ) ;
	//fclose( fp );
	
	if( residual )
	{
		//for the last time
		cholmod_free_dense( &z, p->chol_comm );
	
		//standardize the residual using the detector precision
		for( i=0; i<nd; i++ ) 
		{
			r[i] *= -1.;
			r[i] += ( (double *) p->y->x )[i] ;
			r[i] *= sqrt( ( (double *) p->C->x )[i] ) ;
		}
		
	} 
	
	//if the covariance matrix must be computed we have to do another full
	//  sweep-- this could be time consuming
	
	if( covariance )
	{
		for( k=0; k<g->npts; k++ )
		{
			c = 0;
			for( l=0; l<h->npar; l++ ) 
			{
				if( p->parinfer[l] ) { h->theta[ l ] = g->pts[k][c] ; c++; }
			}
			update_patch( p, h, TRUE, mask ) ;
			
			for( i=0; i<nl; i++ )
			{
				for( j=0; j<nl; j++ ) C[i][j] -= g->weights[k] * ( ( (double *)p->mu->x )[ i ] * x[j] + x[i] * ( ( (double *)p->mu->x )[ j ] ) ) ;
			}
		}
		
		for( i=0; i<nl; i++ ) ( (double *) p->sdev->x )[ i ] = C[i][i]; 
	}
	
	//copy from x into the return
	
	for( i=0; i < nl; i++ ) ( (double *) p->mu->x )[ i ] = x[ i ] ;
	
	//copy from r into the return
	
	if( residual )
	{
		for( i=0; i < nd; i++ ) ( (double *) p->resid->x )[ i ] = r[ i ] ; 
	}

	free( x );
	
	if( covariance )
	{
		for( i=0; i<nl; i++ ) free( C[i] );
		free( C );
	}
	
	if( residual ) free( r ) ;
	
	return;
}


void grid_print( struct grid *g, FILE *fp )
{
	int npar = g->npar, k, l;
	
	for( k=0; k<g->npts; k++ )
	{
		fprintf( fp, "%lf\t", update_sync_link( g->pts[k][0] ) );
		fprintf( fp, "%lf\t", update_gdust_link( g->pts[k][1] ) );
		for( l=2; l<npar; l++ ) fprintf( fp, "%lf\t", exp( -.5 *g->pts[k][l] ) );
		fprintf( fp, "%lf\t%lf\n", g->log_function_values[k], g->weights[k] );
	}
	
	return;

}

/************************** experimental but not used *************************************/

struct grid *grid_construct2( struct super *s, int thread_num )
{
	int id = s->individual_id[ thread_num ] ;
	struct hyperpar *h = s->h[ id ] ;
	int k, d, l, npar = 5, count;
	double delta_z = 1.;
	
	gsl_vector *x = gsl_vector_alloc( npar );
	
	//makes refinement to the Hessian if necessary, gets directions for grid
	struct grid_axes *ga = grid_construct_axes( h, s->p[id] );
	
	struct grid *grid = grid_create( npar, 5 ) ;
	
	//create this for evaluation of the gsl_model formulation
	struct model_gsl *par = ( struct model_gsl *)malloc( sizeof( struct model_gsl ) );
	par->s = s;
	par->thread_num = thread_num;
	
	int *gmax = calloc( npar, sizeof(int) );
	int *gmaxx = calloc( npar, sizeof(int) );
	int *gmin = calloc( npar, sizeof(int) );
	int *gminn = calloc( npar, sizeof(int) );	
	
	double *z = calloc( npar, sizeof(double) );
	int *iz = calloc( npar, sizeof(int) );
				
	for( l=0; l<npar; l++ ) gsl_vector_set( x , l, s->h[ id ]->theta[l] ) ; 
	
	double log_function_mode = -gsl_objective_function( x, par );
	
	for( k=0; k<npar; k++ )
	{
		//change direction in this loop
		for( d=-1; d<2 ; d+=2 )
		{
		
			memset( z, 0, npar * sizeof(double) );
			memset( iz, 0, npar * sizeof(int) ); 
			
			count = 0;
			
			while( grid->npts < MAX_GRID_PTS )
			{
			
				//increment the grid counter
			
				grid->pts[ grid->npts ] = calloc( npar, sizeof(double) ) ;
			
				z[k] += d * delta_z;
				
				iz[k] += d;
				
				//convert the z value to a theta value // done to here...
				grid_z_to_x( grid->pts[ grid->npts ] , npar, s->h[ id ]->theta, z, ga ) ;
				
				for( l=0; l<npar; l++ ) gsl_vector_set( x , l, grid->pts[ grid->npts ][ l ] ) ;
				
				grid->log_function_values[ grid->npts ] = -gsl_objective_function( x, par );
				
				//printf("\n The log function value at mode %lf, at pt %lf, diff %lf", log_function_mode, grid->log_function_values[ grid->npts ] ,log_function_mode - grid->log_function_values[ grid->npts ] );
				
				if( log_function_mode - grid->log_function_values[ grid->npts ] > 2.5  )
				{
					 free( grid->pts[ grid->npts ] ) ;
					 break;
				}
				
				gmax[k] = gmax[k] > iz[k] ? gmax[k] : iz[k] ;
				
				gmin[k] = gmin[k] < iz[k] ? gmin[k] : iz[k] ;
				
				grid->npts += 1;
				
				count++;
				
				//printf("\n in the first part %d", grid->npts );
			
			}
		
		}
	
	}
	
	//for( l=0; l<npar; l++ ) printf("\n Number of points in dir %d is %d", l, gmax[l] - gmin[l] + 1);
	
	
	//fill in the rest of the grid
	
	int *len = calloc( npar, sizeof(int) ), len_length ;
	
	for( k=0, len_length = 1; k<npar; k++)
	{
		len[k] = 1 + gmax[k] - gmin[k];
		len_length *= len[k];
	}
	
	int *iz_axes = calloc( npar, sizeof(int) );
	int *izz = calloc( npar, sizeof(int) );
	memset( izz, 0, npar * sizeof(int) );
	
	for( k=0; k<len_length; k++ )
	{
		
		if( grid->npts == MAX_GRID_PTS ) break;
		
		for( l=0; l<npar; l++ )
		{
			iz[l] = (izz[l] <= gmax[l] ? izz[l] : gmax[l] - izz[l]);
			z[l] = iz[l] * delta_z;
		}
		
		grid->pts[ grid->npts ] = calloc( npar, sizeof(double) ) ;
		
		//convert the z value to a theta value 
		grid_z_to_x( grid->pts[ grid->npts ] , npar, s->h[ id ]->theta, z, ga ) ;
		
		for( l=0; l<npar; l++ ) gsl_vector_set( x , l, grid->pts[ grid->npts ][ l ] ) ;
				
		grid->log_function_values[ grid->npts ] = - gsl_objective_function( x, par ) ;
		
		//check if the difference from the mode is too large- discard point if so
		if( log_function_mode - grid->log_function_values[ grid->npts ] > 2.5 ) 
		{
			free( grid->pts[ grid->npts ] );
		}
		else
		{
			grid->npts += 1;
		}

		//compute the next configuration
		for( l=npar-1; l>=0; l-- )
		{
			if( izz[l] = ( izz[l]+1)%len[l] ) break;
			
		}
		
	}
	
	//now compute the grid weights
	
	double max = -DBL_MAX, nconst = 0. ;
	k = 0;
	while( k < grid->npts )
	{
		if(  grid->log_function_values[ k ] > max )  max = grid->log_function_values[ k ] ;
		k++;
	}
	
	for( k=0; k<grid->npts; k++ ) 
	{
		grid->weights[ k ] = exp( grid->log_function_values[ k ] - max ) ; 
	 	nconst += grid->weights[ k ] ;
	}
	
	for( k=0; k<grid->npts; k++ ) grid->weights[ k ] /= nconst ;
	
	grid_destroy_grid_axes( ga ) ;
	
	free( z );
	free( iz );
	gsl_vector_free( x );
	
	free( iz_axes );
	free( izz );
	free( len );
	free( gmax );
	free( gmaxx );
	free( gmin );
	free( gminn );
	
	return( grid );
}




