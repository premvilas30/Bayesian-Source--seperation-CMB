//patch.c

#include "patch.h"

struct patch * patch_create( int *model, int nrow, int ncol, int nin_map, int nin_template, int nout_map, int *templated, int *mask, int all, int masked_pix_include,  int covariance, int residual, int n_specidx, int *specidx, int *parinfer, int *prior_mu, double *prior_mu_ref_freq, int *spectype )
{
	int k;
	
	struct patch *p = (struct patch *) malloc( sizeof( struct patch ) ) ;
	
	p->model = calloc( nout_map, sizeof(int) );
	
	for( k=0; k<nout_map; k++ ) p->model[k] = model[k];
	
	//printf("\n The value of model is %d", p->model );
	
	p->masked_pix_include  = masked_pix_include ;
	
	p->nin_map = nin_map ;
	
	p->nin_template = nin_template ;
	
	p->ninferred_map = nout_map - nin_template ;
	
	//p->nin_mu_prior = 0;
	//for( k=0; k<p->ninferred_map; k++ ) p->nin_mu_prior += prior_mu[k] ;
	
	p->nout_map = nout_map ;
	
	p->covariance = covariance ; 
	
	p->residual = residual ;

	p->n_specidx = n_specidx;
	
	p->specidx = calloc( 2, sizeof(int) );
	
	p->spectype = calloc( 2, sizeof(int));
	
	p->parinfer = calloc( nout_map + 2 + (nout_map-1) + nin_map , sizeof(int) );
	
	p->n_parinfer = 0;
	
	for( k=0; k< nout_map + 2 + (nout_map-1) + nin_map; k++ )
	{
		p->parinfer[k] = parinfer[k] ;
		p->n_parinfer += parinfer[k] ;
	}
	
	for( k=0; k<2; k++ ) 
	{
		p->specidx[k] = specidx[k] ;
		p->spectype[k] = spectype[k];
	}
	
	p->template_source = calloc( nout_map, sizeof(int) ) ;
	
	p->conv_ant_to_therm = calloc( nin_map, sizeof(double) );
	for( k=0; k<nin_map; k++ ) p->conv_ant_to_therm[k] = 1.;
	
	for( k=0; k<nout_map; k++ ) p->template_source[k] = templated[k] ;
	
	if( FALSE ) 
	{
		p->graph = graph_construct_2d_lattice_graph( nrow, ncol ) ;
	}
	else 
	{
		p->graph_models = graph_construct_mixed_graph_ignoring_masked_pixels( nrow, ncol, mask );
		p->graph = graph_construct_independent_pixel_graph_ignoring_masked_pixels( nrow, ncol, mask ) ; //only need this object for the number of nodes
		//if( p->model == INDPIXEL ) p->graph = graph_construct_independent_pixel_graph_ignoring_masked_pixels( nrow, ncol, mask ) ;
		//printf("\n number of nodes in graph is %d", p->graph->number_nodes ) ;
	}
	
	p->chol_comm = routines_create_cholmod_common() ;
	
	int 	nnode = p->graph->number_nodes, 
			nd  = nnode * nin_map, 
			nl = nnode * p->ninferred_map,
			nt = nnode * p->nin_template ;
	
	if( p->ninferred_map > 0 )
	{
	
		if( p->ninferred_map > 1 )
		{
			p->A = calloc( nin_map , sizeof(double *) );
			for( k=0; k<nin_map; k++ ) p->A[k] =  calloc( p->ninferred_map , sizeof(double) );
		}
		else
		{
			//just to stop any glibc
			p->A = calloc( nin_map , sizeof( double * ) );
			for( k=0; k<nin_map; k++ ) p->A[k] = calloc( 2, sizeof(double) ) ;
		}
	
	}
	
	if( p->nin_template > 0 )
	{
		if( p->nin_template > 1 )
		{
			p->F = calloc( nin_map , sizeof(double*) );
			for( k=0; k<nin_map; k++ ) p->F[k] = calloc( p->nin_template , sizeof(double) ); 
		}
		else
		{
			p->F = calloc( nin_map, sizeof(double*) );
			for( k=0; k<nin_map; k++ ) p->F[k] = calloc( 2, sizeof(double) );
		}
	}

	
	if( p->ninferred_map > 0 )
	{
		p->B = patch_set_up_B( p->graph_models[0], nin_map, p->ninferred_map , p->chol_comm );
	}
	
	if( p->nin_template > 0 )
	{
		p->G = patch_set_up_B( p->graph_models[0], nin_map, p->nin_template, p->chol_comm ) ;
	}
	
	
	p->C = cholmod_speye( nd , nd, CHOLMOD_REAL, p->chol_comm ) ; //alloc memory and initialize later
	
	p->y = cholmod_zeros( nd, 1, CHOLMOD_REAL, p->chol_comm ) ;
	
	p->z = cholmod_zeros( nl, 1, CHOLMOD_REAL, p->chol_comm ) ;
	
	p->z1 = cholmod_zeros( nd, 1, CHOLMOD_REAL, p->chol_comm ) ;
	
	if( p->ninferred_map > 0 )
	{
		p->T = graph_create_sparse_precision_matrix_template_mixed_model( p->model, p->graph_models, p->ninferred_map, masked_pix_include, mask, p->chol_comm ) ;
		p->Q = graph_create_sparse_precision_matrix_from_template( p->T, p->chol_comm ) ;
	}
	
	/*if( p->model = NNPIXEL ) 
	{
		p->T = graph_create_sparse_precision_matrix_template( p->graph, p->ninferred_map, masked_pix_include, mask, p->chol_comm) ;
		p->Q = graph_create_sparse_precision_matrix_from_template( p->T, p->chol_comm ) ;
	}
	if( p->model = INDPIXEL )
	{	
		 int dim = p->graph->number_nodes * p->ninferred_map ;
		 p->T = cholmod_speye( (size_t) dim, (size_t)  dim, CHOLMOD_REAL , p->chol_comm );
		 p->Q = cholmod_speye( (size_t) dim, (size_t)  dim, CHOLMOD_REAL , p->chol_comm );
	}*/
	
	//p->Q = graph_create_sparse_precision_matrix_from_template( p->T, p->chol_comm ) ;
	
	// p->L : symbolic factorization done after everything is first initialized
	
	//initialize this to the identity... it get's free'd later
	
	if( all ) p->Q__ = cholmod_speye( nl , nl, CHOLMOD_REAL, p->chol_comm ) ;
	
	p->mu = cholmod_zeros( nl, 1, CHOLMOD_REAL, p->chol_comm ) ;
	
	p->mu_prior = cholmod_zeros( nl, 1, CHOLMOD_REAL, p->chol_comm ); 
	
	p->mu_prior_ref_freq = calloc( nout_map, sizeof(double) );
	
	p->template_ref_freq = calloc( nout_map, sizeof(double) );
	
	p->prior_mu_source =  calloc( nout_map, sizeof(int) );
	
	for( k=0; k<p->nout_map; k++ ) 
	{
		p->prior_mu_source[k] = prior_mu[k] ;
		p->mu_prior_ref_freq[k] = prior_mu_ref_freq[k];
		p->template_ref_freq[k] = prior_mu_ref_freq[k];
	}
	
	if( p->nin_template > 0 ) p->mu_template = cholmod_zeros( nt, 1, CHOLMOD_REAL, p->chol_comm ) ;
	
	if( covariance ) 
	{
		p->Sig = cholmod_zeros( nl, nl, CHOLMOD_REAL, p->chol_comm ) ;
		//posterior standard deviation of field values
		p->sdev = cholmod_zeros( nl, 1, CHOLMOD_REAL, p->chol_comm ) ;
	}
	
	if( residual )
	{
		p->resid = cholmod_zeros( nd, 1, CHOLMOD_REAL, p->chol_comm ) ;
	}
	
	if( p->n_specidx > 0 )
	{
		if( p->specidx[0] ) p->synch_index = cholmod_allocate_dense( nd, 1, nd, CHOLMOD_REAL, p->chol_comm );
		if( p->specidx[1] ) p->dust_index = cholmod_allocate_dense( nd, 1, nd, CHOLMOD_REAL, p->chol_comm );
	}
	
	return( p ) ;
}

void patch_destroy( struct patch *p, int all  )
{
	graph_destroy( p->graph_models[0] );
	graph_destroy( p->graph_models[1] );
	free( p->graph_models );
	
	graph_destroy( p->graph );
	
	free( p->model );
	
	free( p->spectype );
	
	free( p->template_source );
	
	free( p->conv_ant_to_therm );
	
	free( p->parinfer );
	
	int k;
	for( k=0; k<p->nin_map; k++ )
	{
		if( p->ninferred_map > 0 ) free( p->A[k] );
		if( p->nin_template > 0 ) free( p->F[k] );
	}
	if( p->ninferred_map > 0 ) free( p->A );
	if( p->nin_template > 0 ) free( p->F );
	
	if( p->ninferred_map > 0 ) cholmod_free_sparse( &p->B, p->chol_comm );
	if( p->nin_template > 0 ) cholmod_free_sparse( &p->G, p->chol_comm );
	cholmod_free_sparse( &p->C, p->chol_comm );
	cholmod_free_dense( &p->y, p->chol_comm );
	cholmod_free_dense( &p->z, p->chol_comm );
	cholmod_free_dense( &p->z1, p->chol_comm );
	if( p->ninferred_map > 0 )
	{
		cholmod_free_sparse( &p->T, p->chol_comm );
		cholmod_free_sparse( &p->Q, p->chol_comm );
	}
	cholmod_free_dense( &p->mu, p->chol_comm );
	cholmod_free_dense( &p->mu_prior, p->chol_comm );
	free( p->mu_prior_ref_freq );
	free( p->template_ref_freq );
	free( p->prior_mu_source ) ;
	
	if( p->n_specidx > 0 )
	{
		if( p->specidx[0] ) cholmod_free_dense( &p->synch_index, p->chol_comm );
		if( p->specidx[1] ) cholmod_free_dense( &p->dust_index, p->chol_comm );
	}
	
	if( p->nin_template > 0 ) cholmod_free_dense( &p->mu_template, p->chol_comm );
	if( all )
	{
		cholmod_free_factor( &p->L, p->chol_comm );
		cholmod_free_sparse( &p->Q__, p->chol_comm );
	}
	if( p->covariance ) 
	{
		cholmod_free_dense( &p->Sig, p->chol_comm );
		cholmod_free_dense( &p->sdev, p->chol_comm );
	}
	if( p->residual )
	{
		cholmod_free_dense( &p->resid, p->chol_comm ) ;
	}
	routines_destroy_cholmod_common( p->chol_comm );
	free( p );
	return;
}

cholmod_sparse * patch_set_up_B( struct graph *graph , int  nin_map, int nout_map, cholmod_common *comm_ptr )
{
	int 	k, l, pix, 
			nrow = graph->number_nodes * nin_map, 
			ncol = graph->number_nodes * nout_map, 
			nz = graph->number_nodes * nin_map * nout_map ;
	
	cholmod_triplet *Btrip = cholmod_allocate_triplet( nrow, ncol, nz, 0, CHOLMOD_REAL, comm_ptr );
	
	Btrip->nnz = nz; //number non-zero entries
	
	Btrip->stype = 0; //unsymmetric
	
	int 	*ridx = (int *)Btrip->i,
			*cidx = (int *)Btrip->j,
			c = 0;
			
	//FILE *chk_B = fopen("../PLANCK_data/chk_b.txt","w");
	
	//fprintf(chk_B, "%d\t0\n", graph->number_nodes);
	//fprintf(chk_B, "%d\t%d\n", nin_map, nout_map);
	
	for( k=0; k<nin_map; k++ )
	{
		for( l=0; l<nout_map; l++ )
		{
			for( pix=0; pix<graph->number_nodes; pix++)
			{
				ridx[c] = pix + k * graph->number_nodes;
				cidx[c] = pix + l * graph->number_nodes;
				c++;
			}
		}
	}
	
	//fclose( chk_B) ; 
	
	cholmod_sparse *B = cholmod_triplet_to_sparse( Btrip, Btrip->nnz, comm_ptr );
	
	printf("\n Symmetry type is B->stype = %d ", B->stype);
	
	cholmod_free_triplet( &Btrip, comm_ptr);
	
	return( B );
}

struct patch **patch_create_from_block( int *model, struct block *b, int nin_map, int nin_template, int nout_map, int *templated, int masked_pix_include, int covariance, int residual, int n_specidx, int *specidx, int *parinfer, int *prior_mu, double *prior_mu_ref_freq, int *spectype )
{
	int k;
	
	struct patch **p = (struct patch **)malloc( ( b->n_block + 1 ) * sizeof( struct patch *)) ;
	
	p[0] = patch_create( model, b->nrow, b->ncol, nin_map, nin_template, nout_map, templated, b->mask[0], FALSE, masked_pix_include, FALSE, FALSE, n_specidx, specidx, parinfer, prior_mu, prior_mu_ref_freq, spectype );
	for( k=1; k<b->n_block+1 ; k++ ) p[k] = patch_create( model,b->nrow_subblock, b->ncol_subblock, nin_map, nin_template, nout_map, templated, b->mask[k] , TRUE, masked_pix_include, covariance, residual, n_specidx, specidx, parinfer, prior_mu, prior_mu_ref_freq, spectype ) ;
	
	return( p );
	
}

void patch_destroy_from_block( struct patch **p, struct block *b, int individual )
{
	int k;
	
	if( !individual) patch_destroy( p[0], FALSE );
	for( k=1; k<b->n_block+1; k++ ) patch_destroy( p[k], TRUE );
	free( p );
	
}

void patch_combine( struct patch **p, struct block *b )
{
	//this is a copying function that respects the ordering of the 
	// mask in all
	
	//need to allow for copying in multiple maps
	
	double *mu;
	int k, l, m;
	
	for( k=1; k<b->n_block+1; k++ )
	{
		//printf("\n combining block %d",k);
		
		mu = (double *)p[k]->mu->x;
		//copy over to hold
		for( m=0; m<p[k]->ninferred_map; m++ )
		{
			//printf("\n looping through m = %d",m);
			for( l=0; l<b->n_unmasked[k]; l++ ) 
				b->x[ b->idx[k][l] + m*b->n_pixel ] = mu[ l + m*b->n_unmasked[k] ];
		}
	}
	
	//put mean back into main patch
	mu =  (double *)p[0]->mu->x;
	for( m=0; m<p[0]->ninferred_map; m++)
	{
		for( l=0; l<b->n_pixel; l++ ) 
		{
			if( b->mask[0][l] ) 
			{
				mu[ b->reverse_idx[ l ] + m*b->n_unmasked[0] ] = b->x[ l + m*b->n_pixel ] ;
				//printf("\n mu[ %d ] = %.10f " ,b->reverse_idx[ l ] + m*b->n_unmasked[0],  b->x[ l + m*b->n_pixel ] );
			}
		}
	}
	
	/*FILE *fp;
	fp = fopen( "../WMAP_data/see_source_out_1.txt", "w");
	for( k=0; k<b->n_pixel; k++ ) fprintf( fp, "%lf\n", mu[ k ] );
	fclose( fp );
	fp = fopen( "../WMAP_data/see_source_out_2.txt", "w");
	for( k=0; k<b->n_pixel; k++ ) fprintf( fp, "%lf\n", mu[ b->n_pixel + k ] );
	fclose( fp );
	fp = fopen( "../WMAP_data/see_source_out_3.txt", "w");
	for( k=0; k<b->n_pixel; k++ ) fprintf( fp, "%lf\n", mu[ 2*b->n_pixel + k ] );
	fclose( fp );
	fp = fopen( "../WMAP_data/see_source_out_4.txt", "w");
	for( k=0; k<b->n_pixel; k++ ) fprintf( fp, "%lf\n", mu[ 3*b->n_pixel + k ] );
	fclose( fp );*/
	return;
}




