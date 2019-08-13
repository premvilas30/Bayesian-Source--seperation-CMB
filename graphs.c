// graphs.c : implementation of graph functions and graph operations

#include "graphs.h"


struct graph *graph_create(int max_nodes , int max_num_neighbours )
{
	int i;
	
	struct graph *gr;
	gr = (struct graph *)malloc(sizeof(struct graph));
	
	gr->number_nodes = max_nodes;
	
	gr->number_neighbour_nodes = calloc(max_nodes,sizeof(int));
	
	gr->number_neighbour_assigned = calloc(max_nodes,sizeof(int));
	
	gr->neighbour_nodes = calloc(max_nodes,sizeof(int *));
	
	if( max_num_neighbours > 0 )
	{
		for( i=0; i<max_nodes; i++ ) gr->neighbour_nodes[i] = calloc( max_num_neighbours, sizeof(int) );
	}
	
	gr->allocated = ( max_num_neighbours > 0 ) ? TRUE : FALSE ;
	
	gr->ref_index = calloc( max_nodes, sizeof(int) );
	
	return(gr);
}


void graph_destroy( struct graph *gr )
{
	int i;
	
	//graph_print( gr, "gr_.txt");
	
	for(i=0;i<gr->number_nodes;i++)
	{
		//printf( "\n %d ", gr->number_neighbour_nodes[i] ) ;
		if( gr->number_neighbour_nodes[i] > 0 ) free( gr->neighbour_nodes[i] );
	}
	
	free( gr->neighbour_nodes );
	
	free( gr->number_neighbour_nodes );
	
	free( gr->number_neighbour_assigned );
	
	free( gr->ref_index ); 
	
	free( gr );
	
	return;
}


void graph_add_link(struct graph *gr, int node, int neighbour)
{
	if( !gr->allocated )
	{
		gr->number_neighbour_nodes[ node ] += 1;
	}
	else
	{
		gr->neighbour_nodes[ node ][ gr->number_neighbour_assigned[ node ] ] = neighbour;
		gr->number_neighbour_assigned[ node ] += 1;
	}
}


void graph_allocate_link( struct graph *gr )
{
	if( gr->allocated ) return;
	
	int i;
	for( i=0; i<gr->number_nodes; i++)
	{
		if( gr->number_neighbour_nodes[i] > 0) gr->neighbour_nodes[i] = calloc( gr->number_neighbour_nodes[i] , sizeof(int) );
	}
	
	gr->allocated = TRUE;
	
	return;
}

struct graph *graph_collapse_masked_nodes( struct graph *gr , int *mask )
{
	int unmasked = 0, i , k, n ;
	
	for( i=0; i<gr->number_nodes; i++ )  unmasked += mask[i] ;
	
	struct graph *gr_ = graph_create( unmasked , 0 );
	
	int *new_label = calloc( gr->number_nodes, sizeof( int ) ) ;
	
	i = 0;
	for( k = 0; k<gr->number_nodes; k++ )
	{
		if( mask[k] )
		{
			new_label[ k ] = i;
			i++ ; 
		}
	}
	
	n = 0;
	for( i=0; i<gr->number_nodes; i++ ) 
	{
		if( mask[i] ) 
		{
			gr_->number_neighbour_nodes[n] = gr->number_neighbour_nodes[i];
			
			gr_->neighbour_nodes[n] = calloc( gr_->number_neighbour_nodes[n] , sizeof(int) );
			
			for( k=0; k<gr->number_neighbour_nodes[i]; k++ ) gr_->neighbour_nodes[n][k] = new_label[ gr->neighbour_nodes[i][k] ] ;
			
			gr_->ref_index[n] = i;
			
			n++;
		}
		
	}
	
	gr_->allocated = TRUE ;
	
	graph_destroy( gr );
	free( new_label ) ; 
	
	return( gr_ );
}


struct graph *graph_construct_independent_pixel_graph_ignoring_masked_pixels( int nrow, int ncol, int *mask ) 
{
	int num_pixel = nrow * ncol  ;

	struct graph *gr, *grc;
	
	gr = graph_create( num_pixel , 0 );
	
	//remove masked nodes from structure and store reference to pixel position  in lattice
	grc = graph_collapse_masked_nodes( gr , mask );
	
	return( grc );
	
}






struct graph *graph_construct_exact_graph_ignoring_masked_pixels(int nrow, int ncol , int *mask )
/*constructs a graph for a 2d lattice with possibile masking of pixels. pixels numbered row major format-- pixels are completely ignored
  and a new graph is built*/
{
	
	int site, position, loop, num_pixel = nrow * ncol  ;

	struct graph *gr, *grc;
	
	gr = graph_create( num_pixel , 0 );
	
	loop = 0;
	
	while( loop <  2 )
	{
	
	for( site=0; site<num_pixel ; site++ )
	{
		//decide on where the site is
		
		if( site % ncol == 0 ) position = 0 ;
		else if( ( site + 1 ) % ncol == 0 ) position = 1 ;
		else if( site < ncol ) position = 2 ;
		else if( site > (nrow - 1)*ncol ) position = 3;
		else position = 4;
		
		//printf("\n site %d position is %d",site,position);
		
		switch( position )
		{
			case 0:
				
				if( mask[ site + 1 ]  ) graph_add_link( gr, site, site + 1 ) ;
				
				if( site > 0 && mask[ site - ncol ] ) graph_add_link( gr, site, site - ncol );
				
				if( site < (nrow - 1)*ncol && mask[ site + ncol ] ) graph_add_link( gr, site, site + ncol ); 

			break;
			
			case 1:
				
				if( mask[ site - 1 ] ) graph_add_link( gr, site, site - 1 ) ; 
				
				if( site > ncol-1 && mask[ site - ncol ] ) graph_add_link( gr, site, site - ncol ) ;
				
				if( site < num_pixel - 1 && mask[ site + ncol ] ) graph_add_link( gr, site, site + ncol ) ;
				
			break;
			
			case 2:
				
				if( mask[ site - 1 ]  ) graph_add_link( gr, site, site - 1 ) ;
				
				if( mask[ site + 1 ]  ) graph_add_link( gr, site, site + 1 ) ;
				
				if( mask[ site + ncol ]  ) graph_add_link( gr, site, site + ncol ) ;
				
			break;
			
			case 3:
				
				if( mask[ site - 1 ]  ) graph_add_link( gr, site, site - 1 ) ;
				
				if( mask[ site + 1 ]  ) graph_add_link( gr, site, site + 1 ) ;
				
				if( mask[ site - ncol ]  ) graph_add_link( gr, site, site - ncol ) ;
				
			break;
			
			case 4:
				
				if( mask[ site - 1 ]  ) graph_add_link( gr, site, site - 1 ) ;
				
				if( mask[ site + 1 ]  ) graph_add_link( gr, site, site + 1 ) ;
				
				if( mask[ site - ncol ]  ) graph_add_link( gr, site, site - ncol ) ;
				
				if( mask[ site + ncol ]  ) graph_add_link( gr, site, site + ncol ) ;
			
			break;
			
		}
		
	}
	
	graph_allocate_link( gr ); 
	
	loop++;
	
	}
	
	//remove masked nodes from structure and store reference to pixel position  in lattice
	grc = graph_collapse_masked_nodes( gr , mask );
	
	
	return( grc );

}


//return both the nearest neighbour and independence graph here

struct graph **graph_construct_mixed_graph_ignoring_masked_pixels(int nrow, int ncol , int *mask )
{
		
	struct graph **grphs = (struct graph **)malloc( 2 * sizeof(struct graph *) ) ;
	
	//model 0 is the 
	grphs[0] = graph_construct_independent_pixel_graph_ignoring_masked_pixels( nrow, ncol, mask ) ;
	grphs[1] = graph_construct_exact_graph_ignoring_masked_pixels( nrow, ncol , mask ) ;
	
	return( grphs );

}


struct graph *graph_construct_2d_lattice_graph(int nrow, int ncol )
/*constructs a graph for a 2d lattice with possibile masking of pixels. pixels numbered row major format*/
{
	
	int site, position, loop, num_pixel = nrow * ncol  ;

	struct graph *gr, *grc;
	
	gr = graph_create( num_pixel , 0 );
	
	loop = 0;
	
	while( loop <  2 )
	{
	
	for( site=0; site<num_pixel ; site++ )
	{
		//decide on where the site is
		
		if( site % ncol == 0 ) position = 0 ;
		else if( ( site + 1 ) % ncol == 0 ) position = 1 ;
		else if( site < ncol ) position = 2 ;
		else if( site > (nrow - 1)*ncol ) position = 3;
		else position = 4;
		
		//printf("\n site %d position is %d",site,position);
		
		switch( position )
		{
			case 0:
				
				graph_add_link( gr, site, site + 1 ) ;
				
				if( site > 0 ) graph_add_link( gr, site, site - ncol );
				
				if( site < (nrow - 1)*ncol ) graph_add_link( gr, site, site + ncol ); 

			break;
			
			case 1:
				
				graph_add_link( gr, site, site - 1 ) ; 
				
				if( site > ncol-1 ) graph_add_link( gr, site, site - ncol ) ;
				
				if( site < num_pixel - 1 ) graph_add_link( gr, site, site + ncol ) ;
				
			break;
			
			case 2:
				
				graph_add_link( gr, site, site - 1 ) ;
				
				graph_add_link( gr, site, site + 1 ) ;
				
				graph_add_link( gr, site, site + ncol ) ;
				
			break;
			
			case 3:
				
				graph_add_link( gr, site, site - 1 ) ;
				
			   graph_add_link( gr, site, site + 1 ) ;
				
				graph_add_link( gr, site, site - ncol ) ;
				
			break;
			
			case 4:
				
				graph_add_link( gr, site, site - 1 ) ;
				
				graph_add_link( gr, site, site + 1 ) ;
				
				graph_add_link( gr, site, site - ncol ) ;
				
				graph_add_link( gr, site, site + ncol ) ;
			
			break;
			
		}
		
	}
	
	graph_allocate_link( gr ); 
	
	loop++;
	
	}
	
	//remove masked nodes from structure and store reference to pixel position  in lattice
	//grc = graph_collapse_masked_nodes( gr , mask );
	
	
	return( gr );

}


cholmod_sparse *graph_create_sparse_precision_matrix_template( struct graph *gr, int nout_map , int masked_pix_include, int *mask, cholmod_common *comm_ptr )
{
	
	int 	nz = gr->number_nodes, //number of non-zero  entries in the matrix
			k , l, m;
	
	if( masked_pix_include )
	{
		for( k=0; k<gr->number_nodes; k++ ) nz += mask[k] * gr->number_neighbour_nodes[k] ;
	}
	else
	{
		for( k=0; k<gr->number_nodes; k++ ) nz += gr->number_neighbour_nodes[k] ;
	}
	
	cholmod_triplet *Ttrip =  cholmod_allocate_triplet( gr->number_nodes * nout_map, gr->number_nodes* nout_map, nz * nout_map , 0, CHOLMOD_REAL, comm_ptr ); 
	
	Ttrip->nnz = nz * nout_map;
	
	int	 *ridx = (int *)Ttrip->i,
			 *cidx = (int *)Ttrip->j,
			 c = 0;
	double *x = (double *)Ttrip->x;
	
	for( k=0; k<gr->number_nodes; k++ )
	{
		//put arg in here to handle the mask
		for( m=0; m<nout_map; m++ )
		{
			//diagonal entries
			ridx[c] = k + m * gr->number_nodes;
			cidx[c] = k + m * gr->number_nodes;
			if( masked_pix_include ) 
			{
				if( mask[k] ) x[c] = (double) gr->number_neighbour_nodes[k]; else x[c] = 1000.;
			}
			else
			{
				x[c] = (double) gr->number_neighbour_nodes[k];
			}
			c++;
			//off-diagonal entries
			
			
			for( l=0; l<gr->number_neighbour_nodes[k]; l++ )
			{
				ridx[c] = k + m * gr->number_nodes;
				cidx[c] = gr->neighbour_nodes[k][l] + m * gr->number_nodes ;
				if( masked_pix_include )
				{
					if( mask[k] ) { x[c] =  -1.;  c++; } 
				}
				else
				{
					x[c] = -1.;	c++;
				}
				
			}
			
			
		}
		
	}
	
	//cholmod_print_triplet( Ttrip, "Ttrip", comm_ptr );
	
	cholmod_sparse *Tunsymm = cholmod_triplet_to_sparse( Ttrip, nz, comm_ptr );	
	
	cholmod_free_triplet( &Ttrip, comm_ptr );
	
	cholmod_sparse *T = cholmod_copy( Tunsymm, -1, 1, comm_ptr );
	
	cholmod_free_sparse( &Tunsymm, comm_ptr );
	
	return(T);
}


cholmod_sparse *graph_create_sparse_precision_matrix_template_independent_pixel( struct graph *gr, int nout_map , int masked_pix_include, int *mask, cholmod_common *comm_ptr )
{
	
	int 	nz = gr->number_nodes, //number of non-zero  entries in the matrix
			k , l, m;
	
	cholmod_triplet *Ttrip =  cholmod_allocate_triplet( gr->number_nodes * nout_map, gr->number_nodes* nout_map, nz * nout_map , 0, CHOLMOD_REAL, comm_ptr ); 
	
	Ttrip->nnz = nz * nout_map;
	
	int	 *ridx = (int *)Ttrip->i,
			 *cidx = (int *)Ttrip->j,
			 c = 0;
	double *x = (double *)Ttrip->x;
	
	for( k=0; k<gr->number_nodes; k++ )
	{
		//put arg in here to handle the mask
		for( m=0; m<nout_map; m++ )
		{
			//diagonal entries
			ridx[c] = k + m * gr->number_nodes;
			cidx[c] = k + m * gr->number_nodes;
			if( masked_pix_include ) 
			{
				if( mask[k] ) x[c] = 1.; else x[c] = 1000.;
			}
			else
			{
				x[c] = 1.;
			}
			c++;
		}
		
	}
	
	//cholmod_print_triplet( Ttrip, "Ttrip", comm_ptr );
	
	cholmod_sparse *Tunsymm = cholmod_triplet_to_sparse( Ttrip, nz, comm_ptr );	
	
	cholmod_free_triplet( &Ttrip, comm_ptr );
	
	cholmod_sparse *T = cholmod_speye( gr->number_nodes * nout_map, gr->number_nodes * nout_map, CHOLMOD_REAL , comm_ptr );//cholmod_copy( Tunsymm, -1, 1, comm_ptr );
	
	cholmod_free_sparse( &Tunsymm, comm_ptr );
	
	return(T);
}



cholmod_sparse *graph_create_sparse_precision_matrix_from_template( cholmod_sparse *T, cholmod_common *comm_ptr )
{
	
	//cholmod_sparse *Q = cholmod_speye( T->nrow, T->ncol, CHOLMOD_REAL , comm_ptr );//cholmod_copy( T , -1, 1, comm_ptr );
	cholmod_sparse *Q = cholmod_copy_sparse( T , comm_ptr );
	
	return( Q );
	
}


void graph_print( struct graph *gr , char *file )
{
	int i,k;
	FILE *fp;
	
	fp = fopen( file , "w" );
	//print a list of the links
	for( i=0; i<gr->number_nodes; i++ )
	{
		for( k=0; k<gr->number_neighbour_nodes[i]; k++)
			fprintf( fp, "%d \t %d \n", i, gr->neighbour_nodes[i][k] );
	}
	fclose(fp);
	
	return;
}

//-------------------------------------- construct precision matrix for mixed objects ---------------------------------- //


cholmod_sparse *graph_create_sparse_precision_matrix_template_mixed_model( int *model, struct graph **graph_models , int nout_map, int masked_pix_include, int *mask, cholmod_common *comm_ptr )
{
	
	int 	nz = 0, //initialize the number of non-zero entries in the graph
			nnodes = graph_models[0]->number_nodes, //this is the same regardless the model
			k , l, m;
	
	for( m=0; m<nout_map; m++ )
	{
		if( model[m] == 0 )
		{
			//just the nodes themselves, no neighbours
			nz += nnodes ;
		}
		else
		{
			//add the nodes 
			nz += nnodes ;
			//add the neighbours (off diagonals)
			for( k=0; k<nnodes; k++ ) nz += graph_models[1]->number_neighbour_nodes[k] ;
		}
	}
	
	cholmod_triplet *Ttrip =  cholmod_allocate_triplet( nnodes * nout_map, nnodes * nout_map, nz , 0, CHOLMOD_REAL, comm_ptr ); 
	
	Ttrip->nnz = nz ;
	
	int	 *ridx = (int *)Ttrip->i,
			 *cidx = (int *)Ttrip->j,
			 c = 0;
	double *x = (double *)Ttrip->x;
	
	for( m=0; m<nout_map; m++ )
	{
	
		for( k=0; k<nnodes; k++ )
		{
		
			ridx[c] = k + m * nnodes;
			cidx[c] = k + m * nnodes;
			if( model[m] == 0 ) x[c] = 1.; else x[c] = (double) graph_models[1]->number_neighbour_nodes[k];
			c++;
			
			if( model[m] == 1 )
			{
				//add the off diagonal entries here
				
				for( l=0; l<graph_models[1]->number_neighbour_nodes[k]; l++ )
				{
					ridx[c] = k + m * nnodes;
					cidx[c] = graph_models[1]->neighbour_nodes[k][l] + m * nnodes ;
					x[c] = -1.;
					c++;
				}
	
			}
			
		}
		
	}
	
	//printf("\n ** Number nodes ** %d ", nnodes );
	
	/*FILE *fp;
	fp = fopen("../WMAP_result/PrecisionTestTrip.txt", "w");
	for( k=0; k<nz; k++ )
	{
		fprintf(fp, "%d\t%d\t%lf\n", ridx[k], cidx[k], x[k]);
	}
	fclose( fp );*/
	
	//cholmod_print_triplet( Ttrip, "Ttrip", comm_ptr );
	
	cholmod_sparse *Tunsymm = cholmod_triplet_to_sparse( Ttrip, nz, comm_ptr );	
	
	cholmod_free_triplet( &Ttrip, comm_ptr );
	
	cholmod_sparse *T = cholmod_copy( Tunsymm, -1, 1, comm_ptr );
	
	cholmod_free_sparse( &Tunsymm, comm_ptr );
	
	return(T);
}


