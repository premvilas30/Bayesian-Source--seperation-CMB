// graphs.h
// functions to create and manage graphs and graph operations

#ifndef _GRAPHS_H_
#define _GRAPHS_H_

#include "defs.h"
#include "required_libs.h"

//structures

struct graph
{
	int number_nodes;
	int **neighbour_nodes;
	int *number_neighbour_nodes;
	int *number_neighbour_assigned;
	int *ref_index;
	int allocated;
};



//functions

struct graph *graph_create(int max_nodes, int  max_num_neighbours );

void graph_destroy( struct graph *gr );

void graph_add_link(struct graph *gr, int node, int neighbour);

void graph_allocate_link( struct graph *gr );

struct graph *graph_collapse_masked_nodes( struct graph *gr , int *mask );

struct graph *graph_construct_independent_pixel_graph_ignoring_masked_pixels( int nrow, int ncol, int *mask ) ;

struct graph *graph_construct_exact_graph_ignoring_masked_pixels(int nrow, int ncol , int *mask );

struct graph *graph_construct_2d_lattice_graph(int nrow, int ncol );

struct graph **graph_construct_mixed_graph_ignoring_masked_pixels(int nrow, int ncol , int *mask ) ;

cholmod_sparse *graph_create_sparse_precision_matrix_template( struct graph *gr, int nout_map , int use_mask, int *mask, cholmod_common *comm_ptr );

cholmod_sparse *graph_create_sparse_precision_matrix_template_independent_pixel( struct graph *gr, int nout_map , int masked_pix_include, int *mask, cholmod_common *comm_ptr );

cholmod_sparse *graph_create_sparse_precision_matrix_from_template( cholmod_sparse *T, cholmod_common *comm_ptr );

void graph_print( struct graph *gr , char *file );

cholmod_sparse *graph_create_sparse_precision_matrix_template_mixed_model( int *model, struct graph **graph_models , int nout_map, int masked_pix_include, int *mask, cholmod_common *comm_ptr );

#endif
