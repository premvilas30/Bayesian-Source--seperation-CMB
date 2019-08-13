
#include <R.h>
#include "defs.h"
#include "super.h"
#include "optimizer.h"
#include "omp.h"
#include "wmap.h"
#include "planck.h"
#include "results.h"
#include "required_libs.h"
#include "grid.h"

// function to be called from R to do analysis of a subpatch

/*#define WMAP_7YR 0
#define WMAP_9YR 1
#define PLANCK 2*/

//do not write to file only return

void CMB_analysis_of_patch( int *which_data, int *nside_sub, int *n_blocks, int *row, int *col, int *covariance, int *residual, int *parinfer, double *post_ex_mu, double *resid, int *max_grid_pts, int *used_grid_pts, double *grid_pts , int *num_thread )
{

	int n_freq_band, n_src, nside, n_templates = 0, k, l, i, k_, m, m_, base ;
	char **file, **rfile, *mask_file, **template_file;
		
	//construct the 
	int masked_pix_include = FALSE ;
	
	//analyse the sub-blocks individually
	int individual = TRUE; 

	omp_set_num_threads( *num_thread ) ;

	int w_d = *which_data;

	switch( *which_data )
	{
		case 0:
			wmap_get_7yr_data_dims( &n_freq_band, &n_src, &nside );
			file = wmap_get_7yr_data_file_names();
			rfile = wmap_get_7yr_result_file_names();
			mask_file = wmap_get_7yr_temperature_analysis_mask_name( );
			template_file = wmap_get_7yr_template_file_names();
		break;
		
		case 1:
			wmap_get_9yr_data_dims( &n_freq_band, &n_src, &nside );
			file = wmap_get_9yr_data_file_names();
			rfile = wmap_get_9yr_result_file_names();
			mask_file = wmap_get_9yr_temperature_analysis_mask_name( );
			template_file = wmap_get_9yr_template_file_names();
		break;
		
		case 2:
			planck_get_data_dims( &n_freq_band, &n_src, &nside );
			file = planck_get_data_file_names();
			rfile = planck_get_result_file_names();
			mask_file = planck_get_temperature_analysis_mask_name( );
			//where are the template files for the Planck data?
		break;
		
	}
	
	int N = nside * nside, *prior_mu = calloc( n_src, sizeof(int) ); 
	
	double 	*obs_freq = calloc( n_freq_band, sizeof(double) ), 
				*nullptr=NULL,
				*obs_precision = calloc( n_freq_band, sizeof(double) ),
				*spec = calloc( 2, sizeof(double) ),
				*offset = calloc( n_freq_band, sizeof(double) ),
				*prior_mu_ref_freq = calloc( n_src, sizeof(double) );
	
	//get the windows that are to be analyzed
	int *Windows = calloc( *n_blocks, sizeof(int) ), windows_in_row =  nside / (*nside_sub) ;
	for( k=0; k<*n_blocks; k++ ) Windows[k] = col[k] + windows_in_row * row[k] + 1 ;

	switch( *which_data )
	{
		case 0:
			wmap_get_7yr_band_frequency_and_error_precisions_and_offsets( obs_freq, obs_precision, offset );
			wmap_get_7yr_prior_source_dims( n_src, prior_mu, prior_mu_ref_freq );
		break;
		
		case 1:
			wmap_get_9yr_band_frequency_and_error_precisions_and_offsets( obs_freq, obs_precision, offset );
			wmap_get_9yr_prior_source_dims( n_src, prior_mu, prior_mu_ref_freq );
		break;
		
		case 2:
			planck_get_band_frequency_and_error_precisions( obs_freq, obs_precision );
			//need templates for the planck data
		break;	
	}
	
	//intialize the spectral parameters to midpoint of allowable range 
	spec[0] = -3.;
	spec[1] = 1.7;
	
	
	//optimizer
	int optrun = 3, max_iter = 50, or, window, thread_num;
	struct optpar **opt = (struct optpar **)malloc( (*num_thread) * sizeof( struct optpar *) ) ;
	
	//grid
	struct grid **grid;
	grid = (struct grid **)malloc( (*num_thread) * sizeof( struct grid *)) ;
	
	//results
	struct results *results;
	
	//none of the maps use templates-- put as prior means instead
	int *templated = calloc( 4, sizeof(int) ); n_templates = 0;
	for( k=0; k<4; k++ ) n_templates += templated[k] ; 	
	
	//indicator whether to include a map of the spectral index
	int *specidx = calloc(2,sizeof(int));
	specidx[0] = 0;
	specidx[1] = 0;
	
	int n_specidx = 0;
	
	//infer all the parameters
	//int *parinfer = calloc( 6, sizeof(int) );
	int n_parinfer = 0;
	//for( k=0; k<6; k++ ) parinfer[k] = 1;
	
	for( k=0; k<6; k++ ) n_parinfer += parinfer[k];
	
	int model = NNPIXEL,  *mask;	

	for( base = 0 ; base < 1; base++ )
	{
		
		mask = data_get_mask( base, N, mask_file ) ; 
		
		//create the super structure for this base patch 
		
		struct super *s = super_create( model, base, nside, nside, *nside_sub, *nside_sub, n_freq_band, n_templates, n_src, templated, mask, individual, masked_pix_include, 4, *covariance, *residual,  n_specidx, specidx, parinfer, prior_mu, prior_mu_ref_freq ) ;
		
		free( mask ) ; 
		
		super_initialize( s, file, template_file, obs_freq, obs_precision, offset, spec, 0,  nullptr, nullptr, 0 );
		
		results = results_create( N, n_freq_band, n_src, *covariance, *residual ) ;
		
		if( base == 0 )
		{
			//create the results  files if they don't already exist (these will be overwritten otherwise)
			results_create_fits_results_files( s->b, 21, rfile );
		} 
		
		if( individual )
		{
		
			//optimization (& grid later)
		
			#pragma omp parallel for private( thread_num , or, k, i )
			for( k_=0; k_<*n_blocks; k_++ )
			{
		
				window = Windows[k_] ;
		
				thread_num = omp_get_thread_num() ; 
				
				//printf("\n beginning analysis of block %d on processor %d ", window, thread_num );
	
				s->individual_id[ thread_num ] = window;
			
				// this is the optimizer ( possibly with restarts )
			
				opt[ thread_num ] = optimizer( s, max_iter, thread_num, optrun ) ;  
		
				for(k=0; k< opt[thread_num]->npar; k++) s->h[ window ]->theta[k] = opt[ thread_num ]->opt_theta[k] ; 
				
				gsl_matrix_memcpy( s->h[ window ]->hessian , opt[ thread_num ]->hessian ) ;
				
				optimizer_optpar_destroy( opt[ thread_num ] ) ;	
			
				//construct the grid here...
			
				grid[ thread_num ] = grid_construct( s, thread_num ) ;
				
				//write the grid to the output to be returned to R
				used_grid_pts[ k_ ] =  grid[ thread_num ]->npts;
				
				for( i=0; i<used_grid_pts[ k_ ]; i++ )
				{
					// the extra +1 here in the index is for the objective function
					for( k=0; k< opt[thread_num]->npar; k++) grid_pts[ grid[ thread_num ]->max_pts * k_ + ( opt[thread_num]->npar + 1 ) * i + k  ] =  grid[ thread_num ]->pts[i][k];
					grid_pts[ grid[ thread_num ]->max_pts * k_ + ( opt[thread_num]->npar + 1 ) * i + opt[thread_num]->npar ] = grid[ thread_num ]->log_function_values[i];
				}
			
				// **** NEED TO FIX THE POSTERIOR FUNCTIONALS FOR NEW SET UP ****
				grid_compute_posterior_functionals( grid[ thread_num ], s->p[ window ], s->h[window], s->b->mask[ window ] , *covariance, *residual ) ;
			
				results_write_block_to_input_vector( k_, window, *residual, post_ex_mu, resid, s->p, s->b, s->h[ window ] ) ;
			
				grid_destroy( grid[ thread_num ] ) ; 
						
			}
			
			free( s->p );
		}
		
		//write the result back to input vectors
		
		//if( results_write_patch_result_to_file( base, results, s->b, rfile) ) printf("\n Problem  writing result to file" ) ;
		
		super_destroy( s ) ;
		results_destroy( results );
		
	}

	free( Windows );
	free( obs_freq );
	free( spec );
	free( obs_precision );
	free( offset );
	free( file );
	free( rfile );
	free( templated );
	free( prior_mu );
	free( prior_mu_ref_freq );
	
	return;
}
