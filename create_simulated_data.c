//test program for the seven year WMAP data

#include "defs.h"
#include "super.h"
#include "optimizer.h"
#include "omp.h"
#include "wmap.h"
#include "simulated_data.h"
#include "results.h"
#include "required_libs.h"
#include "grid.h"

#define OMP_DYNAMIC TRUE
#define DYNAMIC_SCHEDULE_BLK_SIZE 5

int main()
{
	int Npatch = 12;

	int n_freq_band = 5, n_src = 4,
		 n = 512, n_sub = 128, N = n*n,  
		 k, l, k_, base, *mask, covariance = FALSE, residual = FALSE, n_templates = 0;
		 
	//timer
	time_t start, end;
	
	start = time(NULL);
	
	//construct the 
	int masked_pix_include = TRUE ;
	
	//analyse the sub-blocks individually
	int individual = TRUE; 
	int max_thread = omp_get_max_threads() ;
	int num_thread = 1;
	omp_set_num_threads( num_thread ) ;//max_thread );
	
	// files to print out to 
	char **file = wmap_get_9yr_data_file_names();
	char **rfile = simulated_data_data_file_names();
	// files giving the sources
	char **template_file = wmap_get_9yr_template_file_names();
	
	double *obs_freq = calloc( n_freq_band, sizeof(double) ), *nullptr=NULL;
	double *obs_precision = calloc( n_freq_band, sizeof(double) );
	double *spec = calloc( 2, sizeof(double) );
	double *offset = calloc( n_freq_band, sizeof(double) );
	
	wmap_get_9yr_band_frequency_and_error_precisions_and_offsets( obs_freq, obs_precision, offset );
	
	//intialize the spectral parameters to midpoint of allowable range 
	spec[0] = -3.0;//SYNC_LOWER + (SYNC_UPPER-SYNC_LOWER)/2.  ;-2.9
	spec[1] = 1.7;//DUST_LOWER + (DUST_UPPER-DUST_LOWER)/2.  ;
	
	
	//optimizer
	int optrun = 3, max_iter = 50, or, window, thread_num;
	struct optpar **opt = (struct optpar **)malloc( num_thread * sizeof( struct optpar *)) ;
	
	//results
	struct results *results;
	
	int *templated = calloc( 4, sizeof(int) );	
	
	templated[1] = 1;
	templated[2] = 1;
	templated[3] = 1;
	
	int *prior_mu = calloc( 4, sizeof(int) );
	prior_mu[0] = 0;// prior_mu[1] = 1; prior_mu[2] = 1; prior_mu[3] = 1;
	double *prior_mu_ref_freq = calloc( 4, sizeof(double) );
	//reference frequencie in units of GHz
	prior_mu_ref_freq[0] = 23.;
	prior_mu_ref_freq[1] = 23;
	prior_mu_ref_freq[2] = 94.;
	prior_mu_ref_freq[3] = 23.;
	
	//printf("\n check: %s", file[0] );
	//printf("\n check: %s", template_file[1] );
	
	//indicator whether to include a map of the spectral index
	int *specidx = calloc(2,sizeof(int));
	specidx[0] = 0;
	specidx[1] = 0;
	
	int n_specidx = 0;
	
	int *parinfer = calloc( 6, sizeof(int) ), n_parinfer = 0;
	for( k=0; k<6; k++ ) parinfer[k] = 1;
	
	/*parinfer[0] = 1;
	parinfer[1] = 0;
	parinfer[2] = 1;*/
	
	for( k=0; k<6; k++ ) n_parinfer += parinfer[k];
	
	int *model = calloc( 4, sizeof(int) );
	model[0] = 1; model[1] = 0; model[2] = 0; model[3] = 0;
	
	//cycle through the 12 base resoulution Healpix pixels
	
	//FILE * outf = fopen( "../WMAP_result/GRID_out.txt",  "w") ;
	
	for( k=0; k<4; k++ ) n_templates += templated[k] ; 
	
	for( base = 0 ; base < Npatch; base++ )
	{
		
		mask = data_get_mask_ones( N ) ; 
		
		//create the super structure for this base patch 
		
		struct super *s = super_create( model, base, n, n, n_sub, n_sub, n_freq_band, n_templates, n_src, templated, mask, individual, masked_pix_include, max_thread, covariance, residual,  n_specidx, specidx, parinfer, prior_mu, prior_mu_ref_freq ) ;
		
		free( mask ) ; 
		
		super_initialize( s, file, template_file, obs_freq, obs_precision, offset, spec, 0,  nullptr, nullptr, 0 );
		
		results = results_create( N, 0 , n_freq_band, FALSE, FALSE ) ;
		
		if( base == 0 )
		{
			//create the results  files if they don't already exist (these will be overwritten otherwise)
			results_create_fits_results_files( s->b, 5, rfile );
		} 
		
		if( individual )
		{
		
			//optimization (& grid later)
			
			//update_patch( s->p[1], s->h[1], TRUE, s->b->mask[1] );
			
			//results_write_block_to_results( 1, results, s->p, s->b, s->h[ 1 ] ) ;
		
			#pragma omp parallel for private( thread_num , or, k )
			for( window = 1; window < s->b->n_block + 1 ; window++ )
			{
		
				thread_num = omp_get_thread_num() ; 
				
				printf("\n beginning analysis of block %d on processor %d ", window, thread_num );
	
				s->individual_id[ thread_num ] = window;
			
				update_simulate_data( s->p[window], s->h[window], mask, 0 );
			
				results_write_block_to_results( window, results, s->p, s->b, s->h[ window ], TRUE );
				
				patch_destroy( s->p[ window ] , TRUE );
				
				printf("\n beginning analysis of block %d on processor %d ", window, thread_num );
				
			}
			
			free( s->p );
		}
		
		
		if( results_write_patch_result_to_file_2( base, results, s->b, rfile ) ) printf("\n Problem  writing result to file" ) ;
		
		super_destroy( s ) ;
		
		printf( "\n ********* finished base patch %d ************", base );
		
	}
	
	results_simulated_data_copy_hitrate( N, n_freq_band, file, rfile );
	
	end = time(NULL);
	
	//for( k=0; k<num_thread; k++ ) grid_destroy( grid[k] );
	//free( grid );
	
	printf("\n\n Time to process all sky on %d processors was %.3f hr \n\n", max_thread,  (double)(end-start)/3600.00);
	
	free( obs_freq );
	free( spec );
	free( obs_precision );
	free( offset );
	free( file );
	free( rfile );
	free( templated );
	free( prior_mu );
	free( prior_mu_ref_freq );
	
	return(1);
}
