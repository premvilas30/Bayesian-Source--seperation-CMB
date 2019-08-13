//test program for the seven year WMAP data

#include "defs.h"
#include "super.h"
#include "model.h"
#include "optimizer.h"
#include "omp.h"
#include "update.h"
#include "wmap.h"
#include "results.h"

#define OMP_DYNAMIC TRUE
#define DYNAMIC_SCHEDULE_BLK_SIZE 3

int main()
{
	int n_freq_band = 5, n_src = 4,
		 n = 512, nsub = 64, N = n*n, 
		 k, l, k_, window, *mask;
	
	//construct the 
	int masked_pix_include = FALSE ;
	
	//analyse the sub-blocks individually
	int individual = TRUE; 
	
	int max_thread = omp_get_max_threads() ;
	omp_set_num_threads( max_thread );
	
	//seven year WMAP_data & mask
	
	char **file = wmap_get_7yr_data_file_names();
	char **rfile = wmap_get_7yr_result_file_names();
	char *mask_file = wmap_get_7yr_temperature_analysis_mask_name( );
	
	int *mask = calloc( N , sizeof(int) );

	data_read_mask(  myid , N, mask, mask_file ) ; 
	
	double *obs_freq = calloc( n_freq_band, sizeof(double) ), *nullptr=NULL;
	double *obs_precision = calloc( n_freq_band, sizeof(double) );
	double *spec = calloc( 2, sizeof(double) );
	
	wmap_get_7yr_band_frequency_and_error_precisions( obs_freq, obs_precision );
	
	//intialize the spectral parameters to midpoint of allowable range 
	spec[0] = SYNC_LOWER + (SYNC_UPPER-SYNC_LOWER)/2.  ;
	spec[1] = DUST_LOWER + (DUST_UPPER-DUST_LOWER)/2.  ;
	
	//create the results  files if they don't already exist (these will be overwritten otherwise)
	results_create_fits_results_files( s->b, 2*n_src+2, rfile );
	
	//optimizer
	int optrun = 2, max_iter = 25, or, window, thread_num;
	struct optpar **opt = (struct optpar **)malloc( max_thread * sizeof( struct optpar *)) ;
	
	//results
	struct results *results;
	
	//cycle through the 12 base resoulution Healpix pixels
	
	for( base = 0; base < 12; base++ )
	{
		
		mask = data_get_mask( base, N, mask_file ) ; 
		
		results = results_create( N, n_src ) ;
		
		//create the super structure for this base patch 
		
		struct super *s = super_create( myid, n, n, n_sub, n_sub, n_freq_band, n_src, mask, individual, masked_pix_include, max_thread ) ;
		
		super_initialize( s, file, obs_freq, obs_precision, src_precision, spec, 0,  nullptr, nullptr );
		
		//optimization (& grid later)
		
		#pragma omp parallel for private( thread_num , or, k )
		for( window = 1; window < s->b->n_block + 1 ; window++ )
		{
		
			thread_num = omp_get_thread_num() ; 
	
			s->individual_id[ thread_num ] = window;
			
			for( or = 0 ; or < optrun; or++ )
			{
				 opt[ thread_num ] = optimizer( s, max_iter, thread_num ) ; 
		
				for(k=0; k< opt[thread_num]->npar; k++) s->h[ window ]->theta[k] = opt[ thread_num ]->opt_theta[k] ; 
		
				optimizer_optpar_destroy( opt[ thread_num ] ) ;
			
				update_patch( s->p[ window ], s->h[ window ], TRUE, s->b->mask[ window ] );
			
				results_write_block_to_results( window, results, s->p, s->b, s->h[ window ] ) ;
			
			}
		
			printf("\n Patch %d : Window  %d completed ",base, window);
		}	
		
		
		char *file_name = "../TEST_data/patch_mean_out.txt";
		out = fopen( file_name, "w" );
		for( k=0; k<N; k++ ) fprintf( out, "%lf\n", results->mu[0][k] );
		fclose( out );
	
		file_name = "../TEST_data/patch_sync_ind_out.txt";
		out = fopen( file_name, "w" );
		for( k=0; k<N; k++ ) fprintf( out, "%lf\n", results->sync_ind[k] );
		fclose( out );
	
		file_name = "../TEST_data/patch_dust_ind_out.txt";
		out = fopen( file_name, "w" );
		for( k=0; k<N; k++ ) fprintf( out, "%lf\n", results->dust_ind[k] );
		fclose( out );
	
		file_name = "../TEST_data/patch_precision_out.txt";
		out = fopen( file_name, "w" );
		for( k=0; k<N; k++ ) fprintf( out, "%lf\n", results->precision[0][k] );
		fclose( out );	
		
		//free things to be allocated fresh for base++
		free( mask ) ; 
		super_destroy( s );
		results_destroy( results );
		
		printf( "\n ********* finished base patch %d ************", base );
		
	}
	
	free( obs_freq );
	free( spec );
	free( obs_precision );
	free( file );
	free( rfile );
	
	return(1);
}
