//test program for the seven year WMAP data

#include "defs.h"
#include "super.h"
#include "optimizer.h"
#include "omp.h"
#include "wmap.h"
#include "results.h"
#include "required_libs.h"
#include "grid.h"

#define OMP_DYNAMIC TRUE
#define DYNAMIC_SCHEDULE_BLK_SIZE 5

int main()
{
	int n_freq_band = 5, n_src = 4,
		 n = 512, n_sub = 64, N = n*n,  
		 k, l, k_, base, *mask, covariance = FALSE, residual = TRUE, n_templates = 0;
		 
	//timer
	time_t start, end;
	
	start = time(NULL);
	
	//construct the 
	int masked_pix_include = TRUE ;
	
	//analyse the sub-blocks individually
	int individual = TRUE; 
	int max_thread = omp_get_max_threads() ;
	int num_thread = 4;
	omp_set_num_threads( num_thread ) ;//max_thread );
	
	//seven year WMAP_data & mask
	
	char **file = wmap_get_7yr_data_file_names();
	char **rfile = wmap_get_7yr_result_file_names();
	char *mask_file = wmap_get_7yr_temperature_analysis_mask_name( );
	char **template_file = wmap_get_7yr_template_file_names();
	
	double *obs_freq = calloc( n_freq_band, sizeof(double) ), *nullptr=NULL;
	double *obs_precision = calloc( n_freq_band, sizeof(double) );
	double *spec = calloc( 2, sizeof(double) );
	double *offset = calloc( n_freq_band, sizeof(double) );
	
	wmap_get_7yr_band_frequency_and_error_precisions_and_offsets( obs_freq, obs_precision, offset );
	
	//intialize the spectral parameters to midpoint of allowable range 
	spec[0] = -3.;//SYNC_LOWER + (SYNC_UPPER-SYNC_LOWER)/2.  ;
	spec[1] = 1.7;//DUST_LOWER + (DUST_UPPER-DUST_LOWER)/2.  ;
	
	
	//optimizer
	int optrun = 3, max_iter = 50, or, window, thread_num;
	struct optpar **opt = (struct optpar **)malloc( num_thread * sizeof( struct optpar *)) ;
	
	//grid
	struct grid **grid;
	grid = (struct grid **)malloc( num_thread * sizeof( struct grid *)) ;
	
	//results
	struct results *results;
	
	int *templated = calloc( 4, sizeof(int) );	
	
	templated[1] = 0;
	templated[2] = 0;
	
	int *prior_mu = calloc( 4, sizeof(int) );
	prior_mu[0] = 1; prior_mu[1] = 1; prior_mu[2] = 1; prior_mu[3] = 1;
	double *prior_mu_ref_freq = calloc( 4, sizeof(double) );
	//reference frequencie in units of GHz
	prior_mu_ref_freq[0] = 23.;
	prior_mu_ref_freq[1] = 23;
	prior_mu_ref_freq[2] = 94.;
	prior_mu_ref_freq[3] = 23.;
	
	printf("\n check: %s", file[0] );
	printf("\n check: %s", template_file[1] );
	
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
	
	int model = NNPIXEL;
	
	//cycle through the 12 base resoulution Healpix pixels
	
	//FILE * outf = fopen( "../WMAP_result/GRID_out.txt",  "w") ;
	
	for( k=0; k<4; k++ ) n_templates += templated[k] ; 
	
	for( base = 0 ; base < 12; base++ )
	{
		
		mask = data_get_mask( base, N, mask_file ) ; 
		
		//create the super structure for this base patch 
		
		struct super *s = super_create( model, base, n, n, n_sub, n_sub, n_freq_band, n_templates, n_src, templated, mask, individual, masked_pix_include, max_thread, covariance, residual,  n_specidx, specidx, parinfer, prior_mu, prior_mu_ref_freq ) ;
		
		free( mask ) ; 
		
		super_initialize( s, file, template_file, obs_freq, obs_precision, offset, spec, 0,  nullptr, nullptr, 0 );
		
		results = results_create( N, n_freq_band, n_src, covariance, residual ) ;
		
		if( base == 0 )
		{
			//create the results  files if they don't already exist (these will be overwritten otherwise)
			results_create_fits_results_files( s->b, 21, rfile );
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
	
				s->individual_id[ thread_num ] = window;
			
				// this is the optimizer ( possibly with restarts )
			
				opt[ thread_num ] = optimizer( s, max_iter, thread_num, optrun ) ;  
		
				for(k=0; k< opt[thread_num]->npar; k++) s->h[ window ]->theta[k] = opt[ thread_num ]->opt_theta[k] ; 
				
				gsl_matrix_memcpy( s->h[ window ]->hessian , opt[ thread_num ]->hessian ) ;
				
		
				optimizer_optpar_destroy( opt[ thread_num ] ) ;
			
				FILE * outf = fopen( "../WMAP_result/THETA_out.txt",  "a") ; 
				fprintf( outf , "%lf\t" , update_sync_link( s->h[ window ]->theta[0] ) );
				fprintf( outf , "%lf\t" , update_gdust_link( s->h[ window ]->theta[1] ) );
				for( k=2; k<s->h[ window ]->npar; k++ ) fprintf( outf , "%lf\t" , s->h[ window ]->theta[k] ); 
				fprintf( outf , "\n" ); 
				fclose( outf );	
			
				//construct the grid here...
			
				grid[ thread_num ] = grid_construct( s, thread_num ) ;
				
				printf("\n Patch %d : Window %d number of grid pts %d", base, window, grid[ thread_num ]->npts ) ;
			
				//printf("\n Going into posterior computations" );
				//printf("\n Going into posterior computations" );
			
			
				// **** NEED TO FIX THE POSTERIOR FUNCTIONALS FOR NEW SET UP ****
				grid_compute_posterior_functionals( grid[ thread_num ], s->p[ window ], s->h[window], s->b->mask[ window ] , covariance, residual ) ;
			
				//printf("\n Leaving posterior computations" );
				//printf("\n Leaving posterior computations" );
			
				//for( k=0; k<n_src*N; k++) fprintf( outf, "\n%lf", ((double *) s->p[window]->sdev)[k] ) ;
			
				//grid_print( grid[ thread_num ], outf );
			
				//update_patch( s->p[ window ], s->h[ window ], TRUE, s->b->mask[ window ] );
			
				results_write_block_to_results( window, results, s->p, s->b, s->h[ window ] ) ;
			
				printf("\n Patch %d : Window  %d completed ",base, window);
			
				grid_destroy( grid[ thread_num ] ) ; 
				
			}
		
		
		}
		else
		{
			thread_num  = 0;
			
			gsl_vector *par = gsl_vector_alloc( 6 );
			 
			gsl_vector_set( par, 0 , 0. ); gsl_vector_set( par, 1 , 0. );
			gsl_vector_set( par, 2 , 0.01 ); gsl_vector_set( par, 3 , 0. );
			gsl_vector_set( par, 4 , 0.01 ); gsl_vector_set( par, 5 , 0.01 );
			
			struct model_gsl *mgsl = malloc(sizeof(struct model_gsl));
			
			mgsl->s = s;
			
			//evaluate the log posterior print and exit
			double ll = gsl_objective_function( par , (void *)mgsl );
			
			printf("\n evaluation of log posterior is %.10f", ll);
			
			gsl_vector_free( par );
			
			exit(1);
			
			opt[ 0 ] = optimizer( s, max_iter, thread_num, optrun ) ;
		
		
		}
		//fclose( outf );
		
		
		if( results_write_patch_result_to_file( base, results, s->b, rfile ) ) printf("\n Problem  writing result to file" ) ;
		
		/*FILE *out;
		
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
		super_destroy( s );
		results_destroy( results );*/
		
		printf( "\n ********* finished base patch %d ************", base );
		
	}
	
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
