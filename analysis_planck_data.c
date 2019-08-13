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
	int Npatch = 12;

	int n_freq_band = 9, n_src = 2,
		 n = 1024, n_sub = 32, N = n*n,  
		 k, l, c, k_, base, *mask, covariance = FALSE, residual = TRUE, monopole = FALSE, n_templates = 0;
		 
	//timer
	time_t start, end;
	
	start = time(NULL);
	
	//construct the 
	int masked_pix_include = FALSE ;
	
	//analyse the sub-blocks individually
	int individual = TRUE; 
	int max_thread = omp_get_max_threads() ;
	int num_thread = 3;
	omp_set_num_threads( num_thread ) ;//max_thread );
	
	//nine year WMAP_data & mask
	
	char **file = planck_get_data_file_names();
	char **rfile = planck_get_result_file_names();
	char *mask_file = planck_get_temperature_analysis_mask_name( );
	char **template_file = planck_get_template_file_names();
	
	double *obs_freq = calloc( n_freq_band, sizeof(double) ), *nullptr=NULL;
	double *obs_precision = calloc( n_freq_band, sizeof(double) );
	double *spec = calloc( 2, sizeof(double) );
	double *offset = calloc( n_freq_band, sizeof(double) );
	
	planck_get_band_frequency_and_error_precisions( obs_freq, obs_precision);
	
	//intialize the spectral parameters to midpoint of allowable range 
	spec[0] = 1.7;//SYNC_LOWER + (SYNC_UPPER-SYNC_LOWER)/2.  ;
	//spec[1] = 1.7;//DUST_LOWER + (DUST_UPPER-DUST_LOWER)/2.  ;
	
	
	//optimizer
	int optrun = 2, max_iter = 50, or, window, thread_num;
	struct optpar **opt = (struct optpar **)malloc( num_thread * sizeof( struct optpar *)) ;

	//par
	struct model_gsl **par = ( struct model_gsl **)malloc( num_thread * sizeof( struct model_gsl *) );	
	for( k=0; k<num_thread; k++ ) par[k] = (struct model_gsl *) malloc( sizeof( struct model_gsl ) ) ;
	
	//grid
	struct grid **grid;
	grid = (struct grid **)malloc( num_thread * sizeof( struct grid *)) ;
	
	//results
	struct results *results;
	
	int *templated = calloc( 4, sizeof(int) );	
	/*templated[0] = 0; 
	templated[1] = 0; 
	templated[2] = 1; 
	templated[3] = 1;*/
	
	//templated[1] = 0;
	
	int *prior_mu = calloc( 4, sizeof(int) );
	prior_mu[0] = 0; 
	prior_mu[1] = 0; 
	prior_mu[2] = 0; 
	prior_mu[3] = 0;
	double *prior_mu_ref_freq = calloc( 4, sizeof(double) );
	//reference frequencie in units of GHz
	prior_mu_ref_freq[0] = 30.;
	prior_mu_ref_freq[1] = 30;
	prior_mu_ref_freq[2] = 30.;
	prior_mu_ref_freq[3] = 30.;
	
	//printf("\n check: %s", file[0] );
	//printf("\n check: %s", template_file[1] );
	
	//indicator whether to include a map of the spectral index
	int *specidx = calloc(2,sizeof(int));
	specidx[0] = 0;
	specidx[1] = 0;
	
	int *spectype = calloc( 2, sizeof(int));
	spectype[0] = SPECPAR_GEN;
	
	int n_specidx = 0;
	
	int *parinfer = calloc( 14, sizeof(int) ), n_parinfer = 0;
	
	//spectral params and CMB precision
	for( k=0; k<3; k++ ) parinfer[k] = 1;
	parinfer[1] = 0;
	
	//precisions for other sources
	for( k=3; k<4; k++ ) parinfer[k] = 1 - templated[k-2];
	
	//overall amplitude  for sync, dust, ffem
	parinfer[6] = 0;
	parinfer[7] = 0;
	parinfer[8] = 0;
	
	
	for( k=0; k<14; k++ ) 
	{
		n_parinfer += parinfer[k];
		printf("\n parinfer[%d] = %d", k, parinfer[k]);
	}
	
	//printf("\n the number  of inferred params is %d", n_parinfer);
	
	int *model = calloc( 4, sizeof(int) );
	model[0] = 1; model[1] = 1; model[2] = 0; model[3] = 0;
	//NNPIXEL;
	
	double log_marg_like = 0., log_marg_1;
	
	//cycle through the 12 base resoulution Healpix pixels
	
	//FILE * outf = fopen( "../WMAP_result/GRID_out.txt",  "w") ;
	
	for( k=0; k<4; k++ ) n_templates += templated[k] ;
	
	//printf("\n **** number templates %d *****\n", n_templates ); 
	
	for( base = 0 ; base < Npatch; base++ )
	{
		
		mask = data_get_mask( base, N, mask_file ) ; 
		
		//create the super structure for this base patch 
		
		struct super *s = super_create( model, base, n, n, n_sub, n_sub, n_freq_band, n_templates, n_src, templated, mask, individual, masked_pix_include, max_thread, covariance, residual,  n_specidx, specidx, parinfer, prior_mu, prior_mu_ref_freq, spectype ) ;
		
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
		
			log_marg_1 = 0;
		
			//optimization (& grid later)
			
			//update_patch( s->p[1], s->h[1], TRUE, s->b->mask[1] );
			
			//results_write_block_to_results( 1, results, s->p, s->b, s->h[ 1 ] ) ;
		
			#pragma omp parallel for private( thread_num , or, k ) reduction( + : log_marg_1 )
			for( window = 1; window < s->b->n_block + 1 ; window++ )
			{
		
				//only analyze the window if there are a sufficient number of points unmasked
				
				if( s->p[window]->graph->number_nodes > 0 )
				{
		
				thread_num = omp_get_thread_num() ; 
				
				//printf("\n beginning analysis of block %d on processor %d ", window, thread_num );
	
				s->individual_id[ thread_num ] = window;
			
				// this is the optimizer ( possibly with restarts )
			
				opt[ thread_num ] = optimizer( s, max_iter, thread_num, optrun ) ;  
		
				c = 0 ;
				for(k=0; k< s->h[window]->npar ; k++) 
				{
					if( s->p[window]->parinfer[k] ) { s->h[ window ]->theta[k] = opt[ thread_num ]->opt_theta[c] ;  c++ ; }
				}
				
				//recompute the Hessian here if including offsets 
				
				if( monopole )
				{
					//do not include the monopoles when building grid
					for( k=9; k<14; k++ ) s->p[window]->parinfer[k] = 0 ;
					s->p[window]->n_parinfer -= 5 ;
					opt[ thread_num ]->npar -= 5 ;
					
					par[ thread_num ]->s = s;
					par[ thread_num ]->thread_num = thread_num;
					
					gsl_matrix_free( s->h[window]->hessian ) ;
					
					s->h[window]->hessian = optimizer_get_hessian( (void *) par[ thread_num ], opt[ thread_num ], -1. ) ;
				
				}else{
				
					gsl_matrix_memcpy( s->h[ window ]->hessian , opt[ thread_num ]->hessian ) ;
				
				}
				
				optimizer_optpar_destroy( opt[ thread_num ] ) ;
			
				FILE * outf = fopen( "../WMAP_result/THETA_out.txt",  "a") ; 
				fprintf( outf , "%lf\t" , update_sync_link( s->h[ window ]->theta[0] ) );
				fprintf( outf , "%lf\t" , update_gdust_link( s->h[ window ]->theta[1] ) );
				for( k=2; k<s->h[ window ]->npar; k++ ) fprintf( outf , "%lf\t" , s->h[ window ]->theta[k] ); 
				fprintf( outf , "\n" ); 
				fclose( outf );	
			
				//construct the grid here...
			
				grid[ thread_num ] = grid_construct( s, thread_num ) ;
				
				log_marg_1 += grid[ thread_num ]->log_marginal_likelihood ; 
				
				printf("\n Patch %d : Window %d number of grid pts %d", base, window, grid[ thread_num ]->npts ) ;
			
			
				// **** NEED TO FIX THE POSTERIOR FUNCTIONALS FOR NEW SET UP ****
				grid_compute_posterior_functionals( grid[ thread_num ], s->p[ window ], s->h[window], s->b->mask[ window ] , covariance, residual ) ;
			
				results_write_block_to_results( window, results, s->p, s->b, s->h[ window ], FALSE ) ;
			
				printf("\n Patch %d : Window  %d completed ",base, window);
			
				grid_destroy( grid[ thread_num ] ) ; 
		
				}
				
			}
			
			free( s->p );
		}
		
		
		if( results_write_patch_result_to_file( base, results, s->b, rfile) ) printf("\n Problem  writing result to file" ) ;
		
		super_destroy( s ) ;
		results_destroy( results );
		
		log_marg_like += log_marg_1; 
		
		printf( "\n ********* finished base patch %d ************ ML: %lf", base, log_marg_like );
		
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
	free( model );
	for( k=0; k<num_thread; k++ ) free( par[k] );
	free( par );
	
	return(1);
}
