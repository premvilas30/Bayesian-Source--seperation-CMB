//test program for super

#include "defs.h"
#include "super.h"
#include "model.h"
#include "optimizer.h"

int main()
{
	int n = 16, N = n*n, k; 
	int *mask = calloc( N , sizeof(int) );
	for( k=0; k<N; k++ ) mask[k] = TRUE;
	
	double *obs_freq = calloc( 6, sizeof(double) ), *nullptr=NULL;
	double *obs_precision = calloc( 6, sizeof(double) );
	double *src_precision  =  calloc(  5, sizeof(double) );
	double *spec = calloc( 2, sizeof(double) );
	
	obs_freq[0] = 30.; obs_freq[1] = 44.; obs_freq[2] = 70. ;
	obs_freq[3] = 100.; obs_freq[4] = 143.;
	
	obs_precision[0] = 793.65; obs_precision[1] = 833.33;
	obs_precision[2] = 884.95; obs_precision[3] = 3571.43;
	obs_precision[4] = 5555.56;
	
	src_precision[0] = 100.00; src_precision[1] = 100.00;
	src_precision[2] = 100.00; src_precision[3] = 100.00;
	
	spec[0] = -2.85 ;
	spec[1] = 1.4 ;
	
	char **file = calloc( 10, sizeof( char *));
	file[0] = "../TEST_data/Y_1_test_4source_256.txt"; 
	file[1] = "../TEST_data/Y_2_test_4source_256.txt";
	file[2] = "../TEST_data/Y_3_test_4source_256.txt"; 
	file[3] = "../TEST_data/Y_4_test_4source_256.txt";
	file[4] = "../TEST_data/Y_5_test_4source_256.txt";
	
	file[5] = "../TEST_data/H_test_256.txt"; 
	file[6] = "../TEST_data/H_test_256.txt";
	file[7] = "../TEST_data/H_test_256.txt"; 
	file[8] = "../TEST_data/H_test_256.txt";
	file[9] = "../TEST_data/H_test_256.txt";
	
	struct super *s = super_create( 0, n, n, n, n, 5, 4, mask, FALSE) ;
	
	super_initialize( s, file, obs_freq, obs_precision, src_precision, spec, 1,  nullptr, nullptr); 
	
	/*update_patch( s->p[0], s->h, FALSE ) ;
	for( k=1; k<s->b->n_block+1; k++) 
	{
		update_patch( s->p[k], s->h, TRUE ) ;
		//cholmod_print_dense( s->p[k]->mu, "mu",  s->p[k]->chol_comm ) ;
	}
	patch_combine( s->p, s->b );
	
	FILE *out;
	
	double *x;
	int ntot;
	
	//non zero values in B
	out = fopen("../TEST_data/out_B_test_256.txt", "w");
	x = (double *) s->p[1]->B->x ;
	ntot = s->p[1]->B->nzmax ;
	for( k=0; k<ntot ; k++ ) fprintf( out, "%.10f\n", x[ k ]) ;
	fclose( out );
	
	//non zero values in C
	out = fopen("../TEST_data/out_C_test_256.txt", "w");
	x = (double *) s->p[1]->C->x ;
	ntot = s->p[1]->C->nzmax ;
	for( k=0; k<ntot ; k++ ) fprintf( out, "%.10f\n", x[ k ]) ;
	fclose( out );
	
	//non zero values in Q
	out = fopen("../TEST_data/out_Q_test_256.txt", "w");
	x = (double *) s->p[1]->Q->x ;
	ntot = (int) s->p[1]->Q->nzmax ;
	for( k=0; k<ntot ; k++ ) fprintf( out, "%.10f\n", x[ k ]) ;
	fclose( out );
	
	//non zero values in Q__
	out = fopen("../TEST_data/out_Q___test_256.txt", "w");
	cholmod_sparse *sp = cholmod_copy( s->p[1]->Q__ , -1, 1, s->p[1]->chol_comm);
	x = (double *) sp->x ;
	ntot = (int) sp->nzmax ;
	//cholmod_print_sparse( sp, "sp", s->p[1]->chol_comm );
	for( k=0; k<ntot ; k++ ) fprintf( out, "%.10f\n", x[ k ]) ;
	cholmod_free_sparse(  &sp, s->p[1]->chol_comm );
	fclose( out );
	
	//non zero values in mu
	out = fopen("../TEST_data/out_mu_test_256.txt", "w");
	x = (double *) s->p[0]->mu->x ;
	ntot = s->p[0]->mu->nrow ;
	for( k=0; k<ntot ; k++ ) fprintf( out, "%.10f\n", x[ k ]) ;
	fclose( out );	*/
	
	if( s->individual ) s->p[0] = s->p[1] ;
	
	double llike = model_log_likelihood( s->p[0]->mu , s->p[0] );
	
	double lprior = model_log_igmrf_prior( s->p[0]->mu , s->p[0] , s->h );
	
	double lpriorhpar = model_log_prior_hyperpar(s->h);
	
	double lmarginalhp = model_log_joint_marginal_hyperpar( s->p, s->h, 1 );
	
	
	out = fopen("../TEST_data/out_likecomp_test_256.txt", "w");
	fprintf( out, "%.20f\n%.20f\n%.20f\n%.20f",llike,lprior,lpriorhpar, lmarginalhp ) ;
	fclose( out );	
	
	
	return( 1 );
	
	printf("\n ... Entering optimization ... ") ; 
	
	struct optpar *opt = optimizer( s ) ; 
	
	out = fopen("../TEST_data/out_optim_comp_256.txt", "w");
	
	for( k=0; k<opt->npar; k++ )
	{
		fprintf( out, "%lf\t%lf\n", opt->opt_theta[k],opt->grad_opt_theta[k]);
	}
	
	fclose(out);
	
	optimizer_optpar_destroy( opt ) ; 
	
	super_destroy( s );
	
	free( mask );
	free( obs_freq );
	free( spec );
	free( obs_precision );
	free( src_precision );
	free( file );
	
	return(1);
}
