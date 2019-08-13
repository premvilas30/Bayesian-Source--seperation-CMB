//update.c

#include "update.h"

double update_sync_link( double x )
{
	 return( ( SYNC_UPPER * exp(x) + SYNC_LOWER )/( 1.+ exp(x) ) ); //transform back to interval [SYNC_LOWER, SYNC_UPPER]
}

double update_gdust_link( double x )
{
	return( ( DUST_UPPER * exp(x) + DUST_LOWER )/( 1.+ exp(x) ) ); //transform back to interval [DUST_LOWER, DUST_UPPER]
}

double update_specpar_link( double x )
{
	return( ( GEN_SPECPAR_UPPER * exp(x) + GEN_SPECPAR_LOWER )/( 1.+ exp(x) ) );
}

double update_precision_link( double x )
{
	return( exp( x ) ) ; //transform back to interval [0,\infty]
}


//temperature conversion factors to be initialized
void update_temperature_conversion_factors( struct patch *p, struct hyperpar *h )
{

	int k;
	double a; 
	
	for( k=0; k<p->nin_map; k++ )
	{
		a  = pow( 10., log10( h->nu[k] ) + 9. + log10PLANCK - log10BOLTZMANN - log10T0 );
		p->conv_ant_to_therm[k] = exp( -  ( 2.*log( a ) + a - 2.* log( gsl_expm1(a) ) ) ) ;
	}

}

// functions for various types of loadings in the matrix A

double update_cmb_tdymc( double nu, double nu_ref, double beta )
{
	return( 1. );
}

double update_pow_law( double nu, double nu_ref, double beta )
{
	return( pow( nu/nu_ref, beta ) ) ;
}

double update_mod_bb_spec( double nu, double nu_ref, double beta )
{
	double a, b;
	a = pow( 10., log10( nu ) + 9. + log10PLANCK - log10BOLTZMANN - log10T1 ) ;
	b = pow( 10., log10( nu_ref ) + 9. + log10PLANCK - log10BOLTZMANN - log10T1  ) ;
	return( exp( log( gsl_expm1(b) ) - log( gsl_expm1(a) ) ) );
}

void update_A( struct patch *p, struct hyperpar *h )
{
	
	int k;
	
	for( k=0; k<p->nin_map; k++ ) p->A[k][0] = 1.;
	
	double beta = update_specpar_link( h->theta[0] );
	
	for( k=0; k<p->nin_map; k++ ) p->A[k][1] = update_pow_law( h->nu[k] , h->nu[0], beta );
	
	
}


void update_A_x( struct patch *p, struct hyperpar *h )
{
	int k, c = 0;
	double *w;//, *ant_to_therm;
	
	int ref_freq = 3;
	
	//must correct for antenna to thermodynamic conversion factor
	
	//modify the structure of this update slightly 
	// include the special case when all sources are inferred and uniqueness contraint needed
	
	w = calloc( p->nin_map + 1, sizeof(double) );

	
	if( !p->template_source[0] )
	{
	
		for( k=0; k<p->nin_map; k++ ) p->A[k][0] = 1.; //frequency independent spectrum in Thermodyanmic temperature
		
		c++;
	
	}
	
	
	if( !p->template_source[1] )
	{
	
		//Synchrotron
		
		double s = update_sync_link( h->theta[0] ), nu_ref = 23.; //Haslam map
		
		for( k=0; k<p->nin_map; k++ ) 
		{
			if( p->prior_mu_source[1] )
			{
				p->A[k][c] = pow( h->nu[k] / p->mu_prior_ref_freq[1],  s ) * p->conv_ant_to_therm[k] ;
			}
			else
			{
				p->A[k][c] = pow( h->nu[k] / h->nu[ref_freq],  s ) * p->conv_ant_to_therm[k] ;
			}
		}
		
		c++;
	
	}
	
	if( !p->template_source[2] )
	{
	
		//Galactic dust: a modified black body spectrum multiplied by conversion factor

		
		double d = update_gdust_link( h->theta[1] );
		ref_freq = 4;
	
		for( k=0; k<p->nin_map; k++) 
		{
			if( p->prior_mu_source[2] ) 
				p->A[k][c] = pow( h->nu[k] / p->mu_prior_ref_freq[2] , 1 + d ) * p->conv_ant_to_therm[k] ;
			else
				p->A[k][c] = pow( h->nu[k] / h->nu[ref_freq], 1 + d ) * p->conv_ant_to_therm[k] ;
		}
		
		for( k=0; k<p->nin_map; k++) 
		{
			w[k] = pow( 10., log10( h->nu[k] ) + 9. + log10PLANCK - log10BOLTZMANN - log10T1 );
		}
	

		if( p->prior_mu_source[2] )	
			w[ p->nin_map ] = pow( 10., log10( p->mu_prior_ref_freq[2] ) + 9. + log10PLANCK - log10BOLTZMANN - log10T1 );
	
		for( k=0; k<p->nin_map; k++) 
		{
			if( p->prior_mu_source[2] ) //p->template_source[1] ) //p->ninferred_map < p->nout_map )
				p->A[k][c] = p->A[k][c] * exp( log( gsl_expm1( w[ p->nin_map ] ) ) - log( gsl_expm1(w[k]) ) );
			else
				p->A[k][c] = p->A[k][c] * exp( log( gsl_expm1(w[ref_freq]) ) - log( gsl_expm1(w[k]) ) );
		} 
		
		c++;
	
	}
	
	
	if( !p->template_source[3] )
	{
	
		//free-free emission
	
		//nu_ref = 23.;
		//ref_freq = 0;
	
		for( k=0; k<p->nin_map; k++) 
		{
			if( p->prior_mu_source[3] ) // p->ninferred_map < p->nout_map )
				p->A[k][c] = pow( h->nu[k] / p->mu_prior_ref_freq[3] , -2.15 ) * p->conv_ant_to_therm[k] ;
			else
				p->A[k][c] = pow( h->nu[k] / h->nu[ref_freq] , -2.15 ) * p->conv_ant_to_therm[k] ;
		}
		
		c++;
	
	}	

	free(w);
	
	return;	
}

void update_A_simulated( struct patch *p, struct hyperpar *h )
{

	int k, c = 0;
	double *w, *v;//, *ant_to_therm;
	
	double nu_ref;
	
	//modify the structure of this update slightly 
	// include the special case when all sources are inferred and uniqueness contraint needed
	
	w = calloc( p->nin_map + 1, sizeof(double) );
	v = calloc( p->nin_map, sizeof(double) );
	//ant_to_therm = calloc( p->nin_map , sizeof(double) );
	
	double arg;
	
	if( !p->template_source[0] ) 
	{
		for( k=0; k<p->nin_map; k++ ) p->A[k][c] = 1.; //frequency independent spectrum in Thermodyanmic temperature
		
		c++;
	}	

	if( !p->template_source[1] )
	{
		//Synchrotron
		
		v[0] = 2.; v[1] = 1.; v[2] = 3.; v[3] = 4; v[4] = 5;
		
		double s = update_sync_link( h->theta[0] );
		
		for( k=0; k<p->nin_map; k++ ) 
		{
			p->A[k][c] = pow( v[k] /*/ v[1]*/ ,  s );
		}
		
		c++;
	}
	
	if( !p->template_source[2] )
	{
		
		v[0] = 5.; v[1] = 4.; v[2] = 3.; v[3] = 2.; v[4] = 1.; 

		double d = update_gdust_link( h->theta[1] );
	
		for( k=0; k<p->nin_map; k++) 
		{
			p->A[k][c] =  pow( v[k]/* / v[4]*/, 1 + d ) ; 
		}
		
		c++;
	
	}
	
	if( !p->template_source[3] )
	{
		//free-free emission
		
		v[0] = 2.; v[1] = 2.; v[2] = 1.; v[3] = 2.; v[4] = 2.;
	
		for( k=0; k<p->nin_map; k++) 
		{
			p->A[k][c] = pow( v[k] /*/ v[2]*/ , -2. ) ;
		}
		
		c++;
	}	
	
	/*int l;
	for( k=0; k<p->nin_map; k++ )
	{
		printf("\n");
		for( l=0; l<c; l++ ) printf("%lf, ", p->F[k][l]);
	}
	
	printf("\n ====== \n");*/


	free(w);
	free(v);
	
	return;

}



//update_B
void update_B_pixel_specific( struct patch *p, struct hyperpar *h )
{
	//update the values in p->A first
	
	update_A( p , h );
	
	int 	k, l, c, nnode = p->graph->number_nodes,
			*colptr = (int *)p->B->p ;
	double *x =  (double *)p->B->x, nu_ref;
	
	for( k=0; k<p->ninferred_map; k++ )
	{
		for( c=k*nnode; c<(k+1)*nnode; c++ )
		{
			for( l=0; l<p->nin_map; l++ ) x[ colptr[ c ] + l ] = p->A[l][k];
		}
		
		//special case for synchrotron if the parameter is to vary as a field
		if( p->specidx[0] &&  k == 1 )
		{
			//Synchrotron
		
			nu_ref = 23.; //Haslam map
		
			/*for( k=0; k<p->nin_map; k++ ) 
			{
				if( 0 )//p->template_source[1] )//p->ninferred_map < p->nout_map )
					p->A[k][c] = pow( h->nu[k] / nu_ref,  s ) ;
				else
					p->A[k][c] = pow( h->nu[k] / h->nu[ref_freq],  s ) ;
			}*/
			
			for( c=k*nnode; c<(k+1)*nnode; c++ )
			{
				for( l=0; l<p->nin_map; l++ ) x[ colptr[ c ] + l ] = pow( h->nu[l] / nu_ref, ((double*)p->synch_index->x)[ c - k*nnode ] ) ;
			}
		
		}
		
	}
	
	return;
}



//next thing to do is to write update_F

void update_F( struct patch *p, struct hyperpar *h )
{

	int k, c = 0;
	double *w;//, *ant_to_therm;
	
	double nu_ref;
	
	//modify the structure of this update slightly 
	// include the special case when all sources are inferred and uniqueness contraint needed
	
	w = calloc( p->nin_map + 1, sizeof(double) );
	//ant_to_therm = calloc( p->nin_map , sizeof(double) );
	
	double arg;
	/*for( k=0; k<p->nin_map; k++ )
	{
		//arg = pow( 10., log10( h->nu[k] ) + 9. + log10PLANCK - log10BOLTZMANN - log10T0 );
		//ant_to_therm[k] = exp( - ( 2.*log( arg ) + arg - 2.* log( gsl_expm1(arg) ) ) ) ; 
	}*/
	
	if( p->template_source[0] ) 
	{
		for( k=0; k<p->nin_map; k++ ) p->F[k][c] = 1.; //frequency independent spectrum in Thermodyanmic temperature
		
		c++;
	}	

	if( p->template_source[1] )
	{
		//Synchrotron
		
		//pow( h->nu[k] / p->mu_prior_ref_freq[1],  s ) * ant_to_therm[k] ;
		
		//in this case the "precision" hyperpar acts as an amplitude multiplier
		
		double s = update_sync_link( h->theta[0] );
		
		for( k=0; k<p->nin_map; k++ ) 
		{
			p->F[k][c] = exp( h->theta[6] ) * pow( h->nu[k] / p->template_ref_freq[1] ,  s ) * p->conv_ant_to_therm[k] ;
		}
		
		c++;
	}
	
	if( p->template_source[2] )
	{
		//Galactic dust (using the FDS template)

		double d = update_gdust_link( h->theta[1] );
	
		for( k=0; k<p->nin_map; k++) 
		{
			p->F[k][c] =  pow( h->nu[k] / p->template_ref_freq[2], 1 + d ) * p->conv_ant_to_therm[k] ; //part with emissivity index
		}
		
		for( k=0; k<p->nin_map; k++) 
		{
			w[k] = pow( 10., log10( h->nu[k] ) + 9. + log10PLANCK - log10BOLTZMANN - log10T1 );
		}
		
		w[ p->nin_map ] = pow( 10., log10( p->template_ref_freq[2] ) + 9. + log10PLANCK - log10BOLTZMANN - log10T1 );
	
		for( k=0; k<p->nin_map; k++) 
		{
				p->F[k][c] = exp( h->theta[7] ) * p->F[k][c] * exp( log( gsl_expm1(w[p->nin_map]) ) - log( gsl_expm1(w[k]) ) );
		} 
		
		c++;
	
	}
	
	if( p->template_source[3] )
	{
		//free-free emission
	
		for( k=0; k<p->nin_map; k++) 
		{
			p->F[k][c] = exp( h->theta[8] ) * pow( h->nu[k] / p->template_ref_freq[3] , -2.14 ) * p->conv_ant_to_therm[k] ;
		}
		
		c++;
	}	


	free(w);
	//free(ant_to_therm);
	/*int l;
	printf("\n******");
	for( k =0 ; k<p->nin_map; k++ )
	{
		for( l=0; l<c; l++ ) printf("\t %.10f", p->F[k][l]);
		printf("\n");
	}
	printf("\n******");*/
	
	return;

}

void update_F_simulated( struct patch *p, struct hyperpar *h )
{

	int k, c = 0;
	double *w, *v;//, *ant_to_therm;
	
	double nu_ref;
	
	//modify the structure of this update slightly 
	// include the special case when all sources are inferred and uniqueness contraint needed
	
	w = calloc( p->nin_map + 1, sizeof(double) );
	v = calloc( p->nin_map, sizeof(double) );
	//ant_to_therm = calloc( p->nin_map , sizeof(double) );
	
	double arg;
	
	if( p->template_source[0] ) 
	{
		for( k=0; k<p->nin_map; k++ ) p->F[k][c] = 1.; //frequency independent spectrum in Thermodyanmic temperature
		
		c++;
	}	

	if( p->template_source[1] )
	{
		//Synchrotron
		
		v[0] = 2.; v[1] = 1.; v[2] = 3.; v[3] = 4; v[4] = 5;
		
		double s = update_sync_link( h->theta[0] );
		
		for( k=0; k<p->nin_map; k++ ) 
		{
			p->F[k][c] = pow( v[k] /*/ v[1]*/ ,  s );
		}
		
		c++;
	}
	
	if( p->template_source[2] )
	{
		
		v[0] = 5.; v[1] = 4.; v[2] = 3.; v[3] = 2.; v[4] = 1.; 

		double d = update_gdust_link( h->theta[1] );
	
		for( k=0; k<p->nin_map; k++) 
		{
			p->F[k][c] =  pow( v[k]/* / v[4]*/, 1 + d ) ; 
		}
		
		c++;
	
	}
	
	if( p->template_source[3] )
	{
		//free-free emission
		
		v[0] = 2.; v[1] = 2.; v[2] = 1.; v[3] = 2.; v[4] = 2.;
	
		for( k=0; k<p->nin_map; k++) 
		{
			p->F[k][c] = pow( v[k] /*/ v[2]*/ , -2. ) ;
		}
		
		c++;
	}	
	
	/*int l;
	for( k=0; k<p->nin_map; k++ )
	{
		printf("\n");
		for( l=0; l<c; l++ ) printf("%lf, ", p->F[k][l]);
	}
	
	printf("\n ====== \n");*/


	free(w);
	free(v);
	
	return;

}



//update_B
void update_B( struct patch *p, struct hyperpar *h )
{
	//update the values in p->A first
	
	if( UPDATE_SIMULATED ) update_A_simulated( p, h ); else update_A( p , h );
	
	//cholmod_print_sparse( p->B,"B", p->chol_comm );
	
	int 	k, l, c, nnode = p->graph->number_nodes,
			*colptr = (int *)p->B->p ;
	double *x =  (double *)p->B->x;
	
	for( k=0; k<p->ninferred_map; k++ )
	{
		for( c=k*nnode; c<(k+1)*nnode; c++ )
		{
			for( l=0; l<p->nin_map; l++ ) x[ colptr[ c ] + l ] = p->A[l][k];
		}
	}
	
	return;
}


//update_G
void update_G( struct patch *p, struct hyperpar *h )
{
	
	//update the values in p->F first
	
	//printf("\n ok to  here b4 F");
	//printf("\n ok to  here 0");
	
	if( UPDATE_SIMULATED ) update_F_simulated( p, h ); else update_F( p, h ); 
	
	//printf("\n ok to  here aftr F");
	//printf("\n ok to  here 1");
	
	int k, l, c, nnode = p->graph->number_nodes, 
		*colptr = (int *)p->G->p;
	double *x = (double *)p->G->x;
	
	for( k=0; k<p->nin_template; k++ )
	{
		//printf("\n ok to  here temp %d",k);
		for( c=k*nnode; c<(k+1)*nnode; c++ )
		{
			for( l=0; l<p->nin_map; l++ ) x[ colptr[c] + l ] = p->F[l][k];
		}
	}
	
	return;

}

//update_Q : precision matrix for IGMRF prior
void update_Q( struct patch *p, struct hyperpar *h )
{
	//use p->T and the values in h->theta
	
	int 	k, l, c, k_ = 0, nnode = p->graph->number_nodes,
			*colptr = (int *)p->T->p;
	double	*q = (double *)p->Q->x,
				*t = (double *)p->T->x;
	
	
	for( k=0; k<p->nout_map; k++ )
	{
		if( !p->template_source[k] )
		{
			for( c=k_*nnode; c<(k_+1)*nnode; c++ )
			{
				for( l=colptr[c]; l<colptr[c+1]; l++ ) 
				{
					q[ l ] = update_precision_link( h->theta[ 2 + k ] ) * t[ l ]; 
					//printf("\n Entry[%d,%d] = %.10f, theta = %lf", l, c, q[l], h->theta[2+k] );
				}
			}
			k_++;
		}
	}

}

//need to fix this 

void update_Q_using_mask( struct patch *p, struct hyperpar *h, int *mask )
{
	//use p->T and the values in h->theta
	
	int 	k, l, c, k_=0, nnode = p->graph->number_nodes,
			*colptr = (int *)p->T->p, node;
	double	*q = (double *)p->Q->x,
				*t = (double *)p->T->x;
	
	for( k=0; k<p->nout_map; k++ )
	{
	
		if( 0 )
		{
			for( c=k_*nnode; c<(k_+1)*nnode; c++ )
			{
				node = c - k*nnode ;
				for( l=colptr[c]; l<colptr[c+1]; l++ ) q[ l ] = update_precision_link( h->theta[ 2 + k ] ) * t[ l ] ;
			
				if( !mask[ node ] ) 
				{
					l = colptr[ c ] ; //access the diagonal entry
					q[ l ] = 1E10; //set to large precision to force to zero
				}
			
			}
			k_++;
		}
	}

}


//a new version of this needs to be written... 

void update_Q__( struct patch *p, struct  hyperpar *h )
{
	//this assumes that the values of p->B and p->Q have
	// already been updated
	
	int k, l, nd = p->nin_map * p->graph->number_nodes  ;
	
	double 	*a = calloc(2,sizeof(double)),
				*b = calloc(2,sizeof(double));
	a[0] = 1.;
	

	
	cholmod_sparse *t_B = cholmod_transpose( p->B, 1, p->chol_comm );
	
	cholmod_sparse *t_BC = cholmod_ssmult( t_B, p->C, 0, TRUE, FALSE, p->chol_comm );
	
	cholmod_sparse *t_BCB = cholmod_ssmult( t_BC, p->B, -1, TRUE, TRUE, p->chol_comm );

	 //compute t_BCy for the mode later also
		
	//include term from prior Q\mu 
		
	cholmod_sdmult( p->Q, 0, a, b, p->mu_prior, p->z, p->chol_comm );
		
		
	//add B^t C y
		
	b[0] = 1.;
	
	cholmod_sdmult( t_BC, 0, a, b, p->y, p->z , p->chol_comm );
	
		
	
	if( p->nin_template > 0 )
	{
		
		a[0] = 1.;
		b[0] = 0.;
		
		cholmod_sdmult( p->G, 0, a, b, p->mu_template, p->z1, p->chol_comm );
		
		//cholmod_print_sparse(p->G, "G", p->chol_comm );
	
		// this line here  needs to be modified for the new set up (**)
	
		a[0] = -1.;
		b[0] = 1 ;
		
		cholmod_sdmult( t_BC, 0, a, b, p->z1, p->z , p->chol_comm );
	
	}

	
	//add a line here for the monopoles
	
/*	for( k=0; k<p->nin_map; k++ )
	{
		for( l=0; l<p->graph->number_nodes; l++ )
		{
			((double*) p->z->x)[ k * p->graph->number_nodes + l ] -= h->theta[ 9 + k ] ; 
		}	
	}*/
	
	//free  what is not needed
	cholmod_free_sparse( &t_B, p->chol_comm );
	
	cholmod_free_sparse( &t_BC, p->chol_comm );

	
	//add the prior precision
	
	a[0] = 1.;
	b[0] = 1.; 
	
	cholmod_free_sparse( &p->Q__, p->chol_comm );
	
	cholmod_sparse *U;
	
	// what happens here if everything is templated?
	
	U = cholmod_add( t_BCB, p->Q, a, b, TRUE, TRUE, p->chol_comm );
	
	p->Q__ = cholmod_copy( U, -1, 1, p->chol_comm );
	
	cholmod_free_sparse( &U, p->chol_comm);
	
	cholmod_free_sparse( &t_BCB , p->chol_comm );
	
	free(a);
	free(b);
	
	return;
}

void update_L_initial( struct patch *p )
{

	p->L = cholmod_analyze( p->Q__, p->chol_comm);
	
	//cholmod_factorize( p->Q__, p->L, p->chol_comm );
	
	//cholmod_print_factor( p->L,"L",p->chol_comm );
	return;
}

void update_mu( struct patch *p )
{
	
	cholmod_factorize( p->Q__, p->L, p->chol_comm );
	
	//cholmod_print_factor( p->L,"L",p->chol_comm );
	
	//get determinant of p->Q__ and store
	int 	ncol = p->L->n, k,
			*colptr = (int *)p->L->p;
	
	double 	*x = (double *) p->L->x, 
				ldet = 0.;
	
	for( k=0; k<ncol; k++ ) ldet += log( x[ colptr[k] ] ); 
	
	p->log_det_Q__ = 2. * ldet;
	
	
	//get the mode of the full conditional
	
	cholmod_free_dense( &p->mu, p->chol_comm );
	
	p->mu = cholmod_solve( CHOLMOD_A, p->L, p->z, p->chol_comm );
	
	/*double norm = 0.;
	double *vv = (double *) p->z->x ;
	
	for( k=0; k<p->z->nrow; k++ ) norm += vv[k] * vv[k] ;
	
	printf("\n The norm of mu is %lf ", sqrt(norm) );*/
	
	return;
}

cholmod_dense *update_linear_predictor( struct patch *p, struct hyperpar *h )
{
	int k, l, nd = p->nin_map  * p->graph->number_nodes ;
	double 	*a = calloc(2,sizeof(double)),
				*b = calloc(2,sizeof(double));
	a[0] = 1.;
	
	cholmod_dense *x = cholmod_zeros( nd, 1, CHOLMOD_REAL, p->chol_comm ) ;
	
	if( p->ninferred_map > 0 )
	{
	
		cholmod_sdmult( p->B, 0, a, b, p->mu, x , p->chol_comm ) ; 
	
	}
	//printf("\n*****1**");
	//for( k=0; k<5; k++ ) printf("%lf\n", ((double*)x->x)[k] );

	
	
	if( p->nin_template > 0 )
	{
		//add the Gs term
		
		//multiply G by mu_template and add this to x
		
		b[0] = 1.;
		
		cholmod_sdmult( p->G, 0, a, b, p->mu_template, x, p->chol_comm ) ; 
	} 
	
	
	for( k=0; k<p->nin_map; k++ )
	{
		for( l=0; l<p->graph->number_nodes; l++ )
		{
			((double*) x->x)[ k * p->graph->number_nodes + l ] += h->theta[ 9 + k ] ; 
		}	
	}
	
	
	//printf("\n***2****");
	//for( k=0; k<5; k++ ) printf("%lf\n", ((double*)x->x)[k] );
	
	free( a ) ; 
	free( b ) ;
	
	return( x ) ;
}

void update_covariance( struct patch *p )
{
	
	//printf("\n Inside update_covariance");
	//printf("\n Inside update_covariance");
	
	int nl = p->nout_map * p->graph->number_nodes, k ;
	
	//form the dense identity matrix
	
	cholmod_dense *Id = cholmod_allocate_dense( nl, nl, nl, CHOLMOD_REAL,  p->chol_comm );
	
	for( k=0; k<nl; k++ ) ( (double *) Id->x )[ k * nl + k ] = 1. ;
	
	cholmod_free_dense( &p->Sig, p->chol_comm ) ;
	
	//compute the inverse of the precision matrix	
	
	p->Sig = cholmod_solve( CHOLMOD_A, p->L, Id, p->chol_comm ) ;
	
	cholmod_free_dense( &Id, p->chol_comm ) ;
	
	//printf("\n Done update_covariance");
	//printf("\n Done update_covariance");
	
	return ;
}

void update_patch( struct patch *p, struct hyperpar *h, int all, int *mask )
{
	
	int k;
	for( k=0; k<h->npar; k++ )
	{
		//printf("\n h->theta[%d] = %.20f ...", k, h->theta[k]);
	}
	
	//update A, B and Q if there are maps being inferred
	if( p->ninferred_map > 0 ) 
	{
		update_B( p, h ); 
		if( !p->masked_pix_include ) update_Q( p, h );
	}
	
	// if all maps are assumed "known" templates
	// just update G
	if( p->nin_template > 0 ) update_G( p, h );

	
	if( all && p->ninferred_map > 0 )
	{
	
		update_Q__( p, h );
	
		update_mu( p );
		
	}
	
	return;
}

void update_initial( struct block *b, struct patch **p, struct hyperpar **h  )
{
	int k;
	
	update_temperature_conversion_factors( p[0], h[0] ) ;
	
	update_patch( p[0], h[0], FALSE,  b->mask[0] ) ;
	
	for( k=1; k< b->n_block+1; k++ ) 
	{
		
		//printf("\n block %d done", k ) ;
		
		if( !UPDATE_SIMULATED ) update_temperature_conversion_factors( p[k], h[k] ) ;
		
		update_patch( p[k], h[k], FALSE, b->mask[k] ) ;	
		
		if( p[k]->ninferred_map > 0 )
		{
			update_Q__( p[k],  h[k] ) ;
	
			//initialize the Cholesky factor
			update_L_initial( p[k] );
		}
	}
	return;
}

void update_simulate_data( struct patch *p, struct hyperpar *h, int *mask, int seed )
{
	
	update_patch( p, h, FALSE, mask );
				
	cholmod_dense *x = update_linear_predictor( p , h );
	
	int nd = p->nin_map * p->graph->number_nodes, k ;
	
	double *y = (double *) p->y->x, *yx = (double *) x->x ;
	
	//set up a RNG with the given seed
	
	const gsl_rng_type *T;
	
 	gsl_rng *r;
  	
  	gsl_rng_env_setup();
  	
  	T = gsl_rng_default;
  	
  	r = gsl_rng_alloc(T);
  
  	gsl_rng_set(r,(int)time(NULL));
  	
  	double sigma ;
	
	for( k=0; k<nd; k++ ) 
	{
		sigma = 1./sqrt( ((double *) p->C->x)[k] );
		y[k] = yx[k] + gsl_ran_gaussian( r , sigma );
		//printf("\n signal = %lf with error %lf ", yx[k], y[k] );
	}
	
	cholmod_free_dense( &x, p->chol_comm );
	
	gsl_rng_free( r );
		
	return;		

}



