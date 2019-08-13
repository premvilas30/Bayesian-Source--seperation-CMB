// code for carrying out MCMC inference for spectral parameters

#include "mcmc.h"

struct mcmc *mcmc_create( int npar, int it, int burn, int nnodes, int nout_map )
{
	struct mcmc *mcmc = ( struct mcmc *)malloc( sizeof(struct mcmc) ) ;
	mcmc->it = it;
	mcmc->burn = burn;
	mcmc->npar = npar;
	mcmc->theta = calloc( npar, sizeof(double) );
	mcmc->theta_prop = calloc( npar, sizeof(double) );
	mcmc->s = calloc( nnode * nout_map, sizeof(double) );
	mcmc->log_posterior = -DBL_MAX;
	mcmc->log_posterior_prop = -DBL_MAX;
	return( mcmc );
}

void mcmc_destroy( struct mcmc *mcmc )
{
	free( mcmc->theta );
	free( mcmc->s );
}
