
# a little program to look at the Jeffreys priors ( up to multiplicative constants )

# rough approximation to channel noise variances

v <- c( 310 * 9.323038*1E-9, 640 * 1.324985*1E-8, 2178 * 1.008426*1E-8, 650 * 3.785305*1E-9, 911 * 6.69351*1E-10, 846 * 1.386299*1E-9, 979 * 1.376246*1E-8 ) *1E6

tau <- 1/v

nu <- c( 30, 44,  70, 100, 143, 217, 353 )

# look at the synchrotron spectral prior first

syncprior <- function( beta )
{
	#param on range -3.3 to -2.7
	if( length(beta) > 1 )
	{
		len <- length(beta)
		st <- numeric( len )
		for( k in 1:len ) st[k] <- syncprior( beta[k] )
		return( st )
	}else{
		u <- tau * ( nu /nu[1] )^(2*beta) * ( log(  nu/nu[1] ) )^2
		return( sqrt( sum(u) ) )
	}
}

log10planck <- -33.1787441114855639057
log10boltz <- -22.859916308396282858 

w <- log10( nu ) + 9 + log10planck - log10boltz - log10( 18.1 )
w <-10**w

dustprior <- function( beta )
{
	#param on range 1. to 2.
	if( length(beta) > 1 )
	{
		len <- length(beta)
		st <- numeric(len)
		for( k in 1:len ) st[k] <- dustprior( beta[k] )
		return( st )
	}else{
		c <- exp( log( expm1( w[4] ) ) - log(  expm1( w )  ) ) * log( nu/nu[4] )  
		u <- tau * c^2  * ( nu / nu[4] )^(2+beta)
		return( sqrt( sum(u) ) )
		
	}
}

nu <- c(30,44,70,100)
sigma.nu  <- 1

b.prior <- function( x )
{
	if( length(x) > 1)
	{
		v <- numeric( length(x) )
		for( k in 1:length(x) )  v[k] <- b.prior(x[k])
		return(v)
	}else{
		u <- ( ( nu/nu[1] )^x * log( nu  / nu[1] )/sigma.nu )^2
		return( sqrt( sum(u) ) )
	}
}









		
