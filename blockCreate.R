
blockCreate <- function( nside, nrowblock, ncolblock, nmapin, data, obsprec, hitrate )
{
	J <- nside^2

	nblock <- nrowblock * ncolblock
	rblock <- nside / nrowblock
	cblock <- nside / ncolblock

	Jblock <-  rblock  *  cblock 

	blocks <- list()

	for( b in 1:nblock )
	{
		blocks[[ b ]] <- list( )
	
		#get out the indexes corresponding to each block (row major)

		i <- 1 : cblock
	
		off <- floor( ( b  - 1 ) / ncolblock ) * ( J / nrowblock ) + ( b - 1 ) * cblock
	
		j <- off + 1 : cblock
	
		a <- NULL
	
		for( rowstart in 0:( rblock  - 1 ) ) a <- c( a, rowstart * nside + j )
		
		blocks[[b]]$J <- Jblock
	
		blocks[[ b ]]$idx <- a
		
		blocks[[b]]$y <- numeric( Jblock * nmapin )
		
		C <- numeric( Jblock * nmapin )
		
		for( m in 1:nmapin )
		{
			blocks[[b]]$y[ a + (m-1)*Jblock ] <- data[ m , a ]
			C[ a + (m-1)*Jblock ] <- obsprec[ m ] * hitrate
		}
	}

}
