#test patches blocks

library(Matrix)

source( "ABC.R" )
source( "simIGMRF.R" )

ns <- 2
nf <- 2
nside <- 16
J <- nside ^ 2

y1 <- ( file = paste( "../../TEST_data/Y1_test_" ,J,".txt",sep="") )
y2 <- ( file = paste( "../../TEST_data/Y2_test_" ,J,".txt",sep="") )
Y <- c(y1,y2)

phi <- c( 500, 600)
tau <- c( 793.65 , 833.33 )#, 884.95, 3571.43, 5555.56, 5555.56 )

sync.u <- -2.3
sync.l <- -3.
dust.u <- 2.
dust.l <- 1.

sync <- -2.85
dust <- 1.4

theta <- numeric(4)
theta[1] <- -log( (sync.u - sync)/(sync - sync.l) )
theta[2] <- -log((dust.u - dust)/(dust - dust.l))
theta[3:4] <- log(phi)

hit <- rep( 1,  J * nf)

B <- getB( J, ns, sync, dust )

C <- getC( J, tau, hit )

Q <- bdiag( phi[1] * precision.IGMRF( J ), phi[2] * precision.IGMRF( J ) )

nrowblock <- 2
ncolblock <- 2
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
	
	blocks[[ b ]]$idx <- a
}

#Qstr <- Q + t(B) %*% C %*% B

#L <- chol( Qstr )

#log.det.Qstr <- 2 * sum( log( diag( L ) ) )

#Qinv <-  chol2inv( L )

#mustar <- Qinv %*%  t(B) %*% C %*% Y
