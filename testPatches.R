#testing script to compare outputs from C and  R

library(Matrix)

source( "AB.R" )
source( "simIGMRF.R" )

ns <- 1
nf <- 1
nside <- 16
J <- nside ^ 2

phi <- 100

sync <- -2.8
dust <- 1.4

B <- getB( J, ns, sync, dust )

S <- simIGMRF( J, ns, phi )


