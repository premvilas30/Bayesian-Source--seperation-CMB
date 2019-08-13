#evaluate the log posterior of the data at mustar

llike <- - .5 *  nf* J * log( 2*pi) +  0.5 * sum( log( diag(C) ) ) - 0.5* t( Y - B %*% mustar) %*% C %*% ( Y - B %*% mustar ) 

lprior <- -.5 * ns  *  (J-1) * log( 2 * pi ) + .5  * (J-1) * sum( log( phi) ) - .5 * t( mustar ) %*% Q %*% mustar

lpriorhpar1 <- log( sync.u*exp(theta[1]) - sync.l ) + theta[1] - 2 * log( 1 + exp(theta[1]) )

lpriorhpar2 <- log( dust.u*exp(theta[2]) - dust.l ) + theta[2] - 2 * log( 1 + exp(theta[2]) )

lpriorhpar3 <- dgamma( exp( theta[3:4] ), shape=1, rate = .000001 , log=T ) + theta[3:4]

lpriorhpar <- lpriorhpar1 + lpriorhpar2 + sum( lpriorhpar3 )

lcondmode <- - .5 * ns * J * log( 2*pi ) + .5 * 2. * sum( log( diag(L) ) )

u <- c( llike@x , lprior@x, lpriorhpar, lcondmode )

#check residual difference in matrices

x <- scan("../../TEST_data/out_B_test_256.txt", quiet=T)
cat( "\n Residual for B = ",  all.equal( x, B@x ) )

x <- scan("../../TEST_data/out_C_test_256.txt", quiet=T)
cat( "\n Residual for C = ", all.equal( x, C@x ) )

Q <- Matrix(Q)
Q <- tril(Q)
x <- scan("../../TEST_data/out_Q_test_256.txt",  quiet=T)
cat( "\n Residual for Q = ", all.equal( x, Q@x ) )

Qstr <-   tril(Qstr)
x <- scan("../../TEST_data/out_Q___test_256.txt", quiet=T)
cat( "\n Residual for Q__ = ", all.equal( x, Qstr@x ) )

x <- scan("../../TEST_data/out_mu_test_256.txt", quiet=T)
cat( "\n Residual for mu = ", all.equal( x, mustar@x ) )

x <- scan("../../TEST_data/out_likecomp_test_256.txt", quiet=T)
cat( "\n Residual for ldens = ", sum( ( x - u )^2 ) )
