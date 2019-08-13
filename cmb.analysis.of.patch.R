
library( Matrix )

cmb.analysis.of.patch <- function( experiment="wmap7yr", Nside.block=32, nblock=1, row=1, col=1, covariance=FALSE, residual=TRUE, par.infer=rep(1,6)  ) 
{
	maxGridPts <-10000

	Experiments <- c( "wmap7yr", "wmap9yr", "planck" )
	ExperimentsNside <- c( 512, 512, 1024 )
	ExperimentsMapIn <- c( 5, 5, 7 )
	ExperimentsMapOut <- c( 4, 4, 4 )
	
	whichdata <- which( Experiments == experiment )
	
	Nside <- ExperimentsNside[ whichdata ]
	MapIn <- ExperimentsMapIn[ whichdata ]
	MapOut <- ExperimentsMapOut[ whichdata ]
	
	NpixelBlock <- Nside.block^2
	
	nBlocks <- length( row ) * length( col )
	
	postExpectMu <- numeric( nBlocks * MapOut * NpixelBlock )
	residualVec <- numeric( nBlocks *  MapIn * NpixelBlock )
	usedGridPts <- numeric( nBlocks )
	
	GridPts <- numeric( nBlocks * maxGridPts * ( sum(par.infer) + 1 ) )
	
	
	w <- .C( "CMB_analysis_of_patch", 	as.integer( whichdata-1 ), 						as.integer( Nside.block ) 
													as.integer( nBlocks ), 								as.integer( row ),
													as.integer( col ),									as.integer( covariance ),
													as.integer( residual ),								as.integer( par.infer ),
													postExpectMu = as.double( postExpectMu ),		residual = as.double( residualVec ),
													as.integer( maxGridPts ),							usedGridPts = as.integer( usedGridPts ),
													GridPts  = 	as.double( GridPts ),				as.integer( 1 )
			)
	
	Ret <- list()
	Ret$postExpectMu <- list()
	if( residual ) Ret$residual <- list()
	Ret$GridPts <- list()
	
	for(  b in 1:nBlocks )
	{	
		x <- w$postExpectMu[ ( (b-1)*NpixelBlock*MapOut + 1 ) : ( b * NpixelBlock * MapOut)  ] 
		Ret$postExpectMu[[b]] <- list()
		for( m in 1:MapOut )
		{
			Ret$postExpectMu[[b]][[m]] <- Matrix( x[ ( (m-1)*NpixelBlock + 1 ) : ( m*NpixelBlock ) ] , nrow=Nside.block, byrow=T  )  
		}
		#residual maps
		if( residual )
		{
			x <- w$residual[ ( (b-1)*NpixelBlock*MapIn + 1 ) : ( b * NpixelBlock * MapIn)  ] 
			Ret$residual[[b]] <- list()
			for( m in 1:MapIn )
			{
				Ret$residual[[b]][[m]] <- Matrix( x[ ( (m-1)*NpixelBlock + 1 ) : ( m*NpixelBlock ) ], nrow=Nside.block, byrow=T )
			}
		}
		#grid points
		g <- cumsum( w$usedGridPts )
		if( b == 1 ) offset <- 0 else offset <- g[b-1]
		l <- sum( par.infer ) + 1
		x <- w$GridPts[ ( offset + 1 ) : w$usedGridPts[b]  ] 
		Ret$GridPts[[b]] <- matrix( x, ncol = l, byrow=T  )
	}

	return( Ret )
}
