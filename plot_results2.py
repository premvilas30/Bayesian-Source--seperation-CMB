# script for making diagnostics of results

# required packages

import matplotlib.pyplot as plt
import numpy as np
import pyfits as pf
import healpy as hp


dn = [ 'K', 'Ka', 'Q', 'V', 'W' ]

resn = [ 'K_residual', 'Ka_residual', 'Q_residual', 'V_residual', 'W_residual' ]

residnames = [ '../WMAP_result/WMAP_result1/' + m + '_9yr.fits' for m in resn ]
datnames = [ '../WMAP_data/wmap_band_imap_r9_9yr_' + m + '_v5.fits' for m in dn ]
plotnames = [ '../../Dropbox/AstronomyPaper/Plots/studentized_' + m  + '_pyr.jpgw' for m in resn ]

n = len( dn )

msk = pf.open( '../WMAP_data/wmap_temperature_kq75_analysis_mask_r9_9yr_v5.fits' )
z = msk[1].data['TEMPERATURE']
idx = np.where( z == 1 )[0]

for k in range( 0, n ) :
	td = pf.open( datnames[k] )
	x = td[1].data['TEMPERATURE']
	x1 = x[ idx ]
	tr = pf.open( residnames[k] )
	y = tr[1].data['SOURCE']
	y1 = y[ idx ]
	plt.plot( x1, y1, 'k.' )
	plt.plot([x1.min(), x1.max()], [0, 0], 'r-', lw=2 )
	plt.xlabel( 'observed value' , fontsize = 18 )
	plt.ylabel( 'standardized residual', fontsize = 18 )
	plt.title( dn[k] + ' channel ' , fontsize = 20 )
	plt.savefig( plotnames[ k ] )
