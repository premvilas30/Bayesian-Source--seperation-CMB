import matplotlib.pyplot as plt
import numpy as np
import pyfits as pf
import healpy as hp

# just take care with the directory paths below-- you will have to 
# change these specific to your own code

mnames = [ 'Ka', 'K', 'Q', 'V', 'W' ]
prefix = '../WMAP_data/wmap_band_imap_r9_9yr_'
fnames = [ prefix + m + '_v5.fits' for m in mnames ]

n = len( mnames )

# outfile names
folder = #enter a directory here'/Plots/'
plotnames = [ folder + m + '_9yr.pdf' for m in mnames ]

from matplotlib import cm

cmap = cm.bwr
cmap.set_under('w')
cmap.set_bad('gray')

mask = hp.read_map('../WMAP_data/wmap_temperature_kq75_analysis_mask_r9_9yr_v5.fits').astype(np.bool)

for k in range( 0 , n ) :
	mapin = hp.read_map( fnames[ k ] )
	mapin_masked = hp.ma( mapin )
	mapin_masked.mask  = np.logical_not( mask )
	hp.mollview( mapin_masked.filled(), title = mnames[k] + ' band map', unit='mK, thermodynamic', cmap=cmap )
	plt.savefig( plotnames[ k ] )


cmap = cm.binary
cmap.set_under('green')
cmap.set_bad('gray')
maskname = '../WMAP_data/wmap_temperature_kq75_analysis_mask_r9_9yr_v5.fits'
mapin = hp.read_map( maskname )
hp.mollview( mapin, title = "Temperature analysis mask", cmap=cmap )
plt.savefig( '../../Dropbox/AstronomyPaper/Plots/Mask.pdf' )
