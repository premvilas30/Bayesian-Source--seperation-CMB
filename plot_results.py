# script for making plots of results

# required packages

import matplotlib.pyplot as plt
import numpy as np
import pyfits as pf
import healpy as hp
from matplotlib import cm

wmap = 1

# string descriptors

if wmap == 1:
	#mnames = [ 'cmb', 'sync', 'dust', 'ffem', 'cmb_precision', 'sync_precision', 'dust_precision', 'ffem_precision', 'sync_ind', 'dust_ind', 'K_residual', 'Ka_residual', 'Q_residual', 'V_residual', 'W_residual' ]
	mnames = [ 'cmb', 'sync', 'dust', 'ffem', 'cmb_precision', 'sync_precision','dust_precision', 'ffem_precision','sync_ind', 'dust_ind', 'K_residual', 'Ka_residual', 'Q_residual', 'V_residual', 'W_residual' ]
	folder = '../WMAP_result/'
	#minv  = [ -.4, -.1, -.1, -.1 ]
	#maxv = [ .4, .8, .2, .75]
	mask = hp.read_map('../WMAP_data/wmap_temperature_kq75_analysis_mask_r9_9yr_v5.fits').astype(np.bool)
	fitsnames = [ folder + m + '_9yr.fits' for m in mnames ]
	ninmap = 5
else:
	mnames = [ 'cmb', 'sync', 'dust', 'ffem', 'cmb_precision', 'sync_precision', 'dust_precision', 'ffem_precision', 'sync_ind', 'dust_ind', 'SkyMap_030_residual', 'SkyMap_044_residual', 'SkyMap_070_residual', 'SkyMap_100_residual', 'SkyMap_143_residual', 'SkyMap_217_residual', 'SkyMap_353_residual' ]
	folder = '../PLANCK_result/'
	mask = hp.read_map('../PLANCK_data/PlanckMask2.fits').astype(np.bool)
	fitsnames = [ folder + m + '.fits' for m in mnames ]
	ninmap = 7


n = len( mnames )

#colour map
cmap = cm.gnuplot
cmap.set_under('w')
cmap.set_bad('gray')


# outfile names
plotnames = [ folder + 'Plots/' +  m + '_9yr.pdf' for m in mnames ]

nicenames  = ['CMB', 'Synchrotron', 'Galactic dust','Free free' , 'CMB precision', 'Synchrotron precision', 'Dust precision', 'Free free precision', 'Synchrotron index', 'Galactic dust index', 'K band residual', 'Ka band residual', 'Q band residual', 'V band residual', 'W band residual']


#colour map bwr for plotting the sources 
cmap = cm.bwr
cmap.set_under('w')
cmap.set_bad('gray')

for k in range( 0 , 1 ) :
	mapin = hp.read_map( fitsnames[ k ] )
	mapin_masked = hp.ma( mapin )
	mapin_masked.mask  = np.logical_not( mask )
	hp.mollview( mapin_masked.filled(), title = nicenames[k]) #cmap=cmap ) #, min=-.35, max=.35 )
	plt.savefig( plotnames[ k ] )


for k in range( 1 , 4 ) :
	mapin = hp.read_map( fitsnames[ k ] )
	mapin_masked = hp.ma( mapin )
	mapin_masked.mask  = np.logical_not( mask )
	hp.mollview( mapin_masked.filled(), title = nicenames[k]) #cmap=cmap )
	plt.savefig( plotnames[ k ] )

#colour map gnuplot for precisions
cmap = cm.gnuplot
cmap.set_under('w')
cmap.set_bad('gray')

for k in range( 4 , 9 ) :
	mapin = hp.read_map( fitsnames[ k ] )
	mapin_masked = hp.ma( mapin )
	mapin_masked.mask  = np.logical_not( mask )
	hp.mollview( mapin_masked.filled(), title = nicenames[k]) #cmap=cmap )
	plt.savefig( plotnames[ k ] )

#colour map gnuplot for precisions
cmap = cm.bwr
cmap.set_under('w')
cmap.set_bad('gray')

for k in range( 9 , 15 ) :
	mapin = hp.read_map( fitsnames[ k ] )
	mapin_masked = hp.ma( mapin )
	mapin_masked.mask  = np.logical_not( mask )
	hp.mollview( mapin_masked.filled(), title = nicenames[k] ) #cmap=cmap )
	plt.savefig( plotnames[ k ] )
	
#write_map( "haslam_nested.fits", mapin, nest = True, fits_IDL=False)


