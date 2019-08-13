# Author: Jake VanderPlas <vanderplas@astro.washington.edu>
# License: BSD
#   The figure is an example from astroML: see http://astroML.github.com
import numpy as np
from matplotlib import pyplot as plt

# warning: due to a bug in healpy, importing it before pylab can cause
#  a segmentation fault in some circumstances.
import healpy as hp

#from astroML.datasets import fetch_wmap_temperatures


#------------------------------------------------------------
# Fetch the data
#wmap_unmasked = fetch_wmap_temperatures(masked=False)
#wmap_masked = fetch_wmap_temperatures(masked=True)
#white_noise = np.ma.asarray(np.random.normal(0, 0.062, wmap_masked.shape))

#------------------------------------------------------------
# plot the unmasked map
#fig = plt.figure(1)
#hp.mollview(wmap_unmasked, min=-1, max=1, title='Unmasked map',
#            fig=1, unit=r'$\Delta$T (mK)')

#------------------------------------------------------------
# plot the masked map
#  filled() fills the masked regions with a null value.
#fig = plt.figure(2)
#hp.mollview(wmap_masked.filled(), title='Masked map',
#            fig=2, unit=r'$\Delta$T (mK)')
            

#-----------------------------------------------------------
#plot the CMB reconstructed map
#             

mask = hp.read_map('../WMAP_data/wmap_temperature_kq75_analysis_mask_r9_9yr_v5.fits').astype(np.bool)

fname1 = '../TEMPLATES/wmap_ilc_9yr_v5.fits'
fname2 = '../WMAP_result/WMAP_result1/cmb_9yr.fits'

fig = plt.figure(1)
mapin = hp.read_map( fname1 )
mapin_masked = hp.ma( mapin )
mapin_masked.mask  = np.logical_not( mask )
hp.mollview( mapin_masked.filled(), title = 'Plot', fig=1 )

hp.mollview( mapin, title = 'Plot', fig=1 )
#fig = plt.figure(2)
#mapin = hp.read_map( fname2 )
#mapin_masked = hp.ma( mapin )
#mapin_masked.mask  = np.logical_not( mask )
#hp.mollview( mapin_masked.filled(), title = 'Plot', fig=2 )
            


#------------------------------------------------------------
# compute and plot the power spectrum
cl = hp.anafast(mapin, lmax=1024)
ell = np.arange(len(cl))

#cl_white = hp.anafast(white_noise, lmax=1024)

fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.scatter(ell, ell * (ell + 1) * cl / (2*pi),
           s=4, c='black', lw=0,
           label='data')
#ax.scatter(ell, ell * (ell + 1) * cl_white,
#           s=4, c='gray', lw=0,
 #          label='white noise')

ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$\ell(\ell+1)C_\ell/2\pi$')
ax.set_title('Angular Power (not mask corrected)')
ax.legend(loc='upper right')
ax.grid()
ax.set_xlim(0, 400)
ax.set_ylim(0,0.003)

plt.show()


