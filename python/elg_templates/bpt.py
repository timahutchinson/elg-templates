from os.path import join, basename
from os import environ
from sys import stderr

import numpy as np
from SetCoverPy import setcover, mathutils
from matplotlib import pyplot as plt
import seaborn as sns
from astropy.io import fits
from scipy.optimize import curve_fit

from redmonster.datamgr.io2 import read_ndArch

# Dict of plates and corresponding MJDS
plates = {8123:56931}

# version of redmonster reductions to use
rmver = 'v1_1_0'

rmdir = join( environ['REDMONSTER_SPECTRO_REDUX'], environ['RUN2D'],
             '%s' % rmver)

residuals = []
zs = []
x = []
y = []

# loop over all plates
for plate in plates.keys():
    # open redmonster file
    hdu = fits.open( join( rmdir, '%s' % plate, 'redmonster-%s-%s.fits'
                          % (plate, plates[plate]) ) )
    # open data file
    idlpath = join( environ['BOSS_SPECTRO_REDUX'], environ['RUN2D'],
                   '%s' % plate, 'spPlate-%s-%s.fits' % (plate, plates[plate]) )
    hduidl = fits.open(idlpath)
    
    # wavelength solution
    wave = 10**(hduidl[0].header['COEFF0'] + hduidl[0].header['COEFF1'] *
                np.arange(hduidl[0].header['NAXIS1']))

    # loop over fibers
    for i,spec in enumerate(hduidl[0].data):
        if hdu[1].data.ZWARNING[i] == 0:
            if hdu[1].data.CLASS1[i] == 'ssp_galaxy_glob':
                residual = spec - hdu[2].data[i,0]
                this_wave = wave / (1 + hdu[1].data.Z1[i])
                nii = np.abs(wave - 6548.86).argmin()
                ha = np.abs(wave - 6564.61).argmin()
                oiii = np.abs(wave - 5008.24).argmin()
                hb = np.abs(wave - 4862.68).argmin()
                x.append( np.log10( np.sum( residual[nii-3:nii+4] ) /\
                                  np.sum( residual[ha-3:ha+4] ) ) )
                y.append( np.log10( np.sum( residual[oiii-3:oiii+4] ) /\
                                  np.sum( residual[hb-3:hb+4] ) ) )





def gaussian(A, sigma, lambda0, x):
    return A * np.exp(-(x-lambda0)**2/sigma**2)











