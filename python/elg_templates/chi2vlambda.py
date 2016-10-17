from os import environ
from os.path import exists, join

import numpy as np
from astropy.io import fits

from redmonster.physics.misc import poly_array

plate = 8123
mjd = 56931
npoly = 4
rmpath = join( environ['REDMONSTER_SPECTRO_REDUX'], environ['RUN2D'], environ['REDMONSTER_VER'], '%s' % plate, 'redmonster-%s-%s.fits' % (plate,mjd) )
platepath = join( environ['BOSS_SPECTRO_REDUX'], environ['RUN2D'], '%s' % plate, 'spPlate-%s-%s.fits' % (plate,mjd) )

hdu = fits.open(rmpath)
hduplate = fits.open(platepath)
npix = hdu[2].data.shape[1]
nfibers = hdu[0].header['NFIBERS']

wave = 10**(hdu[0].header['COEFF0'] + np.arange(npix)*hdu[0].header['COEFF1'])

for i in xrange(nfibers):
    if hdu[1].data.ZWARNING == 0:
        pmat = np.zeros(npix, nfibers)
        ninv = np.diagonal(hduplate[1].data[i])
