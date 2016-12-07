from os import environ
from os.path import join

import numpy as np
from astropy.io import fits
from scipy.linalg import solve

from redmonster.datamgr.io2 import read_ndArch
from redmonster.physics.misc import poly_array

data,parlists,infodict = read_ndArch(join(environ['REDMONSTER_TEMPLATES_DIR'],
                                          'ndArch-ssp_galaxy_noemit-v000.fits'))

temp = data[10,2]
tempwave = 10**(infodict['coeff0'] + 0.0001 * np.arange(temp.shape[-1]))

hdu = fits.open('/Users/timhutchinson/Desktop/test.fits')
spec = hdu[0].data[0]
ivar = hdu[1].data[0]
wave = 10**(hdu[0].header['COEFF0'] + 0.0001 * \
                np.arange(hdu[0].data.shape[-1]))

ind = np.abs(tempwave - wave[0]).argmin()
print ind
temp = temp[ind:ind+spec.shape[0]]

pmat = np.zeros((temp.shape[0],2))
ninv = np.diag(ivar)
pmat[:,0] = temp
pmat[:,1] = np.ones(temp.shape[0])

a = solve(np.dot(np.dot(np.transpose(pmat), ninv), pmat), np.dot(np.dot(np.transpose(pmat), ninv), spec))
model = np.dot(pmat, a)
