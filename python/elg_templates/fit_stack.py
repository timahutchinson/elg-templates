from os import environ
from os.path import join

import numpy as np
from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel
from scipy.linalg import solve
from matplotlib import pyplot as plt
import seaborn as sns
plt.interactive(True)

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
temp = temp[ind:ind+spec.shape[0]]
tempwave = tempwave[ind:ind+spec.shape[0]]

pmat = np.zeros((temp.shape[0],2))
ninv = np.diag(ivar)
pmat[:,0] = temp
pmat[:,1] = np.ones(temp.shape[0])

a = solve(np.dot(np.dot(np.transpose(pmat), ninv), pmat),
          np.dot(np.dot(np.transpose(pmat), ninv), spec))
model = np.dot(pmat, a)

sns.set_palette(sns.color_palette("hls", 8))
sns.set_style('white')
f = plt.figure()
ax = f.add_subplot(211)
plt.plot(wave, convolve(spec, Box1DKernel(5)), color='black')
plt.plot(wave, model, color=sns.color_palette("hls", 8)[0])
plt.ylabel(r'Flux (arbitrary)', size=14)
plt.axis([wave[0], wave[-1], -0.5,4])
ax = f.add_subplot(212)
plt.plot(wave, convolve(spec-model, Box1DKernel(5)), color='black')
plt.ylabel(r'Flux (arbitrary)', size=14)
plt.xlabel(r'Rest-frame wavelength ($\AA$)', size=14)
plt.axis([wave[0], wave[-1], -1,3])
f.tight_layout()
f.savefig('/Users/timhutchinson/compute/repos/elg-templates/plots/stackfits.pdf')

#fit OII line at XXX A and plot from 3700 to 3760
