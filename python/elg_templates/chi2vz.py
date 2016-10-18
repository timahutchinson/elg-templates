from os import environ
from os.path import exists, join

import numpy as np
from scipy import linalg
from astropy.io import fits
import matplotlib.pyplot as p

from redmonster.physics.misc import poly_array

plate = 8123
mjd = 56931
npoly = 4
w = 10
rmpath = join( environ['REDMONSTER_SPECTRO_REDUX'], environ['RUN2D'], environ['REDMONSTER_VER'], '%s' % plate, 'redmonster-%s-%s.fits' % (plate,mjd) )
platepath = join( environ['BOSS_SPECTRO_REDUX'], environ['RUN2D'], '%s' % plate, 'spPlate-%s-%s.fits' % (plate,mjd) )

hdu = fits.open(rmpath)
hduplate = fits.open(platepath)
npix = hdu[2].data.shape[-1]
nfibers = hdu[0].header['NFIBERS']

wave = 10**(hdu[0].header['COEFF0'] + np.arange(npix)*hdu[0].header['COEFF1'])

pwave = 10**(3.0 + np.arange(10000)*0.0001)
count = np.zeros(10000)
#chi2 = np.zeros(10000)
chi2 = []
z = []

for i in xrange(1000):
    print i
    if hdu[1].data.ZWARNING[i] == 0:
        if hdu[1].data.CLASS1[i] == 'ssp_galaxy_glob':
            z.append(hdu[1].data.Z1[i])
            spec = hduplate[0].data[i]
            ivar = hduplate[1].data[i]
            thiswave = wave / (1+hdu[1].data.Z1[i])
            #pT = poly_array(npoly, npix)
            #ninv = np.diag(ivar)
            #a = linalg.solve( np.dot(np.dot(pT,ninv),np.transpose(pT)), np.dot(np.dot(pT,ninv), spec) )
            #pmod = np.dot(np.transpose(pT),a)
            #for j in xrange(w/2,npix-w/2):
                #pmodchi2 = np.sum(((spec[j-w/2:j+w/2+1]-pmod[j-w/2:j+w/2+1])**2)*ivar[j-w/2:j+w/2+1]) / float(w)
            j = np.abs(thiswave-3950).argmin()
            w = 500
            tmodchi2 = np.sum(((spec[j-w/2:j+w/2+1]-hdu[2].data[i,0][j-w/2:j+w/2+1])**2)*ivar[j-w/2:j+w/2+1]) / float(w)
                #ind = np.abs(thiswave[j] - pwave).argmin()
                #chi2[ind] += (pmodchi2 - tmodchi2)
                #count[ind] += 1
            chi2.append(tmodchi2)

#y = chi2/count
#for i,val in enumerate(y):
#    if np.isnan(val): y[i] = 0
#p.plot(pwave, y)
#p.axis([pwave[0], pwave[-1], min(y)*2, max(y)*1.2])
p.scatter(z, chi2, marker='.')
p.axis([0,1.4,0,3])
p.xlabel(r'$z$', size=14)
p.ylabel(r'$\chi_\mathrm{Ca H&K}^2$', size=14)
p.savefig('../../plots/chi2vz.pdf')
