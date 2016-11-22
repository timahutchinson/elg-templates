from os.path import join, basename
from os import environ
from sys import stderr

import numpy as np
import SetCoverPy as SCP
from matplotlib import pyplot as plt
import seaborn as sns
from astropy.io import fits

from redmonster.datamgr.io2 import read_ndArch

# Dict of plates and corresponding MJDS
plates = {8123:56931}

# version of redmonster reductions to use
rmver = 'v1_1_0'

# chi2 threshold to use in set cover distance matrix
mindist = 1

rmdir = join( environ['REDMONSTER_SPECTRO_REDUX'], environ['RUN2D'],
             '%s' % rmver)

# loop over all plates
for plate in plates.keys():
    # list of fibers classified as galaxy with zwarning=0
    fibers = []
    # open redmonster file
    hdu = fits.open( join( rmdir, '%s' % plate, 'redmonster-%s-%s.fits'
                          % (plate, plates[plate]) ) )
    # open data file
    idlpath = join( environ['BOSS_SPECTRO_REDUX'], environ['RUN2D'],
                   '%s' % plate, 'spPlate-%s-%s.fits' % (plate, plates[plate]) )
    hduidl = fits.open(idlpath)

    # create wavelength solution
    wave = 10**(hduidl[0].header['COEFF0'] + hduidl[0].header['COEFF1'] *
                np.arange(hduidl[0].header['NAXIS1']))

    # loop over fibers and keep those classed as gal with zwarn=0
    for i,fiber in enumerate(hdu[1].data.FIBERID):
        if hdu[1].data.ZWARNING[i] == 0:
            if hdu[1].data.CLASS1[i] == 'ssp_galaxy_glob':
                fibers.append(fiber)

    # nfibers x nfibers matrix to hold distances
    distmat = np.zeros( (len(fibers), len(fibers)) )

    # binary distance matrix
    #binmat = np.zeros( distmat.shape )


    # calculate distance from each fiber to every other fiber
    for i,fiber1 in enumerate(fibers):
        for j,fiber2 in enumerate(fibers):
            stderr.write('\r Total: %s i: %s j:%s' % (len(fibers), i, j) )
            # convert inverse variance to sigma squared and add in quadrature
            sigma2 = (1/hduidl[1].data[fiber1]) + (1/hduidl[1].data[fiber2])
            # calculate chi2 as distance
            distmat[i][j] = np.sum( (hduidl[0].data[fiber1] *
                                     hduidl[0].data[fiber2]) / sigma2 )
            #if distmat[i][j] < mindist:
                #binmat[i][j] = 1
    print ' '
    binmat = distmat < mindist




















