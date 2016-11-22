from os.path import join, basename
from os import environ
from sys import stderr

import numpy as np
from SetCoverPy import setcover, mathutils
from matplotlib import pyplot as plt
import seaborn as sns
from astropy.io import fits

from redmonster.datamgr.io2 import read_ndArch

# Dict of plates and corresponding MJDS
plates = {8123:56931}

# version of redmonster reductions to use
rmver = 'v1_1_0'

# chi2 threshold to use in set cover distance matrix
mindist = 0.9

rmdir = join( environ['REDMONSTER_SPECTRO_REDUX'], environ['RUN2D'],
             '%s' % rmver)

# loop over all plates
for plate in plates.keys():
    # list of fibers classified as galaxy with zwarning=0
    fibers = []
    zs = []
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
        if fiber != 516:
            if hdu[1].data.ZWARNING[i] == 0:
                if hdu[1].data.CLASS1[i] == 'ssp_galaxy_glob':
                    fibers.append(fiber)
                    zs.append(hdu[1].data.Z1[i])

    # nfibers x nfibers matrix to hold distances
    distmat = np.zeros( (len(fibers), len(fibers)) )

    # binary distance matrix
    #binmat = np.zeros( distmat.shape )


    # calculate distance from each fiber to every other fiber
    for i in range(len(fibers)):
        # set data and such for this fiber
        fiber1 = fibers[i]
        data1 = hduidl[0].data[fiber1]
        sigma1 = hduidl[1].data[fiber1]
        wave1 = wave / (1+zs[i])
        # normalize slow spectrum to mean value between 4500 and 4700 A, and
        # scale errors accordingly
        high1 = np.abs(wave1-4700).argmin()
        low1 = np.abs(wave1-4500).argmin()
        mean1 = np.mean(data1[low1:high1+1])
        data1 /= mean1
        sigma1 *= mean1
        # loop over upper half of matrix and fill bottom half with symmetry
        for j in range(i, len(fibers)):
            fiber2 = fibers[j]
            stderr.write("\rTotal: %s i: %s j:%s" % (len(fibers), i+1, j+1) )
            # convert to rest-frame wavelength
            wave2 = wave / (1+zs[j])
            # crop spectra to only overlapping region
            if wave1[0] < wave2[0]:
                low_bound = np.abs(wave1 - wave2[0]).argmin()
                data1_2 = data1[low_bound:]
                sigma1_2 = sigma1[low_bound:]
                wave1_2 = wave1[low_bound:]
                high_bound = np.abs(wave2 - wave1[-1]).argmin()
                data2 = hduidl[0].data[fiber2][:high_bound+1]
                sigma2 = hduidl[1].data[fiber2][:high_bound+1]
                wave2 = wave2[:high_bound+1]
            elif wave1[0] > wave2[0]:
                low_bound = np.abs(wave2 - wave1[0]).argmin()
                data2 = hduidl[0].data[fiber2][low_bound:]
                sigma2 = hduidl[1].data[fiber2][low_bound:]
                wave2 = wave2[low_bound:]
                high_bound = np.abs(wave1 - wave2[-1]).argmin()
                data1_2 = data1[:high_bound+1]
                sigma1_2 = sigma1[:high_bound+1]
                wave1_2 = wave1[:high_bound+1]
            else: # same redshift (most likely same spectrum)
                data1_2 = data1
                sigma1_2 = sigma1
                data2 = hduidl[0].data[fiber2]
                sigma2 = hduidl[1].data[fiber2]
            # normalize fast spectrum to mean value between 4500 and 4700 A, and
            # scale errors accordingly
            low2 = np.abs(wave2-4500).argmin()
            high2 = np.abs(wave2-4700).argmin()
            mean2 = np.mean(data2[low2:high2+1])
            data2 /= mean2
            sigma2 *= mean2
            # convert inverse variance to sigma squared and add in quadrature
            variance = (1/sigma1_2) + (1/sigma2)
            # calculate reduced chi2 as distance
            distmat[i][j] = distmat[j][i] = np.sum( (data1_2 - data2)**2 /\
                                                   variance ) / data1.shape[0]
            #if distmat[i][j] < mindist:
                #binmat[i][j] = 1
    print " "
    binmat = distmat < mindist
    cost = np.ones(binmat.shape[0])

    g = setcover.SetCover(binmat, cost)
    g.greedy()
    #g.SolveSCP()

    # Get the archetype indices
    iarchetype = np.nonzero(g.s)[0]
    print "Number of archetypes: %s\n" % iarchetype.shape

    # How many are covered by each archetype?
    n_rep = np.sum(binmat[:, iarchetype], axis=0)
    for i,arch in enumerate(iarchetype):
        print "Archtype #%s represents %s spectra." % (arch, n_rep[i])





















