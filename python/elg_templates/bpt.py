from os.path import join, basename
from os imort environ
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

rmdir = join( environ['REDMONSTER_SPECTRO_REDUX'], environ['RUN2D'],
             '%s' % rmver)

# loop over all plates
for plate in plates.keys():
    # open redmonster file
    hdu = fits.open( join( rmdir, '%s' % plate, 'redmonster-%s-%s.fits'
                          % (plate, plates[plate]) ) )
    # open data file
    idlpath = join( environ['BOSS_SPECTRO_REDUX'], environ['RUN2D'],
                   '%s' % plate, 'spPlate-%s-%s.fits' % (plate, plates[plate]) )
    hduidl = fits.open(idlpath)

    # loop over fibers
    for i,spec in enumerate(hduidl[0].data):
        
