from os import environ
from os.path import join

import numpy as np

from redmonster.datamgr.io2 import read_ndArch
from redmnster.physics.misc import poly_array

data,parlists,infodict = read_ndArch(join(environ['REDMONSTER_TEMPLATES_DIR'],
                                          'ndArch-ssp_galaxy_noemit-v000.fits'))

temp = data[10,2]

