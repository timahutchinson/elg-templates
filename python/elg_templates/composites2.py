from os.path import join, basename
from os import environ
from glob import iglob

import numpy as np
from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel
from scipy.ndimage.filters import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from redmonster.datamgr.io import read_ndArch

version = 'v5_10_0'
platedir = join( environ['REDMONSTER_SPECTRO_REDUX'], '%s_poly1' % version, '*' )

composites = []
waves = []
counts = []
#inds = [9,10,11,12,13,14]
inds = [9,10,11,12,14]
#minspectra = [50,200,1000,1000,1000,1000]
minspectra = [50,1000,1000,1000,1000]
#ages = [.56, 1.00, 1.78, 3.16, 5.62,10.0]
ages = [.56, 1.00, 1.78,3.16,10.0]

for j,ind in enumerate(inds):
    count = np.zeros(20000)
    spectrum = np.zeros(20000)
    wave = 10**(2.000 + np.arange(20000) * .0001)
    for path in iglob(platedir):
        plate = basename(path)
        print plate
        for filename in iglob( join( path, version, '*' ) ):
            if len(basename(filename)) == 26:
                hdu = fits.open(filename)
                #hdu_idl= fits.open( join( environ['BOSS_SPECTRO_REDUX'], version, plate,
                #                          'spPlate%s' % basename(filename)[10:] ) )
                hdu_idl = fits.open(join(environ['BOSS_SPECTRO_REDUX'], 'test/bautista/test_dr14', plate,
                                         'spPlate%s' % basename(filename)[10:]))
                thiswave = 10**(hdu_idl[0].header['COEFF0'] + np.arange(hdu_idl[0].header['NAXIS1']) * .0001)
                for i,zwarn in enumerate(hdu[1].data.ZWARNING):
                    if zwarn == 0:
                        if hdu[1].data.CLASS1[i] == 'ssp_galaxy_glob':
                            if eval(hdu[1].data.MINVECTOR1[i])[0] == ind:
                                try:
                                    offset = np.abs(thiswave[0]/(1+hdu[1].data.Z1[i]) - wave).argmin()
                                    spectrum[offset:offset+hdu_idl[0].header['NAXIS1']] += hdu_idl[0].data[hdu[1].data.FIBERID[i]]
                                    count[offset:offset+hdu_idl[0].header['NAXIS1']] += 1.
                                    print "Spectrum added!"
                                except ValueError:
                                    pass
    spectrum = np.delete(spectrum, np.where(count < minspectra[j]))
    waves.append(np.delete(wave, np.where(count < minspectra[j]) ) )
    counts.append(np.delete(count, np.where(count < minspectra[j])))
    composites.append(spectrum / counts[-1])


x,y,z = read_ndArch( join( environ['REDMONSTER_TEMPLATES_DIR'], 'ndArch-ssp_galaxy_glob-v000.fits' ) )
tempwave = 10**(z['coeff0'] + np.arange(z['nwave']) * z['coeff1'])
sns.set_style('white')
sns.set_palette('muted')
sns.set_context('talk')
f = plt.figure()
f.set_figheight(10)
for i,ind in enumerate(inds):
    offset = np.abs(waves[i][0] - tempwave).argmin()
    this_temp = x[ind][offset:offset+waves[i].shape[0]]
    weights = 1. / counts[i]
    a = np.sum((composites[i]*this_temp)/weights) / np.sum((this_temp**2)/weights) # with weights
    #a = np.sum((composites[i]*this_temp)) / np.sum((this_temp**2)) # without weights)
    model = a * this_temp

    ax = f.add_subplot(len(inds),1,i+1)
    #plt.plot(waves[i], composites[i], label='Composite', color=sns.diverging_palette(10, 220, sep=80, n=7)[-1], linewidth=2)
    #plt.plot(waves[i], composites[i], label='Composite', color=sns.color_palette("coolwarm", 7)[0], linewidth=2)
    firstind, lastind = np.abs(2000-np.array(waves[i])).argmin(), np.abs(5400-np.array(waves[i])).argmin()
    plt.plot(waves[i][firstind:lastind], convolve(composites[i][firstind:lastind], Box1DKernel(5)), label='Composite', color='black')
    #plt.plot(waves[i][firstind:lastind], composites[i][firstind:lastind], label='Composite', color='black')
    with sns.axes_style(rc={"lines.linewidth": 5}):
        #plt.plot(waves[i], model, label='Model', color=sns.diverging_palette(10, 220, sep=80, n=7)[0], linewidth=1)
        #plt.plot(waves[i], model, label='Model', color=sns.color_palette("coolwarm", 7)[-1], linewidth=1)
        #plt.plot(waves[i][firstind:lastind], model[firstind:lastind], label='Model', color=sns.light_palette((210, 90, 60), input="husl")[-3], linewidth=1)
        #plt.plot(waves[i][firstind:lastind], convolve(model[firstind:lastind], Box1DKernel(5)), label='Model', color=sns.light_palette((210, 90, 60), input="husl")[-3], linewidth=1)
        plt.plot(waves[i][firstind:lastind], gaussian_filter(model[firstind:lastind],5), label='Model', color=sns.light_palette((210, 90, 60), input="husl")[-3], linewidth=1)
    ymin = np.sort(composites[i][firstind:lastind])[np.round(composites[i][firstind:lastind].shape[0]*.01)]
    ymax = np.sort(composites[i][firstind:lastind])[np.round(composites[i][firstind:lastind].shape[0]*.99)]*1.1
    plt.axis([2500, waves[i][firstind:lastind][-1], 0, ymax]) # changed ymin to 0
    plt.text(.6*(waves[i][firstind:lastind][-1]-waves[i][firstind:lastind][0])+waves[i][firstind:lastind][0], .05*(ymax-ymin)+ymin, r'%s Gyr template' % ages[i])
    if i == 2:
        plt.ylabel('$f_\lambda$ $(10^{-17}$ erg cm$^{-2}$ s$^{-1}$ $\AA^{-1})$')
    if i == (len(inds)-1):
        plt.xlabel(r'Rest-frame Wavelength ($\AA$)')
    plt.legend(loc=2)
plt.subplots_adjust(hspace = .3)
plt.savefig('/uufs/astro.utah.edu/common/home/u0814744/compute/scratch/comps.pdf')
#    plt.clf()
#    plt.plot(waves[i], counts[i])
#    plt.xlim(waves[i][0], waves[i][-1])
#    plt.savefig('/uufs/astro.utah.edu/common/home/u0814744/compute/scratch/counts%s.pdf' % inds[i])
#    plt.clf()
for j in xrange(len(counts)):
    print "Counts: %s" % np.max(counts[j])
