#!/usr/bin/env python

# Author: Leah Morabito (June. 2014)
# must be run as: casapy --nologger -c convolve.py
#     there are memory issues with too many files being open

import glob
import numpy as np
import os

# make a directory for the images
os.mkdir('channel_images')
os.mkdir('channel_images/restoredcorr')

# get a list of images
casaimages = sorted(glob.glob('*restored.corr'))

# save fits files to the images directory
for image in casaimages:

	ia.open(image)
        imagebeam = ia.restoringbeam()
        print image

        print 'imagebeam [major, minor, bpa]: ', imagebeam['major']['value'], imagebeam['minor']['value'],imagebeam['positionangle']['value']
        bmaj = imagebeam['major']['value']
        bmin = imagebeam['minor']['value']
        bpa  = imagebeam['positionangle']['value']

        to_bmaj = '400.0arcsec'
        to_bmin = '400.0arcsec'
        to_bpa  = '0deg'
        print, 'to beam: ', to_bmaj, to_bmin, to_bpa

        # source --> three quantities that describe the source maj axis, min axis, position angle
        # beam --> three quantities that describe the beam maj axis, min axis, position angle
        # convol beam will return TRUE if point source, and the fit paramters
        convolbeam = ia.deconvolvefrombeam(source=[to_bmaj, to_bmin, to_bpa],beam=[str(bmaj)+'arcsec',str(bmin)+'arcsec',str(bpa)+'deg'])

        bmaj = str(float(convolbeam['fit']['major']['value'])) + 'arcsec'
        bmin = str(float(convolbeam['fit']['minor']['value'])) + 'arcsec'
        bpa  = str(float(convolbeam['fit']['pa']['value'])   ) + 'deg'
        print 'covolving beam: ', bmaj, bmin, bpa

        if to_bmaj != bmaj and to_bmin != bmin:
            print 'Convolving.....'
            # axes: the pixel axes to be convolved
            # type: the type of convolution
            # major, minor, pa = the convolving kernal parameters
            s = ia.convolve2d(outfile = image + '.convolve', axes=[0,1], type='gaussian', major=bmaj, minor=bmin, pa=bpa, scale = -1, overwrite=True)
            exportfits(imagename = image + '.convolve', stokeslast=False, fitsimage = 'channel_images/'+image + '.cnv.fits', overwrite=True )
        else:
            print 'Cannot convolve'

        ia.close()

os.system('mv L*restored.corr channel_images/restoredcorr/')

