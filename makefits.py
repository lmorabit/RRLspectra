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
    print 'saving fits file '+image+'.fits'
    exportfits(imagename = image, stokeslast=False, fitsimage = 'channel_images/'+image+'.fits', overwrite=True )

os.system('mv L*restored.corr channel_images/restoredcorr/')

