# encoding: utf-8
#
# @Author:    Adam Mendenhall
# @Date:      June 6, 2019
# @Filename:  genimg.py
# @License:   BSD 3-Clause
# @Copyright: Adam Mendenhall
#

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from astropy.io import fits
import numpy
import math
import time
from .utils import opticsmath
from .utils import ioutils

def genimg():
    elapsedTime = time.time()
    
    #parameters and variables
    width = 8192
    height = 5120
    fpslayoutfilename = 'fps_RTConfig.txt'
    outfilename = 'simulatedwarpedfiducials.fits'
    bgIntensity = 3.
    bgGaussMean = 0. #in addition to bgIntensity, set to 0 to add and subtract equally likely and overall add no intensity
    bgGaussStdDev = 0. #bigger means bigger spread (means more random), set to 0 to have no randomness
    units2Pixels = min(width,height)/2 #radius not diameter...
    coefs = opticsmath.warpcoefs()
    
    
    
    print('Loading fiducials from %s'%fpslayoutfilename)
    maxR = 0
    #allInstruments = specificInstrumentsXYsFromFile(fpslayoutfilename, lambda x: True)
    fidxys = ioutils.fiducial_xys_from_file('fps_RTConfig.txt')
    fidxys = opticsmath.unitize_xys(fidxys, None)
    fidGaussPeaks = numpy.random.uniform(190,210, int(len(fidxys)/2))
    fidGaussHWHMs = numpy.random.uniform(4,6, int(len(fidxys)/2))
    fidRelevantRadii = numpy.ceil(fidGaussHWHMs*numpy.sqrt(numpy.log2(fidGaussPeaks))).astype('uint32')
    print('We have %d fiducials'%(len(fidxys)/2))
    
    
    
    print('Computing optical distortions')
    coefs.randomizetransform(4)
    
    
    
    print('Creating image of size %dx%d'%(width,height))
    data = numpy.full((width, height), float(bgIntensity))
    
    
    
    print('Drawing warped fiducials...')
    warpedfidxys = coefs.applytransform(fidxys)
    xcs = (numpy.multiply(warpedfidxys[0::2],units2Pixels)+(width /2)).astype('int32')
    ycs = (numpy.multiply(warpedfidxys[1::2],units2Pixels)+(height/2)).astype('int32')
    xmins = numpy.maximum(0,numpy.minimum(width -1,xcs-fidRelevantRadii)).astype('uint32')
    xmaxs = numpy.maximum(0,numpy.minimum(width -1,xcs+fidRelevantRadii)).astype('uint32')
    ymins = numpy.maximum(0,numpy.minimum(height-1,ycs-fidRelevantRadii)).astype('uint32')
    ymaxs = numpy.maximum(0,numpy.minimum(height-1,ycs+fidRelevantRadii)).astype('uint32')
    for i in range(0,len(xcs)):
        dot = numpy.reshape([fidGaussPeaks[i]*math.pow(2,-float((x+xmins[i]-xcs[i])*(x+xmins[i]-xcs[i])+(y+ymins[i]-ycs[i])*(y+ymins[i]-ycs[i]))/fidGaussHWHMs[i]/fidGaussHWHMs[i]) for x,y in numpy.ndindex((xmaxs[i]-xmins[i], ymaxs[i]-ymins[i]))], (xmaxs[i]-xmins[i], ymaxs[i]-ymins[i]))
        data[xmins[i]:xmaxs[i], ymins[i]:ymaxs[i]] += dot
    
    
    print('Adding Gaussian noise...')
    if bgGaussStdDev>0:
        data += numpy.random.normal(bgGaussMean, bgGaussStdDev, width*height).reshape((width, height))
    elif bgGaussMean!=0:
        data += bgGaussMean
    print('Adding Poisson noise...')
    data = numpy.maximum(data, 0)
    data = numpy.random.poisson(data.reshape(width*height), width*height).reshape((width, height))
    
    
    
    print('Writing data to %s'%outfilename)
    hdu = fits.PrimaryHDU(data.T) #transposed since that's the way the axes go
    hdu.writeto(outfilename, overwrite=True)
    
    
    
    elapsedTime = time.time()-elapsedTime
    print('Finished in %f seconds'%elapsedTime)