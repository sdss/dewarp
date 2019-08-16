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

def genimg(width=8192, height=5120, radius=None, fpslayoutfilename='fps_RTConfig.txt', outfilename='simulatedwarpedfiducials.fits', bgIntensity=3, bgGaussMean=2, bgGaussStdDev=1, fidSuperGaussPeak_min=190, fidSuperGaussPeak_max=210, fidSuperGaussAWAM_min=1, fidSuperGaussAWAM_max=2, fidSuperGaussHWHM_dmin=1, fidSuperGaussHWHM_dmax=4, maxn=5):
    """Draws a warped theoretical image of fiducials
    
    Parameters:
        width,height (int):
            dimensions of the generated image in pixels
        radius (float):
            the radius within which all fiducials reside (in units of fiducials, probably mm)
            if None, we choos a radius so that all of them fit
        fpslayoutfilename (str):
            a path to a config file containing fiducial positions (probably has '.txt' extension)
        outfilename (str):
            a path to a (possibly nonexistent, overwrites of exists) fits image file (must have '.fits' extension)
        bgIntensity (float):
            the pixel value (possibly noninteger) of the background, before noise -- make nonzero so that Poisson noise can act on it
        bgGaussMean (float):
            the mean value for Gaussian noise in addition to bgIntensity, applied before Poisson noising
        bgGaussStdDev (float):
            the standard deviation for Gaussian noise in addition to bgIntensity, applied before Poisson noising
            mustn't be negative
        fidSuperGaussPeak_min (float):
            the greatest lower bound of the uniform distribution of peak intensities for each fiducial randomized
            must not be larger than fidSuperGaussPeak_max
        fidSuperGaussPeak_max (float):
            the least upper bound of the uniform distribution of peak intensities for each fiducial randomized
            must not be smaller than fidSuperGaussPeak_min
        fidSuperGaussAWAM_min (float):
            the least upper bound of the uniform distribution of almost-width-almost-maximum intensities (how many pixels to the side is the supergaussian 255/256 times maximum intensity) for each fiducial randomized
            must not be larger than fidSuperGaussAWAM_max
        fidSuperGaussAWAM_max (float):
            the greatest lower bound of the uniform distribution of almost-width-almost-maximum intensities (how many pixels to the side is the supergaussian 255/256 times maximum intensity) for each fiducial randomized
            must not be smaller than fidSuperGaussAWAM_min
        fidSuperGaussHWHM_dmin (float):
            the least upper bound of the uniform distribution of dhalf-width-half-maximum intensities (how many pixels to the side of the AWAM is the supergaussian half the maximum intensity) for each fiducial randomized
            must not be larger than fidSuperGaussHWHM_max
        fidSuperGaussHWHM_dmax (float):
            the greatest lower bound of the uniform distribution of dhalf-width-half-maximum intensities (how many pixels to the side of the AWAM is the supergaussian half the maximum intensity) for each fiducial randomized
            must not be smaller than fidSuperGaussHWHM_dmin
        maxn (int):
            all basis functions (orthogonal gradient/curl zernikes) with degree maxn or lower (except the trivial piston) are included in warping with random magnitudes
    """
    elapsedTime = time.time()
    
    #parameters and variables
    width = int(width)
    height = int(height)
    bgIntensity = float(bgIntensity)
    bgGaussMean = float(bgGaussMean) #in addition to bgIntensity, set to 0 to add and subtract equally likely and overall add no intensity
    bgGaussStdDev = float(bgGaussStdDev) #bigger means bigger spread (means more random), set to 0 to have no randomness
    units2Pixels = min(width,height)/2 #radius not diameter...
    coefs = opticsmath.warpcoefs()
    
    
    
    print('Loading fiducials from %s'%fpslayoutfilename)
    maxR = 0
    fidxys = ioutils.fiducial_xys_from_file(fpslayoutfilename)
    fidxys = opticsmath.unitize_xys(fidxys, radius)
    fidGaussPeaks = numpy.random.uniform(fidSuperGaussPeak_min,fidSuperGaussPeak_max, int(len(fidxys)/2))
    fidGaussAWAMs = numpy.random.uniform(fidSuperGaussAWAM_min,fidSuperGaussAWAM_max, int(len(fidxys)/2))
    fidGaussHWHMs = fidGaussAWAMs+numpy.random.uniform(fidSuperGaussHWHM_dmin,fidSuperGaussHWHM_dmax, int(len(fidxys)/2))
    fidGaussExponents = numpy.divide(numpy.log2(8-numpy.log2(3)-numpy.log2(5)-numpy.log2(17)),numpy.log2(fidGaussAWAMs)-numpy.log2(fidGaussHWHMs))
    fidRelevantRadii = numpy.ceil(numpy.multiply(fidGaussHWHMs,numpy.power(numpy.log2(255*fidGaussPeaks),numpy.divide(1,fidGaussExponents)))).astype('uint32')
    print('We have %d fiducials'%(len(fidxys)/2))
    
    
    
    print('Computing optical distortions')
    coefs.randomizetransform(int(maxn))
    
    
    
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
        dot = [fidGaussPeaks[i]*math.pow(2,-float(numpy.power(((x+xmins[i]-xcs[i])*(x+xmins[i]-xcs[i])+(y+ymins[i]-ycs[i])*(y+ymins[i]-ycs[i]))/fidGaussHWHMs[i]/fidGaussHWHMs[i],fidGaussExponents[i]/2))) for x,y in numpy.ndindex((xmaxs[i]-xmins[i], ymaxs[i]-ymins[i]))]
        dot = numpy.reshape(dot, (xmaxs[i]-xmins[i], ymaxs[i]-ymins[i]))
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