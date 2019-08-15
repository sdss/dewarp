# encoding: utf-8
#
# @Author:    Adam Mendenhall
# @Date:      June 6, 2019
# @Filename:  genimg.py
# @License:   BSD 3-Clause
# @Copyright: Adam Mendenhall
#
# References:
#   Chunyu Zhao and James H. Burge, Part I https://www.osapublishing.org/oe/abstract.cfm?uri=oe-15-26-18014, Part II https://www.osapublishing.org/oe/abstract.cfm?uri=oe-16-9-6586
#

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from astropy.io import fits
import numpy
import math
import time
import opticsmath
from utils import ioutils

elapsedTime = time.time()

#parameters and variables
width = 1920#8192
height = 1080#5120
fpslayoutfilename = 'fps_RTConfig.txt'
outfilename = 'img.fits'
bgIntensity = 3.
bgGaussMean = 0. #in addition to bgIntensity, set to 0 to add and subtract equally likely and overall add no intensity
bgGaussStdDev = 0. #bigger means bigger spread (means more random), set to 0 to have no randomness
entireRadius = -1 #if positive, all fiducials are assumed to be within a circle of this radius, otherwise it is computed
units2Pixels = min(width,height)/2 #radius not diameter...
coefs = opticsmath.warpcoefs()



print('Loading fiducials from %s'%fpslayoutfilename)
maxR = 0
#allInstruments = specificInstrumentsXYsFromFile(fpslayoutfilename, lambda x: True)
fidxys = ioutils.fiducial_xys_from_file('fps_RTConfig.txt')
fidGaussPeaks = numpy.full((int(len(fidxys)/2)), 200.) #TODO: randomize
fidGaussHWHMs = numpy.full((int(len(fidxys)/2)), 5.) #TODO: randomize TODO: make these units not pixels, but fiducial units
fidRelevantRadii = numpy.ceil(fidGaussHWHMs*numpy.sqrt(numpy.log2(fidGaussPeaks))).astype('uint32')
if entireRadius<=0:
    #circular margin of 30 units
    entireRadius = math.sqrt(numpy.nanmax(numpy.square(fidxys[0::2])+numpy.square(fidxys[1::2])))+30
print('We have %d fiducials within a radius (in fiducial units) of %f'%(len(fidxys)/2,entireRadius))



print('Computing optical distortions')
coefs.randomizetransform(4)



print('Creating image of size %dx%d'%(width,height))
data = numpy.full((width, height), float(bgIntensity))



print('Drawing warped fiducials...')
fidxys = numpy.divide(fidxys,entireRadius)
warpedfidxys = coefs.applytransform(fidxys)
xcs = (numpy.multiply(warpedfidxys[0::2],units2Pixels)+(width /2)).astype('int32')
ycs = (numpy.multiply(warpedfidxys[1::2],units2Pixels)+(height/2)).astype('int32')
xmins = numpy.maximum(0,numpy.minimum(width -1,xcs-fidRelevantRadii)).astype('uint32')
xmaxs = numpy.maximum(0,numpy.minimum(width -1,xcs+fidRelevantRadii)).astype('uint32')
ymins = numpy.maximum(0,numpy.minimum(height-1,ycs-fidRelevantRadii)).astype('uint32')
ymaxs = numpy.maximum(0,numpy.minimum(height-1,ycs+fidRelevantRadii)).astype('uint32')
#perhaps there is a fast numpy way to add each of these fiducial dot blocks into the image.  The following commented lines are not it TODO: settle for for or strive for greatness
#dot = numpy.reshape([numpy.multiply(numpy.power(2,-(numpy.square(x+xmins-xcs)+numpy.square(y+ymins-ycs)).astype('float32').divide(numpy.square(fidGaussHWHMs))),fidGaussPeaks) for x,y in numpy.ndindex((xmaxs-xmins, ymaxs-ymins))], (xmaxs-xmins, ymaxs-ymins))
#data[xmins:xmaxs, ymins:ymaxs] += dot
for i in range(0,len(xcs)):
    dot = numpy.reshape([fidGaussPeaks[i]*math.pow(2,-float((x+xmins[i]-xcs[i])*(x+xmins[i]-xcs[i])+(y+ymins[i]-ycs[i])*(y+ymins[i]-ycs[i]))/fidGaussHWHMs[i]/fidGaussHWHMs[i]) for x,y in numpy.ndindex((xmaxs[i]-xmins[i], ymaxs[i]-ymins[i]))], (xmaxs[i]-xmins[i], ymaxs[i]-ymins[i])) #casts to float needed to get precision
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
# Creates large fits image without storing entire array in memory, TODO: perhaps implement this (or decide not to)
#dummyData = numpy.zeros((1,1), dtype=numpy.float64)
#hdu = fits.PrimaryHDU(data=dummyData)
#header = hdu.header
#while len(header) < (36*4-1):
#    header.append()
#header['NAXIS1'] = width
#header['NAXIS2'] = height
#header.tofile(outfilename, overwrite=True)
#with open(outfilename, 'rb+') as f:
#    f.seek(len(header.tostring()) + (width*height*numpy.abs(header['BITPIX']//8)) - 1)
#    f.write(b'\0')



elapsedTime = time.time()-elapsedTime
print('Finished in %f seconds'%elapsedTime)