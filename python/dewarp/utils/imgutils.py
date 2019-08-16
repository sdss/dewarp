# encoding: utf-8
#
# @Author:    Adam Mendenhall
# @Date:      August 16, 2019
# @Filename:  imgutils.py
# @License:   BSD 3-Clause
# @Copyright: Adam Mendenhall
#

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from astropy.io import fits
import numpy
from dewarp.utils import opticsmath

def writetoimage(filename, data, clobber=True, data_is_rowmajorx=True):
    """Writes data to an image
    
    Parameters:
        filename (str):
            the path to a (possible nonexistent) fits image file (must have extension '.fits')
            if clobber is False, may not write if the file already exists
        data (ndarray):
            a 2d numpy array containing image data
        clobber (bool):
            True if whether or not the file already exists should be ignored (if False and file already exists, does not write)
        data_is_rowmajorx (bool):
            if True, indicates that data[x][0] corresponds to a pixel with x-coordinate x and y-coordinate 0 -- fits images use y row major so transposes before writing
            if False, indicates that data[0][x] corresponds to a pixel with x-coordinate x and y-coordinate 0 -- fits images use this
    """
    if data_is_rowmajorx:
        hdu = fits.PrimaryHDU(data.T) #transposed since that's the way the axes go
    else:
        hdu = fits.PrimaryHDU(data)
    hdu.writeto(outfilename, overwrite=clobber)

def readimage(filename):
    """Reads a fits image as a numpy array
    
    Parameters:
        filename (str):
            the path to a fits image file (must have extension '.fits')
    
    Returns:
        data (ndarray):
            a 2d numpy array containing image data
    """
    data = None
    with fits.open(filename) as hdul:
        data = hdul[0].data
    return data

def centroids_hullbijected(data, expectedxys):
    """Finds centroids in an image and matches them to the expectedxys
    
    if the amount of centroids founds is not the same as the amount of points in expectedxys, None is returned

    Parameters:
        data (ndarray):
            a 2d numpy array containing image data, assumed to be x row major
        expectedxys (list):
            a list of interleaved xy coordinates
            these are the centroids that are expected in the image

    Returns:
        xys (list):
            a list of interleaved xy positions, the same length as expectedxys
            these are the positions of the centroids in the image, ordered 'just like' the expectedxys
    """
    copy_centroidxys = centroids(filename)
    if len(copy_centroidxys) is not len(expectedxys):
        return None
    copy_expectedxys = [a for a in expectedxys]
    centroidxys = [None]*len(copy_centroidxys)
    count = 0
    while count<len(copy_centroidxys):
        expectedhullidxs = opticsmath.untrueconvexhull(copy_expectedxys)
        centroidhullidxs = opticsmath.untrueconvexhull(copy_centroidxys)
        if len(expectedhullidxs)!=len(centroidhullidxs):
            return None
        for i in range(len(expectedhullidxs)):
            centroidxys[expectedhullidxs[i]] = copy_centroidxys[centroidhullidxs[i]]
            centroidxys[expectedhullidxs[i]+1] = copy_centroidxys[centroidhullidxs[i]+1]
            copy_centroidxys[centroidhullidxs[i]] = None
            copy_centroidxys[centroidhullidxs[i]+1] = None
            copy_expectedxys[expectedhullidxs[i]] = None
            copy_expectedxys[expectedhullidxs[i]+1] = None
            count += 1
    return centroidxys

def centroids(data):
    """Finds centroids in an image

    Parameters:
        data (ndarray):
            a 2d numpy array containing image data, assumed to be x row major
    
    Returns:
        xys (list):
            a list of interleaved xy positions
            these are the positions of the centroids in the image
    """
    pass

def genimg(unitary_xys, width=8192, height=5210, outfilename='simulatedwarpedfiducials.fits', bgIntensity=3, bgGaussMean=2, bgGaussStdDev=1, superGaussPeak_min=190, superGaussPeak_max=210, superGaussAWAM_min=1, superGaussAWAM_max=2, superGaussHWHM_dmin=1, superGaussHWHM_dmax=4, maxn=5):
    """Draws a warped theoretical image of dots
    
    assumes the list unitary_xys contains interleaved xy coordinates of dots to draw, all within the unit square

    the dots are warped with random coefficients before drawing
    for each dot, draws a supergaussian with random peak intensity, HWHM and plateau width
    background is sum of constant and gaussian distribution
    entire image is poisson noised after everything

    Parameters:
        unitary_xys (list):
            an interleaved array of xy coordaintes, each point must be in the range [-1,1] cross [-1,1]
        width,height (int):
            dimensions of the generated image in pixels
        outfilename (str):
            a path to a (possibly nonexistent, overwrites of exists) fits image file (must have '.fits' extension)
        bgIntensity (float):
            the pixel value (possibly noninteger) of the background, before noise -- make nonzero so that Poisson noise can act on it
        bgGaussMean (float):
            the mean value for Gaussian noise in addition to bgIntensity, applied before Poisson noising
        bgGaussStdDev (float):
            the standard deviation for Gaussian noise in addition to bgIntensity, applied before Poisson noising
            mustn't be negative
        superGaussPeak_min (float):
            the greatest lower bound of the uniform distribution of peak intensities for each dot randomized
            must not be larger than superGaussPeak_max
        superGaussPeak_max (float):
            the least upper bound of the uniform distribution of peak intensities for each dot randomized
            must not be smaller than superGaussPeak_min
        superGaussAWAM_min (float):
            the least upper bound of the uniform distribution of almost-width-almost-maximum intensities (how many pixels to the side is the supergaussian 255/256 times maximum intensity) for each dot randomized
            must not be larger than superGaussAWAM_max
        superGaussAWAM_max (float):
            the greatest lower bound of the uniform distribution of almost-width-almost-maximum intensities (how many pixels to the side is the supergaussian 255/256 times maximum intensity) for each dot randomized
            must not be smaller than superGaussAWAM_min
        superGaussHWHM_dmin (float):
            the least upper bound of the uniform distribution of dhalf-width-half-maximum intensities (how many pixels to the side of the AWAM is the supergaussian half the maximum intensity) for each dot randomized
            must not be larger than superGaussHWHM_max
        superGaussHWHM_dmax (float):
            the greatest lower bound of the uniform distribution of dhalf-width-half-maximum intensities (how many pixels to the side of the AWAM is the supergaussian half the maximum intensity) for each dot randomized
            must not be smaller than superGaussHWHM_dmin
        maxn (int):
            all basis functions (orthogonal gradient/curl zernikes) with degree maxn or lower (except the trivial piston) are included in warping with random magnitudes
    """
    import time
    elapsedTime = time.time()

    #parameters and variables
    width = int(width)
    height = int(height)
    bgIntensity = float(bgIntensity)
    bgGaussMean = float(bgGaussMean) #in addition to bgIntensity, set to 0 to add and subtract equally likely and overall add no intensity
    bgGaussStdDev = float(bgGaussStdDev) #bigger means bigger spread (means more random), set to 0 to have no randomness
    units2Pixels = min(width,height)/2 #radius not diameter...
    coefs = opticsmath.warpcoefs()



    print('Computing dot gaussians')
    maxR = 0
    gaussPeaks = numpy.random.uniform(superGaussPeak_min,superGaussPeak_max, int(len(unitary_xys)/2))
    gaussAWAMs = numpy.random.uniform(superGaussAWAM_min,superGaussAWAM_max, int(len(unitary_xys)/2))
    gaussHWHMs = gaussAWAMs+numpy.random.uniform(superGaussHWHM_dmin,superGaussHWHM_dmax, int(len(unitary_xys)/2))
    gaussExponents = numpy.divide(numpy.log2(8-numpy.log2(3)-numpy.log2(5)-numpy.log2(17)),numpy.log2(gaussAWAMs)-numpy.log2(gaussHWHMs))
    relevantRadii = numpy.ceil(numpy.multiply(gaussHWHMs,numpy.power(numpy.log2(255*gaussPeaks),numpy.divide(1,gaussExponents)))).astype('uint32')
    print('We have %d dots'%(len(unitary_xys)/2))



    print('Computing optical distortions')
    coefs.randomizetransform(int(maxn))



    print('Creating image of size %dx%d'%(width,height))
    data = numpy.full((width, height), float(bgIntensity))



    print('Drawing warped dots...')
    warpedunitary_xys = coefs.applytransform(unitary_xys)
    xcs = (numpy.multiply(warpedunitary_xys[0::2],units2Pixels)+(width /2)).astype('int32')
    ycs = (numpy.multiply(warpedunitary_xys[1::2],units2Pixels)+(height/2)).astype('int32')
    xmins = numpy.maximum(0,numpy.minimum(width -1,xcs-relevantRadii)).astype('uint32')
    xmaxs = numpy.maximum(0,numpy.minimum(width -1,xcs+relevantRadii)).astype('uint32')
    ymins = numpy.maximum(0,numpy.minimum(height-1,ycs-relevantRadii)).astype('uint32')
    ymaxs = numpy.maximum(0,numpy.minimum(height-1,ycs+relevantRadii)).astype('uint32')
    for i in range(0,len(xcs)):
        dot = [gaussPeaks[i]*numpy.power(2,-float(numpy.power(((x+xmins[i]-xcs[i])*(x+xmins[i]-xcs[i])+(y+ymins[i]-ycs[i])*(y+ymins[i]-ycs[i]))/gaussHWHMs[i]/gaussHWHMs[i],gaussExponents[i]/2))) for x,y in numpy.ndindex((xmaxs[i]-xmins[i], ymaxs[i]-ymins[i]))]
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
    imgutils.writetoimage(outfilename, data)



    elapsedTime = time.time()-elapsedTime
    print('Finished in %f seconds'%elapsedTime)
