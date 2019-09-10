# encoding: utf-8
#
# @Author:    Adam Mendenhall
# @Date:      August 21, 2019
# @Filename:  dewarp.py
# @License:   BSD 3-Clause
# @Copyright: Adam Mendenhall
#

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from dewarp.utils import opticsmath
from dewarp.utils import ioutils
from dewarp.utils import imgutils
import pkg_resources

def dewarp(fiducialsimgfile, fibersimgfile, fpslayoutfilename=None, fpsradius=350, genmask=True):
    """Computes the optical warp based on an image of fiducials and applys it to an image of fibers
    
    Parameters:
        fiducialsimgfile (str):
            path to a fits image file containing illuminated fiducials (must have '.fits' extension)
        fibersimgfile
            path to a fits image file containing illuminated metrology fibers (must have '.fits' extension)
        fpslayoutfilename (str):
            path to a config file containing fiducial positions (probably has '.txt' extension)
        fpsradius (float):
            the radius (nonnegative) that all the dots are within, in whatever units the file would like (probably mm)

    Returns:
        xys (list):
            a list of interleaved xy coordinates of the dots found in the image, unwarped based on coefs
    """
    if fpslayoutfilename is None:
        fpslayoutfilename = pkg_resources.resource_filename('dewarp', 'etc/fps_RTConfig.txt')
    coefs = detectwarp(fpslayoutfilename, fpsradius, fiducialsimgfile, genmask=genmask)
    return applywarp(coefs, fibersimgfile)

def detectwarp(fpslayoutfilename=None, fpsradius=350, infilename=None, genmask=True, maskdilation=80):
    """Computes the optical warp of an image based on fiducials
    
    Parameters:
        fpslayoutfilename (str):
            path to a config file containing fiducial positions (probably has '.txt' extension)
        fpsradius (float):
            the radius (nonnegative) that all the dots are within, in whatever units the file would like (probably mm)
        infilename (str):
            path to a fits image file (must have '.fits' extension)

    Returns:
        coefs (warpcoefs):
            an object containing coefficents and functions to dewarp points, probably to be used in applywarp()
    """
    if fpslayoutfilename is None:
        fpslayoutfilename = pkg_resources.resource_filename('dewarp', 'etc/fps_RTConfig.txt')
    if infilename is None:
        infilename = pkg_resources.resource_filename('dewarp', 'etc/simulatedwarpedfiducials.fits')
    ideal_xys = ioutils.fiducial_xys_from_file(fpslayoutfilename)
    ideal_xys = opticsmath.unitize_xys(ideal_xys, fpsradius)
    imgdata = imgutils.readimage(infilename)
    if genmask:
        mask = numpy.ones(imgdata.shape)
        for i in range(0,len(ideal_xys),2):
            x = ideal_xys[i]
            y = ideal_xys[i+1]
            mask[max(0,x-maskdilation):min(mask.shape[0]-1,x+maskdilation):1,max(0,y-maskdilation):min(mask.shape[1]-1,y+maskdilation):1] = 1
        imgutils.writetoimage(pkg_resources.resource_filename('dewarp', 'etc/mask.fits'), mask)
    else:
        mask = imgutils.readimage(pkg_resources.resource_filename('dewarp', 'etc/mask.fits'))
    if mask.shape[0]!=imgdata.shape[0] or mask.shape[1]!=imgdata.shape[1]:
        mask = None
    observed_xys = imgutils.unitdiskcentroids(imgdata, mask)
    observed_xys = opticsmath.sort_closest(observed_xys, ideal_xys)
    coefs = opticsmath.warpcoefs(observed_xys, ideal_xys, .01)
    return coefs

def applywarp(coefs, infilename=None):
    """Applys optical warp to an image
    
    Parameters:
        coefs (warpcoefs):
            an object containing coefficients and functions to dewarp points, probably computed using detectwarp()
        infilename (str):
            path to a fits image file (must have '.fits' extension)

    Returns:
        xys (list):
            a list of interleaved xy coordinates of the dots found in the image, unwarped based on coefs
    """
    if infilename is None:
        infilename = pkg_resources.resource_filename('dewarp', 'etc/simulatedwarpedfiducials.fits')
    mask = imgutils.readimage(pkg_resources.resource_filename('dewarp', 'etc/mask.fits'))
    imgdata = imgutils.readimage(infilename)
    if mask.shape[0]!=imgdata.shape[0] or mask.shape[1]!=imgdata.shape[1]:
        mask = None
    xys = imgutils.unitdiskcentroids(imgdata, mask)
    return coefs.applytransform(xys)

def fakewarp(fpslayoutfilename=None, fpsradius=350, whichinstruments=lambda x:x.lower()=='fiducial', outfilename=None, width=8192, height=5210, bgIntensity=3, bgGaussMean=2, bgGaussStdDev=1, superGaussPeak_min=190, superGaussPeak_max=210, superGaussAWAM_min=1, superGaussAWAM_max=2, superGaussHWHM_dmin=1, superGaussHWHM_dmax=4, coefs=None):
    """Generates a warped image based on a configuration file
    
    Parameters:
        fpslayoutfilename (str):
            the path to a config file containing positions of instruments (probably has '.txt' extension)
        fpsradius (float):
            the radius (nonnegative) that all the dots are within, in whatever units the file would like (probably mm)
        whichinstruments (lambda function of str):
            returns True if evaluated at the name of an instrument to include
        outfilename (str):
            the path to a (possibly nonexistent, overwrites of exists) fits image file (must have '.fits' extension)
    """
    if fpslayoutfilename is None:
        fpslayoutfilename = pkg_resources.resource_filename('dewarp', 'etc/fps_RTConfig.txt')
    if outfilename is None:
        outfilename = 'img.fits'
    xys = ioutils.specific_instrument_entries_from_file(fpslayoutfilename, whichinstruments, [2,3])
    xys = opticsmath.unitize_xys(xys, fpsradius)
    imgutils.genimg(xys, width=8192, height=5210, outfilename=outfilename, bgIntensity=bgIntensity, bgGaussMean=bgGaussMean, bgGaussStdDev=bgGaussStdDev, superGaussPeak_min=superGaussPeak_min, superGaussPeak_max=superGaussPeak_max, superGaussAWAM_min=superGaussAWAM_min, superGaussAWAM_max=superGaussAWAM_max, superGaussHWHM_dmin=superGaussHWHM_dmin, superGaussHWHM_dmax=superGaussHWHM_dmax, coefs=coefs)
