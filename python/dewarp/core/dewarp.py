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
from dewarp.core.exceptions import *
from dewarp import log
import numpy
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

def detectwarp(fpslayoutfilename=None, fpsradius=350, infilename=None, maskfilename=None, readmask=False, genmask=False, maskdilation=200):
    """Computes the optical warp of an image based on fiducials

    readmask overrides genmask -- if both are specified, maskfilename is first read, used for computation, then a new mask is written to maskfilename
    
    Parameters:
        fpslayoutfilename (str):
            path to a config file containing fiducial positions (probably has '.txt' extension)
        fpsradius (float):
            the radius (nonnegative) that all the dots are within, in whatever units the file would like (probably mm)
        infilename (str):
            path to a fits image file (must have '.fits' extension)
        maskfilename (str):
            path to a fits image file (must have '.fits' extension), assumed to be binary (1 is don't consider 0 is do)
        readmask (bool):
            whether or not to read maskfilename for a mask file and use it
        genmask (bool):
            if true, creates a new mask (if none is provided) or modifies the supplied mask based on found dots
        maskdilation (float):
            the radius around the dots to generate a new mask

    Returns:
        coefs (warpcoefs):
            an object containing coefficents and functions to dewarp points, probably to be used in applywarp()
    """
    if fpslayoutfilename is None:
        fpslayoutfilename = pkg_resources.resource_filename('dewarp', 'etc/fps_RTConfig.txt')
    if infilename is None:
        infilename = pkg_resources.resource_filename('dewarp', 'etc/simulatedwarpedfiducials.fits')
    ideal_xys = ioutils.fiducial_xys_from_file(fpslayoutfilename)
    log.info('got configuration fiducials from '+fpslayoutfilename+', there are '+str(int(len(ideal_xys)/2))+' of them')
    ideal_xys = opticsmath.unitize_xys(ideal_xys, fpsradius)
    imgdata = imgutils.readimage(infilename)
    mask = None
    if readmask:
        if maskfilename is not None:
            mask = imgutils.readimage(maskfilename)
        else:
            log.warn('detectwarp() called and readmask=True but maskfilename is None (where to read mask?)')
    elif genmask:
        if maskfilename is not None:
            log.info('creating mask file at '+maskfilename)
            mask = imgutils.genmask(ideal_xys, imgdata.shape[0], imgdata.shape[1], maskdilation, maskfilename)
        else:
            log.warn('detectwarp() called and genmask=True but maskfilename is None (where to put mask?)')
    if mask is not None and (mask.shape[0]!=imgdata.shape[0] or mask.shape[1]!=imgdata.shape[1]):
        log.warn('detectwarp() called and we have a mask with dims '+str(mask.shape[0])+'x'+str(mask.shape[1])+' which differs from the image '+infilename+' with dims '+str(imgdata.shape[0])+'x'+str(imgdata.shape[1]))
        mask = None
    if mask is None:
        log.debug('finding points unmasked')
    else:
        log.debug('finding points with mask')
    observed_xys = imgutils.centroids(imgdata, mask)
    log.debug('found '+str(int(len(observed_xys)/2))+' points in '+infilename)
    if len(observed_xys)!=len(ideal_xys):
        if mask is not None:
            log.warn('ideally there are '+str(int(len(ideal_xys)/2))+' dots but found '+str(int(len(observed_xys)/2))+', trying again with no mask')
            observed_xys = imgutils.centroids(imgdata, None)
            if len(observed_xys)!=len(ideal_xys):
                raise DewarpCandFindProperDotsError('tried without mask, expected '+str(int(len(ideal_xys)/2))+', found '+str(int(len(observed_xys)/2)))
            else:
                log.warn('succeeded with no mask')
        raise DewarpCandFindProperDotsError('tried without mask, expected '+str(int(len(ideal_xys)/2))+', found '+str(int(len(observed_xys)/2)))
    orig_observed_xys = [a for a in observed_xys]
    observed_xys = opticsmath.unitize_xys(observed_xys, min(imgdata.shape[0],imgdata.shape[1])/2)
    observed_xys = opticsmath.center_xys(observed_xys, -.5, -.5)
    observed_xys = opticsmath.sort_closest(observed_xys, ideal_xys)
    log.debug('sorted and matched ideal to observed, computing coefficients')
    coefs = opticsmath.warpcoefs(observed_xys, ideal_xys, .01)
    log.debug('done computing coefficents')
    if genmask:
        if maskfilename is not None:
            log.info('creating mask file at '+maskfilename)
            imgutils.genmask(orig_observed_xys, imgdata.shape[0], imgdata.shape[1], maskdilation, maskfilename)
        else:
            log.warn('detectwarp() called and genmask=True but maskfilename is None (where to put mask?)')
    return coefs

def applywarp(coefs, infilename=None, maskfilename=None, readmask=False, genmask=False, maskdilation=200):
    """Applys optical warp to an image
    
    Parameters:
        coefs (warpcoefs):
            an object containing coefficients and functions to dewarp points, probably computed using detectwarp()
        infilename (str):
            path to a fits image file (must have '.fits' extension)
        maskfilename (str):
            path to a fits image file (must have '.fits' extension), assumed to be binary (1 is don't consider 0 is do)
        readmask (bool):
            whether or not to read maskfilename for a mask file and use it
        genmask (bool):
            if true, creates a new mask (if none is provided) or modifies the supplied mask based on found dots
        maskdilation (float):
            the radius around the dots to generate a new mask

    Returns:
        xys (list):
            a list of interleaved xy coordinates of the dots found in the image, unwarped based on coefs
    """
    if infilename is None:
        infilename = pkg_resources.resource_filename('dewarp', 'etc/simulatedwarpedfiducials.fits')
    imgdata = imgutils.readimage(infilename)
    mask = None
    if readmask:
        if maskfilename is not None:
            mask = imgutils.readimage(maskfilename)
        else:
            log.warn('applywarp() called and readmask=True but maskfilename is None (where to read mask?)')
    if mask is not None and (mask.shape[0]!=imgdata.shape[0] or mask.shape[1]!=imgdata.shape[1]):
        log.warn('applywarp() called and we have a mask with dims '+str(mask.shape[0])+'x'+str(mask.shape[1])+' which differs from the image '+infilename+' with with dims '+str(imgdata.shape[0])+'x'+str(imgdata.shape[1]))
        mask = None
    if mask is None:
        log.debug('finding points unmasked')
    else:
        log.debug('finding points with mask')
    xys = imgutils.centroids(imgdata, mask)
    log.debug('found '+str(int(len(xys)/2))+' points in '+infilename)
    orig_xys = [a for a in xys]
    xys = opticsmath.unitize_xys(xys, min(imgdata.shape[0],imgdata.shape[1])/2)
    xys = opticsmath.center_xys(xys, -.5, -.5)
    if genmask:
        if maskfilename is not None:
            log.info('creating mask file at '+maskfilename)
            imgutils.genmask(orig_xys, imgdata.shape[0], imgdata.shape[1], maskdilation, maskfilename)
        else:
            log.warn('applywarp() called and genmask=True but maskfilename is None (where to put mask?)')
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
    log.info('found '+str(int(len(xys)/2))+' dots in config file '+fpslayoutfilename)
    xys = opticsmath.unitize_xys(xys, fpsradius)
    log.info('generating ideal/fake image at '+outfilename)
    imgutils.genimg(xys, width=8192, height=5210, outfilename=outfilename, bgIntensity=bgIntensity, bgGaussMean=bgGaussMean, bgGaussStdDev=bgGaussStdDev, superGaussPeak_min=superGaussPeak_min, superGaussPeak_max=superGaussPeak_max, superGaussAWAM_min=superGaussAWAM_min, superGaussAWAM_max=superGaussAWAM_max, superGaussHWHM_dmin=superGaussHWHM_dmin, superGaussHWHM_dmax=superGaussHWHM_dmax, coefs=coefs)
