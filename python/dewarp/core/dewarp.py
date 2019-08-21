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

def dewarp(fiducialsimgfile, fibersimgfile, fpslayoutfilename=None):
    """Computes the optical warp based on an image of fiducials and applys it to an image of fibers
    
    Parameters:
        fiducialsimgfile (str):
            path to a fits image file containing illuminated fiducials (must have '.fits' extension)
        fibersimgfile
            path to a fits image file containing illuminated metrology fibers (must have '.fits' extension)
        fpslayoutfilename (str):
            path to a config file containing fiducial positions (probably has '.txt' extension)

    Returns:
        xys (list):
            a list of interleaved xy coordinates of the dots found in the image, unwarped based on coefs
    """
    if fpslayoutfilename is None:
        fpslayoutfilename = pkg_resources.resource_filename('dewarp', 'etc/fps_RTConfig.txt')
    coefs = detectwarp(fpslayoutfilename, radius, fiducialsimgfile)
    return applywarp(coefs, fibersimgfile)

def detectwarp(fpslayoutfilename=None, fiducialradius=350, infilename=None, imageradius=2500):
    """Computes the optical warp of an image based on fiducials
    
    Parameters:
        fpslayoutfilename (str):
            path to a config file containing fiducial positions (probably has '.txt' extension)
        fiducialradius (float):
            the radius (nonnegative) that all the dots are within, in whatever units the file would like (probably mm)
        infilename (str):
            path to a fits image file (must have '.fits' extension)
        imageradius (float):
            the radius (nonnegative) that all the dots in the image are within, in pixels

    Returns:
        coefs (warpcoefs):
            an object containing coefficents and functions to dewarp points, probably to be used in applywarp()
    """
    if fpslayoutfilename is None:
        fpslayoutfilename = pkg_resources.resource_filename('dewarp', 'etc/fps_RTConfig.txt')
    if infilename is None:
        infilename = pkg_resources.resource_filename('dewarp', 'etc/simulatedwarpedfiducials.fits')
    ideal_xys = ioutils.fiducial_xys_from_file(fpslayoutfilename)
    ideal_xys = opticsmath.unitize_xys(ideal_xys, fiducialradius)
    imgdata = imgutils.readimage(infilename)
    observed_xys = imgutils.centroids_hexpeelbijected(imgdata, ideal_xys)
    observed_xys = opticsmath.unitize_xys(observed_xys, imageradius)
    return warpcoefs(observed_xys, ideal_xys, .01)

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
    imgdata = imgutils.readimage(infilename)
    xys = imgutils.centroids(imgdata)
    return coefs.applytransform(xys)

def fakewarp(fpslayoutfilename=None, radius=350, whichinstrument='fiducial', outfilename=None):
    """Generates a warped image based on a configuration file
    
    Parameters:
        fpslayoutfilename (str):
            the path to a config file containing positions of instruments (probably has '.txt' extension)
        radius (float):
            the radius (nonnegative) that all the dots are within, in whatever units the file would like (probably mm)
        whichinstrument (str):
            the (case ignorant) name of the instruments to draw
        outfilename (str):
            the path to a (possibly nonexistent, overwrites of exists) fits image file (must have '.fits' extension)
    """
    if fpslayoutfilename is None:
        fpslayoutfilename = pkg_resources.resource_filename('dewarp', 'etc/fps_RTConfig.txt')
    if infilename is None:
        infilename = pkg_resources.resource_filename('dewarp', 'etc/simulatedwarpedfiducials.fits')
    xys = ioutils.specific_instrument_entries_from_file(fpslayoutfilename, lambda a: a.lower==whichinstrument.lower(), [2,3])
    xys = opticsmath.unitize_xys(xys, radius)
    imgutils.genimg(xys, outfilename=outfilename)