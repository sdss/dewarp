# encoding: utf-8

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from .utils import get_config
from .utils import get_logger

from .utils import opticsmath
from .utils import ioutils
from .utils import imgutils

NAME = 'dewarp'

# Loads config
config = get_config(NAME)

# Inits the logging system. Only shell logging, and exception and warning catching.
# File logging can be started by calling log.start_file_logger(path).
log = get_logger(NAME)

__version__ = '0.0.0'

def dewarp(fpslayoutfilename='fps_RTConfig', fiducialsimgfile, fibersimgfile):
    """Computes the optical warp based on an image of fiducials and applys it to an image of fibers
    
    Parameters:
        fpslayoutfilename (str):
            path to a config file containing fiducial positions (probably has '.txt' extension)
        fiducialsimgfile (str):
            path to a fits image file containing illuminated fiducials (must have '.fits' extension)
        fibersimgfile
            path to a fits image file containing illuminated metrology fibers (must have '.fits' extension)

    Returns:
        xys (list):
            a list of interleaved xy coordinates of the dots found in the image, unwarped based on coefs
    """
    coefs = detectwarp(fpslayoutfilename, fiducialsimgfile)
    return applywarp(coefs, fibersimgfile)

def detectwarp(fpslayoutfilename='fps_RTConfig.txt', infilename='simulatedwarpedfiducials.fits'):
    """Computes the optical warp of an image based on fiducials
    
    Parameters:
        fpslayoutfilename (str):
            path to a config file containing fiducial positions (probably has '.txt' extension)
        infilename (str):
            path to a fits image file (must have '.fits' extension)

    Returns:
        coefs (warpcoefs):
            an object containing coefficents and functions to dewarp points, probably to be used in applywarp()
    """
    ideal_xys = ioutils.fiducial_xys_from_file(fpslayoutfilename)
    imgdata = imgutils.readimage(infilename)
    observed_xys = imgutils.centroids_hexpeelbijected(imgdata, ideal_xys)

def applywarp(coefs, infilename='simulatedwarpedfiducials.fits'):
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
    imgdata = imgutils.readimage(infilename)
    xys =imgutils.centroids(imgdata)
    return coefs.applytransform(xys)

def fakewarp(fpslayoutfilename='fps_RTConfig.txt', radius=350, whichinstrument='fiducial', outfilename='simulatedwarpedfiducials.fits'):
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
    xys = ioutils.specific_instrument_entries_from_file(fpslayoutfilename, lambda a: a.lower==whichinstrument.lower(), [2,3])
    xys = opticsmath.unitize_xys(xys, radius)
    imgutils.genimg(xys, outfilename=outfilename)
