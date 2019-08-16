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

def centroids_hexpeelbijected(filename, expectedhexs):
    """Finds centroids in a fits image and matches them to the expectedhexs
    
    if the amount of centroids founds is not the same as the amount of points in expectedhexs, None is returned

    Parameters:
        filename (str):
            the path to a fits image file (must have extension '.fits')
        expectedhexs (list):
            a list of interleaved hex circumradii and 'angle' positions
            these are the centroids that are expected in the image

    Returns:
        xys (list):
            a list of interleaved xy positions, the same length as expectedhexs
            these are the positions of the centroids in the image, ordered just like the expectedhexs
    """
    centroids = centroids(filename)
    if len(centroids) is not len(expectedhexs):
        return None
    centroidhexs = []
    for i in range(0,len(expectedhexs),2):
        r,a = opticsmath.transform_cart2hex_xy2ra(expectedhexs[i],expectedhexs[i+1])
        centroidhexs.append(r)
        centroidhexs.append(a)

def centroids(filename):
    pass