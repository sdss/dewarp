# encoding: utf-8
#
# @Author:    Adam Mendenhall
# @Date:      August 21, 2019
# @Filename:  __main__.py
# @License:   BSD 3-Clause
# @Copyright: Adam Mendenhall
#

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import argparse
import os
import sys

from .core.dewarp import *
from . import log
from .utils.configuration import get_config
from .utils import opticsmath

if __name__ == '__main__':
    config = get_config('dewarp')
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='computes, applies and/or fakes warp of fiducial/metrology fiber dots')
    
    parser.add_argument('-ofid', '--observefiducialimage', type=str, help='path to fits image containing illuminated fiducials')
    parser.add_argument('-mfid', '--maskfiducial', type=str, help='')
    parser.add_argument('-ofib', '--observefiberimage', type=str, help='path to fits image containing illuminated fibers')
    parser.add_argument('-mfib', '--maskfiber', type=str, help='')

    parser.add_argument('-ifid', '--idealfiduciallayout', type=str, help='path to config file containing fiducials')
    parser.add_argument('-irad', '--idealradius', type=float, help='the radius of the focal plane system, all instruments in the config file must be within this radius')

    parser.add_argument('-iimg',        '--idealimagegenerate',                            dest='iimg',         type=str, help='path to fits image file (potentially nonexistent) to write generated image to')
    parser.add_argument('-iimgw',       '--idealimagewidth',                               dest='iimgw',        type=float, help='width of generated image', default=8192)
    parser.add_argument('-iimgh',       '--idealimageheight',                              dest='iimgh',        type=float, help='height of generated image', default=5120)
    parser.add_argument('-iimgbg',      '--idealimagebackground',                          dest='iimgbg',       type=float, help='background intensity of generated image', default=3)
    parser.add_argument('-iimgbgmean',  '--idealimagebackgroundgaussmean',                 dest='iimgbgmean',   type=float, help='mean for gaussian noise added to background of generated image', default=5)
    parser.add_argument('-iimgbgstddv', '--idealimagebackgroundgaussstddev',               dest='iimgbgstddv',  type=float, help='standard deviation for gaussian noise added to background of generated image', default=2)
    parser.add_argument('-iimgsgpl',    '--idealimagedotsupergausspeakglb',                dest='iimgsgpl',     type=float, help='lower bound for peaks of dots drawn on generated image', default=150)
    parser.add_argument('-iimgsgpu',    '--idealimagedotsupergausspeaklub',                dest='iimgsgpu',     type=float, help='upper bound for peaks of dots drawn on generated image', default=210)
    parser.add_argument('-iimgsgawaml', '--idealimagedotsupergaussalmostwidthalmostmaxglb',dest='iimgsgawaml',  type=float, help='lower bound for AlmostWidthAlmostMaximums (how far away is 255/256 of peak) of dots drawn on generated image', default=6)
    parser.add_argument('-iimgsgawamu', '--idealimagedotsupergaussalmostwidthalmostmaxlub',dest='iimgsgawamu',  type=float, help='upper bound for AlmostWidthAlmostMaximums (how far away is 255/256 of peak) of dots drawn on generated image', default=10)
    parser.add_argument('-iimgsghwhml', '--idealimagedotsupergaussahlfwidthhalfmaxglb',    dest='iimgsghwhml',  type=float, help='lower bound for HalfWidthHalfMaximums (how far away is half of peak) of dots drawn on generated image', default=2)
    parser.add_argument('-iimgsghwhmu', '--idealimagedotsupergaussahlfwidthhalfmaxlub',    dest='iimgsghwhmu',  type=float, help='upper bound for HalfWidthHalfMaximums (how far away is half of peak) of dots drawn on generated image', default=5)
    parser.add_argument('-iimgdegu',    '--idealimagewarpcoefsdegreelub',                  dest='iimgdegu',     type=float, help='all basis functions of this degree and lower are randomized to warp the generated image', default=4)

    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='sets verbose mode')
    parser.add_argument('-d', '--debug', action='store_true', default=False, help='sets debug mode')
    
    args = parser.parse_args()
    
    if args.debug:
        log.setLevel(-1)
    elif args.verbose:
        log.setLevel(10)

    coefs = None
    if args.iimg is not None and args.idealfiduciallayout is not None:
        if args.idealradius is None:
            args.idealradius = 350
        coefs = opticsmath.warpcoefs()
        log.info('randomizing warp coefficients, up to degree '+args.iimgdegu)
        coefs.randomizetransform(args.iimgdegu)
        log.info('generating image from layout '+args.idealfiduciallayout+' and putting it into '+args.iimg)
        fakewarp(fpslayoutfilename=args.idealfiduciallayout, fpsradius=args.idealradius, whichinstruments=lambda a:a.lower()=='fiducial', outfilename=args.iimg, width=args.iimgw, height=args.iimgh, bgIntensity=args.iimgbg, bgGaussMean=args.iimgbgmean, bgGaussStdDev=args.iimgbgstddv, superGaussPeak_min=args.iimgsgpl, superGaussPeak_max=args.iimgsgpu, superGaussAWAM_min=args.iimgsgawaml, superGaussAWAM_max=args.iimgsgawamu, superGaussHWHM_dmin=args.iimgsghwhml, superGaussHWHM_dmax=args.iimgsghwhmu, coefs=coefs)
    elif args.observefiducialimage is not None and args.idealfiduciallayout is not None:
        if args.idealradius is None:
            args.idealradius = 350
        readmask = False
        genmask = False
        if args.maskfiducial is not None:
            genmask = True
            if os.path.exists(args.maskfiducial) and os.path.isfile(args.maskfiducial):
                readmask = True
        log.info('detecting warp in '+args.observefiducialimage+' with reference to config file '+args.idealfiduciallayout)
        coefs = detectwarp(args.idealfiduciallayout, args.idealradius, args.observefiducialimage, args.maskfiducial, readmask, genmask)

    if coefs is not None:
        if args.observefiberimage is not None:
            readmask = False
            genmask = False
            if args.maskfiber is not None:
                genmask = True
                if os.path.exists(args.maskfiber) and os.path.isfile(args.maskfiber):
                    readmask = True
            log.info('applying computed warp to '+args.observefiberimage)
            hopefullytruexys = applywarp(coefs, args.observefiberimage, args.maskfiber, readmask, genmask)
            log.info('done applying warp, displaying')
            import matplotlib.pyplot as p
            p.scatter(hopefullytruexys[0::2],hopefullytruexys[1::2])
            p.show()
