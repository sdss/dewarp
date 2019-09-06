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
from .utils import get_config

if __name__ == '__main__':
    config = get_config('dewarp')
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='computes, applies and/or fakes warp of fiducial/metrology fiber dots')
    
    parser.add_argument('-ofid', '--observefiducialimage', type=str, help='')
    #parser.add_argument('-mfid', '--maskfiducial', type=str, help='')
    parser.add_argument('-ofib', '--observefiberimage', type=str, help='')
    #parser.add_argument('-mfib', '--maskfiber', type=str, help='')
    parser.add_argument('-or', '--observeradius', type=float, help='')
    parser.add_argument('-t', '--truelayout', type=str, help='')
    parser.add_argument('-tr', '--trueradius', type=float, help='')
    parser.add_argument('-gi', '--generateimage', type=str, help='')
    
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='sets verbose mode')
    
    args = parser.parse_args()

    coefs = None
    if args.generateimage is not None and args.truelayout is not None:
        if args.trueradius is None:
            args.trueradius = 350
        if args.verbose:
            print('generating image from layout',args.truelayout,'and putting it into',args.generateimage)
        fakewarp(fpslayoutfilename=args.truelayout, radius=args.trueradius, whichinstrument='fiducial', outfilename=args.generateimage)
    elif args.observefiducialimage is not None and args.truelayout is not None:
        if args.trueradius is None:
            args.trueradius = 350
        if args.observeradius is None:
            args.observeradius = 2500
        if args.verbose:
            print('detecting warp in',args.observefiducialimage,'with reference to config file',args.truelayout)
        coefs = detectwarp(args.truelayout, args.trueradius, args.observefiducialimage, args.observeradius)

    if coefs is not None:
        if args.observefiberimage is not None:
            hopefullytruexys = applywarp(coefs, args.observefiberimage, args.observeradius)
            print('dewarped')
            print(hopefullytruexys[0::2])
            print(hopefullytruexys[1::2])
