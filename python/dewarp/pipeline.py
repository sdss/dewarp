# encoding: utf-8
#
# @Author:    Adam Mendenhall
# @Date:      August 1, 2019
# @Filename:  pipeline.py
# @License:   BSD 3-Clause
# @Copyright: Adam Mendenhall
#

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from utils import opticsmath
from utils import ioutils

#detect warp
ideal = ioutils.fiducial_xys_from_file('fps_RTConfig.txt')
ideal = opticsmath.unitize_xys(ideal, None)
import random
fakeactual = [a+random.uniform(-.05,.05) for a in ideal]
coefs = opticsmath.warpcoefs()
print('average error: ',coefs.computetransform(fakeactual, ideal, .1))
dewarped = coefs.applytransform(fakeactual)

#visualize stuff
opticsmath.unittest_plotzernikes(4)
opticsmath.unittest_plotphifields(4)
