# 
# SDSS5 Software at UW
# Author:      Adam Mendenhall
# Creation:    26/6/2019
# 
# This program generates an image of fiducials.  Gaussian and Shot noise are introduced.
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

#functions#####
'''
returns a tuple containing the x positions and y positions (resp.) of the instruments in the file
only contains information for instruments which satisfy True==whichInstruments(name), where name is the name of the instrument as it appears in the file (e.g. 'BOSS' 'BA' or 'Fiducial')
assumes each instrument is on its own line
ignores (only) lines that begin with '#'
assumes attribute fields for each instrument are separated (only) by white space
assumes that the third and fourth field on a line are the x and y attributes of the instrument
assumes that the fifth field on a line is the name of the instrument
'''
def specificInstrumentsXYsFromFile(filenameWithExtension, whichInstruments):
    f = open(filenameWithExtension, 'r')
    xs = []
    ys = []
    for line in f:
        if len(line)>0 and line[0]!='#':
            attributes = line.split()
            if len(attributes)>4:
                if whichInstruments(attributes[4]):
                    try:
                        xs.append(float(attributes[2]))
                        ys.append(float(attributes[3]))
                    except ValueError:
                        pass
    return (xs,ys)

#using OSA/ANSI standard indicies, NOT Noll's indicies
def index_j2n(j):
    return int(math.floor(math.sqrt(2*j+.25)-.5))
def index_j2m(j):
    n = index_j2n(j)
    return int(2*j-n*(n+2))
def index_nm2j(n,m):
    return int(.5*(n*(n+2)+m))

#These are not strictly zernike coefficients, those for (n,m) are the same as those for (n,-m) instead of negatives of each other -- typically one would multiply by sgn(m*(2*(n%2)-1))
#When n<0 or m<0, outcome element is emtpy array
#When |m|>n, outcome element is emtpy array
#When n%2!=m%2, outcome element is emtpy array
#When n==m==0, outcome element is 1 (singleton array)
def zernike_coefs_list(nms):
    out = []
    for i in range(len(nms)):
        out.append(zernike_coefs(nms[i][0], nms[i][1]))
    return out

#These are not strictly zernike coefficients, those for (n,m) are the same as those for (n,-m) instead of negatives of each other -- typically one would multiply by sgn(m*(2*(n%2)-1))
#When n<0 or m<0, outcome is emtpy array
#When |m|>n, outcome is emtpy array
#When n%2!=m%2, outcome is emtpy array
#When n==m==0, outcome is 1 (singleton array)
def zernike_coefs(n, m):
    out = []
    if n<0 or n%2!=m%2 or abs(m)>n: #invalid
        return out
    if n==0:
        out.append(1)
        return out
    sgn = 1
    maxIdx = int((n-abs(m))/2)
    for k in range(maxIdx+1):
        out.append(combinations(n-k,k)*combinations(n-2*k,maxIdx-k)*sgn)
        sgn = -sgn
    return out

#When n<0 or m<0, outcome element is emtpy array
#When |m|>n, outcome element is emtpy array
#When n%2!=m%2, outcome element is emtpy array
#When n==m==0, outcome element is 0 (singleton array)
def orthoGradZernike_coefs_list(nms):
    out = []
    for i in range(len(nms)):
        out.append(orthoGradZernike_coefs(nms[i][0], nms[i][1]))
    return out

#If arguments are valid, returns a list with two elements: a list containing gradX coefficients and a list containing gradY coefficients
#When n<0 or m<0, outcome is emtpy array
#When |m|>n, outcome is emtpy array
#When n%2!=m%2, outcome is emtpy array
#When n==m==0, outcome is [[0],[0]]
def orthoGradZernike_coefs(n, m):
    out = []
    if n<0 or n%2!=m%2 or abs(m)>n: #invalid
        return out
    #the following logic is correct but stupid (meaning nonmathematical, inelegant) TODO: please the math gods
    if n==0: #piston
        out.append([0])
        out.append([0])
        return out
    if n==1: #tilt
        if m==1:
            out.append([zernike_coefs(0,0), 0])
        else:
            out.append([0, zernike_coefs(0,0)])
    if m==0:
        x_minus = [0]
    else:
        x_minus = zernike_coefs(n-1,m-1)
    y_minus = zernike_coefs(n-1,-m-1)
    if m==-1:
        x_plus = [0]
    else:
        x_plus = zernike_coefs(n-1,m+1)
    if m==1 or m==0:
        y_plus = [0]
    else:
        y_plus = [-a for a in zernike_coefs(n-1,-m+1)] #always negative (when nonzero, doesn't matter when zero)
    if m%n==0:
        c = math.sqrt(.5)
    else:
        c = .5
    #up to here the lists are of type integer, TODO: perhaps preserve that and multiply by fudge factor later
    x_plus = [a*c for a in x_plus]
    x_minus = [a*c for a in x_minus]
    y_plus = [a*c for a in y_plus]
    y_minus = [a*c for a in y_minus]
    if m==1:
        root2 = math.sqrt(2)
        x_minus = [a*root2 for a in x_minus]
    if m==-1:
        root2 = math.sqrt(2)
        y_minus = [a*root2 for a in y_minus]
    #add the minus and plus lists, making sure to line up from the right (the right is the r^0 or r^1 coefficient)
    x_out = [0]*max(len(x_plus),len(x_minus))
    y_out = [0]*max(len(y_plus),len(y_minus))
    for i in range(len(x_plus)-1,-1,-1):
        x_out[i] += x_plus[i]
    for i in range(len(x_minus)-1,-1,-1):
        x_out[i] += x_minus[i]
    for i in range(len(y_plus)-1,-1,-1):
        y_out[i] += y_plus[i]
    for i in range(len(y_minus)-1,-1,-1):
        y_out[i] += y_minus[i]
    return [x_out, y_out]

#computes 'a choose b' or nCr(a,b) = a!/b!(a-b)!
def combinations(a, b):
    if b<0 or b>a:
        return 0
    if b==0 or b==a:
        return 1
    if b>a/2:
        return combinations(int(a),int(a-b))
    out = 1
    for x in range(int(a-b),int(a)):
        out *= (x+1)#not incuding b-a, including b-a+1 up through including a
    for x in range(2,b+1):
        out = int(out/x)#including 2 up through b
    return int(out)

def distort(x, y, coefs):
    a = numpy.arctan2(y, x)#safe (=0) in the case of x=y=0, results are in range [-pi,pi]
    r = numpy.square(x) + numpy.square(y)
    return (numpy.multiply(x,r),numpy.multiply(y,r))
#end functions#

elapsedTime = time.time()

#parameters and variables
width = 1920#8192
height = 1080#5120
fpslayoutfilename = 'fps_RTConfig.txt'
outfilename = 'img.fits'
bgIntensity = 10
bgGaussMean = 0 #in addition to bgIntensity, set to 0 to add and subtract equally likely and overall add no intensity
bgGaussStdDev = .0 #bigger means bigger spread (means more random), set to 0 to have no randomness
opticalDistortCoefs = [1] #each coefficient scales contribution the jth orthonormalized gradient Zernike polynomial, make the first few bigger and have them trail off, units are roughly in pixels
entireRadius = -1 #if positive, all fiducials are assumed to be within a circle of this radius, otherwise it is computed
units2Pixels = min(width,height)/2 #radius not diameter...



print('Loading fiducials from %s'%fpslayoutfilename)
maxR = 0
#allInstruments = specificInstrumentsXYsFromFile(fpslayoutfilename, lambda x: True)
fidXs, fidYs = specificInstrumentsXYsFromFile(fpslayoutfilename, lambda x: x.lower()=='fiducial')
fidGaussPeaks = numpy.full_like(fidXs, 200.) #TODO: randomize
fidGaussHWHMs = numpy.full_like(fidXs, 5.) #TODO: randomize TODO: make these units not pixels, but fiducial units
fidRelevantRadii = numpy.ceil(fidGaussHWHMs*numpy.sqrt(numpy.log2(fidGaussPeaks))).astype('uint32')
if entireRadius<=0:
    #circular margin of 30 units
    entireRadius = math.sqrt(numpy.nanmax(numpy.multiply(fidXs,fidXs)+numpy.multiply(fidYs,fidYs)))+30
print('We have %d fiducials within a radius (in fiducial units) of %f'%(len(fidXs),entireRadius))



print('Computing optical distortions, the first %d of them'%len(opticalDistortCoefs))

#take a coefficient for index j, which defines m,n
#we have Z_n^m(\rho)=\rho^n sum_{k=0}^{{n-m}\over2}\mathcal K(m,n,k)\rho^{-2k}
#we have S_n^m(\rho)=\hat i\stackrel{\hat i}{C_n^m}Z_{n_{\hat i}}^{m_{\hat i}}(\rho)+\hat j\stackrel{\hat j}{C_n^m}Z_{n_{\hat j}}^{m_{\hat j}}(\rho)
#given (x,y), we want to move (x,y) by sum()


print('Creating image of size %dx%d'%(width,height))
data = numpy.full((width, height), float(bgIntensity))



print('Drawing warped fiducials...')
if len(opticalDistortCoefs)>0:
    xcs = numpy.divide(fidXs,entireRadius)
    ycs = numpy.divide(fidYs,entireRadius)
    dx,dy = distort(xcs, ycs, opticalDistortCoefs)
    xcs = (numpy.multiply(xcs+dx,units2Pixels)+(width /2)).astype('int32')
    ycs = (numpy.multiply(ycs+dy,units2Pixels)+(height/2)).astype('int32')
else:
    xcs = (numpy.multiply(numpy.divide(fidXs,entireRadius),units2Pixels)+(width /2)).astype('int32')
    ycs = (numpy.multiply(numpy.divide(fidYs,entireRadius),units2Pixels)+(height/2)).astype('int32')
xmins = numpy.maximum(0,numpy.minimum(width -1,xcs-fidRelevantRadii)).astype('uint32')
xmaxs = numpy.maximum(0,numpy.minimum(width -1,xcs+fidRelevantRadii)).astype('uint32')
ymins = numpy.maximum(0,numpy.minimum(height-1,ycs-fidRelevantRadii)).astype('uint32')
ymaxs = numpy.maximum(0,numpy.minimum(height-1,ycs+fidRelevantRadii)).astype('uint32')
#perhaps there is a fast numpy way to add each of these fiducial dot blocks into the image.  The following commented lines are not it TODO: settle for for or strive for greatness
#dot = numpy.reshape([numpy.multiply(numpy.power(2,-(numpy.square(x+xmins-xcs)+numpy.square(y+ymins-ycs)).astype('float32').divide(numpy.square(fidGaussHWHMs))),fidGaussPeaks) for x,y in numpy.ndindex((xmaxs-xmins, ymaxs-ymins))], (xmaxs-xmins, ymaxs-ymins))
#data[xmins:xmaxs, ymins:ymaxs] += dot
for i in range(0,len(fidXs)):
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