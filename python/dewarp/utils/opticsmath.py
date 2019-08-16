# encoding: utf-8
#
# @Author:    Adam Mendenhall
# @Date:      August 8, 2019
# @Filename:  opticsmath.py
# @License:   BSD 3-Clause
# @Copyright: Adam Mendenhall
#

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import numpy
import math

class warpcoefs:
    
    def __init__(self, fromxys=None, toxys=None, maxerr=None):
        if fromxys is not None and toxys is not None and maxerr is not None:
            self.computeweights(fromxys, toxys, maxerr)
    
    #for key (x,y) we have bivariate polynomial for zernike, increasing row is x and increasing column is y
    zernike_xypolycoefs = {(0,0):numpy.array([[1.]]),(1,-1):numpy.array([[0.,1],[0,0]]),(1,1):numpy.array([[0.,0],[1,0]]),(2,2):numpy.array([[0.,0,-1],[0,0,0],[1,0,0]]),(2,-2):numpy.array([[0.,0,0],[0,2,0],[0,0,0]]),(2,0):numpy.array([[-1.,0,2],[0,0,0],[2,0,0]])}
    tx_xypolycoefs = numpy.array([[0.,0],[1,0]]) #a square 2d array of scalar coefficients
    ty_xypolycoefs = numpy.array([[0.,1],[0,0]]) #a square 2d array of scalar coefficients
    
    def get_zernike_xypolycoefs(self, n, m):
        """Gets the (n,m)th zernike polynomial (Noll numbering)
        
        returns a square 2d numpy array of scalar coefficients for the zernike polynomial of index n,m (Noll numbering), a bivariate polynomial in x and y
        
        coefficient order conforms to default numpy 2d polynomial evaluation:
        increasing row index    corresponds to increasing powers of x
        increasing column index corresponds to increasing powers of y
        
        each time this function is called on this object, the result is stored by the object
        subsequent calls use this memory when applicable to speed up computation (formulae are inherently recursive)
        
        domain of returned function is the unit disk (points with magnitude larger than 1 may evaluate to garbage)
        n%2 must be equal to m%2,             otherwise behavior is unspecified and unsafe
        n must be nonnegative,                otherwise behavior is unspecified and unsafe
        the magnitude of m must not exceed n, otherwise behavior is unspecified and unsafe
        
        Parameters:
            n,m (int):
                the indicies of the desired polynomial
        
        Returns:
            result (ndarray):
                a square 2d numpy array of coefficients
        """
        out = None
        if n<abs(m): #invalid zernikes are 0 (useful for general recursion)
            return numpy.zeros((1,1))
        elif (n,m) in self.zernike_xypolycoefs:
            out = self.zernike_xypolycoefs[(n,m)]
        elif n>abs(m): #general depth 2 recursion
            lastlast = self.get_zernike_xypolycoefs(n-4, m)
            last = self.get_zernike_xypolycoefs(n-2, m)
            out = numpy.zeros((max(lastlast.shape[0],last.shape[0]+2),max(lastlast.shape[1],last.shape[1]+2)))
            out[2:2+last.shape[0]:1,0:last.shape[1]:1] += 4*n*(n-1)*(n-2)/(n+m)/(n-m)/(n-2)*last
            out[0:last.shape[0]:1,2:2+last.shape[1]:1] += 4*n*(n-1)*(n-2)/(n+m)/(n-m)/(n-2)*last
            out[0:last.shape[0]:1,0:last.shape[1]:1] -= 2*(n-1)*(m**2+n*(n-2))/(n+m)/(n-m)/(n-2)*last
            out[0:lastlast.shape[0]:1,0:lastlast.shape[1]:1] -= n*(n+m-2)*(n-m-2)/(n+m)/(n-m)/(n-2)*lastlast
            self.zernike_xypolycoefs[(n,m)] = out
        else: #along the edge of abs(m)=n, still depth 2 recursion
            lastlast = self.get_zernike_xypolycoefs(n-2, int(m-2*m/abs(m)))
            last = self.get_zernike_xypolycoefs(n-1, int(m-m/abs(m)))
            out = numpy.zeros((max(lastlast.shape[0]+2,last.shape[0]+1),max(lastlast.shape[1]+2,last.shape[1])))
            out[0:lastlast.shape[0]:1,2:2+lastlast.shape[1]:1] -= lastlast
            out[2:2+lastlast.shape[0]:1,0:lastlast.shape[1]:1] -= lastlast
            out[1:1+last.shape[0]:1,0:last.shape[1]:1] += 2*last
            self.zernike_xypolycoefs[(n,m)] = out
        return out
    
    def get_orthophi_xypolycoefs(self, n, m):
        """Gets the (n,m)th phi polynomial (Noll numbering)
        
        returns a square 2d numpy array of scalar coefficients for the phi polynomial of index n,m (Noll numbering), a bivariate polynomial in x and y
        
        the set containing the divergences and curls of the phi polynomials is mutually orthogonal and a basis for all vector fields (in the unit disk)
        
        coefficient order conforms to default numpy 2d polynomial evaluation:
        increasing row index    corresponds to increasing powers of x
        increasing column index corresponds to increasing powers of y
        
        this function calls get_zernike_xypolycoefs
        
        domain of returned function is the unit disk (points with magnitude larger than 1 may evaluate to garbage)
        n%2 must be equal to m%2,             otherwise behavior is unspecified and unsafe
        n must be nonnegative,                otherwise behavior is unspecified and unsafe
        the magnitude of m must not exceed n, otherwise behavior is unspecified and unsafe
        
        Parameters:
            n,m (int):
                the indicies of the desired polynomial
        
        Returns:
            result (ndarray):
                a square 2d numpy array of coefficients
        """
        if n==abs(m):
            return self.get_zernike_xypolycoefs(n, m)
        else:
            out = numpy.copy(self.get_zernike_xypolycoefs(n, m))
            last = self.get_zernike_xypolycoefs(n-2, m)*math.sqrt(1+2/(n-1))
            out[:last.shape[0]:,:last.shape[1]:] -= last
            return out
    
    def extendbasis(self, degree):
        """Computes (all) the <phi, orthogonal grad/curl zernike, T and S> basis functions with n=degree
        
        returns an array of tuples, whose entries are square 2d numpy arrays of scalar coefficients for bivariate polynomials in x and y
        each tuple contains a polynomial for dx, then another polynomial for dy (theoretically vector argument, vector valued)
        for phi functions with equivalent divergence and curl (laplacian, occurs waow n=|m|), only one tuple is appended
        for phi functions with distinct divergence and curl (most of them local rotation is different from local scaling), two tuples are appended
        
        coefficient order conforms to default numpy 2d polynomial evaluation:
        increasing row index    corresponds to increasing powers of x
        increasing column index corresponds to increasing powers of y
        
        this function calls get_orthophi_xypolycoefs
        
        domain of returned functions is the unit disk (points with magnitude larger than 1 may evaluate to garbage)
        degree must be nonnegative, otherwise behavior is unspecified and unsafe
        
        Parameters:
            degree (int):
                the n index for each of the returned tuples
        
        Returns:
            result (list):
                a list of tuples of square 2d numpy arrays of coefficients
        """
        if degree<0:
            return None
        if degree==0:
            return [(numpy.zeros((1,1)),numpy.zeros((1,1)))]
        out = []
        for m in range(degree,-1,-2):
            if m==0: #no distinction between -m and +m (no sin nor cos factor)
                phi = self.get_orthophi_xypolycoefs(degree,0)
                pderx = polypder2d(phi,True)
                pdery = polypder2d(phi,False)
                out.append((pderx, pdery)) #div
                out.append((pdery,-pderx)) #curl
            else: #+m is distinct from -m, one phi is times sin, the other phi is times cos
                phicos = self.get_orthophi_xypolycoefs(degree,+m)
                phisin = self.get_orthophi_xypolycoefs(degree,-m)
                cospderx = polypder2d(phicos,True)
                cospdery = polypder2d(phicos,False)
                sinpderx = polypder2d(phisin,True)
                sinpdery = polypder2d(phisin,False)
                out.append((cospderx,cospdery)) #div
                out.append((sinpderx,sinpdery)) #div
                if m!=degree: #pair is not laplacian, divergence is distinct from curl
                    out.append((cospdery,-cospderx)) #curl
                    out.append((sinpdery,-sinpderx)) #curl
        return out
    
    def computetransform(self, fromxys, toxys, maxerr):
        """Computes the optimal (acceptable) coefficients to transform fromxys to toxys
        
        fromxys must be an interleaved array of NORMAL x y coordinates within the UNIT circle (magnitude no larger than 1)
        toxys must be an interleaved array of NORMAL x y coordinates within the UNIT circle (magnitude no larger than 1)
        
        calls extendbasis as many times as necessary until the optimal transform for a set of basis functions yields no error larger than maxerr
        all errors are measured measured along the x axis, then the y axis (not radially)
        
        if numerical instability begins occurring, stops computing and stores best so far (some errors may be larger than maxerr)
        
        applytransform(fromxys) should yield about toxys
        
        Parameters:
            fromxys (list):
                a list of interleaved NORMAL x y coordinates within the UNIT circle, presumably the locations of the observed fiducials
            toxys (list):
                a list of interleaved NORMAL x y coordinates within the UNIT circle, presumably the locations of the actual fiducials in real life
            maxerr (float):
                barring numerical instability, no errors (x and y measured independently) of the reconstruction of toxys are larger than maxerr
        
        Returns:
            result (float):
                the average error (of all the x errors and y errors) of the final computed transform (-1 if numerical instability occurs)
        """
        err = None
        avgerr = None
        degree = 0
        errxys = toxys-fromxys
        eval_matrix = numpy.empty((len(fromxys),0))
        square_matrix = numpy.empty((0,0))
        basiselements = []
        while err is None or err>maxerr:
            degree += 1
            newbasiselements = self.extendbasis(degree)
            basiselements.extend(newbasiselements)
            n = len(square_matrix)
            eval_matrix = numpy.hstack((eval_matrix,numpy.empty((len(fromxys),len(newbasiselements)))))
            for i in range(len(newbasiselements)):
                eval_matrix[0::2,n+i] = numpy.polynomial.polynomial.polyval2d(fromxys[0::2], fromxys[1::2], newbasiselements[i][0])
                eval_matrix[1::2,n+i] = numpy.polynomial.polynomial.polyval2d(fromxys[0::2], fromxys[1::2], newbasiselements[i][1])
            square_matrix = numpy.hstack((square_matrix, numpy.empty((n,len(newbasiselements))))) #append columns to square matrix
            square_matrix = numpy.vstack((square_matrix, numpy.empty((len(newbasiselements),n+len(newbasiselements))))) #append longer rows to matrix to return to 
            for i in range(n+len(newbasiselements)): #the $ij$th entry of square_matrix is $\sum_{x\in\vec x}T_i(x)T_j(x)$
                for j in range(max(n,i),n+len(newbasiselements),1):
                    square_matrix[i][j] = numpy.sum(numpy.multiply(eval_matrix[:,i],eval_matrix[:,j]))
                    if i!=j:
                        square_matrix[j][i] = square_matrix[i][j]
            try:
                inv_matrix = numpy.linalg.inv(square_matrix)
            except numpy.linalg.LinAlgError:
                continue
            newweights = inv_matrix.dot(eval_matrix.T).dot(errxys)
            dewarpedxys = fromxys+eval_matrix.dot(newweights)
            newerr = numpy.max(numpy.absolute(dewarpedxys-toxys))
            newavgerr = numpy.average(numpy.absolute(dewarpedxys-toxys))
            if avgerr is not None and newavgerr>avgerr: #numerical instability is occuring
                err = -1
                avgerr = -1
                for a in newbasiselements:
                    try: #this line might throw a ValueError about the truthiness of a nonsingleton array for some reason (maybe because of numpy vs python int types?)
                        basiselements.remove(a)
                    except ValueError:
                        pass
            else:
                weights = newweights
                err = newerr
                avgerr = newavgerr
        self.tx_xypolycoefs = numpy.zeros((max([c[0].shape[0] for c in basiselements]),max([c[0].shape[1] for c in basiselements])))
        self.ty_xypolycoefs = numpy.zeros((max([c[1].shape[0] for c in basiselements]),max([c[1].shape[1] for c in basiselements])))
        for i in range(len(basiselements)):
            c = basiselements[i]
            self.tx_xypolycoefs[:c[0].shape[0]:,:c[0].shape[1]:] += weights[i]*c[0]
            self.ty_xypolycoefs[:c[1].shape[0]:,:c[1].shape[1]:] += weights[i]*c[1]
        return avgerr
    
    def randomizetransform(self, degree, coefsrange=.01):
        """Assigns random coefficients to each basis function with n not more than degree
        
        calls extendbasis to retrive basis functions
        
        applytransform(fromxys) should yield a slightly distorted version of fromxys
        
        Parameters:
            degree (int):
                the maximum n of the basis functions included in warping
        """
        basiselements = []
        weights = []
        for n in range(1,degree+1,1):
            newbasiselements = self.extendbasis(n)
            basiselements.extend(newbasiselements)
            weights.extend(numpy.random.uniform(-coefsrange,coefsrange,(len(newbasiselements))))
        self.tx_xypolycoefs = numpy.zeros((max([c[0].shape[0] for c in basiselements]),max([c[0].shape[1] for c in basiselements])))
        self.ty_xypolycoefs = numpy.zeros((max([c[1].shape[0] for c in basiselements]),max([c[1].shape[1] for c in basiselements])))
        for i in range(len(basiselements)):
            c = basiselements[i]
            self.tx_xypolycoefs[:c[0].shape[0]:,:c[0].shape[1]:] += weights[i]*c[0]
            self.ty_xypolycoefs[:c[1].shape[0]:,:c[1].shape[1]:] += weights[i]*c[1]
    
    def applytransform(self, xys):
        """Hopefully dewarps xys
        
        xys must be an interleaved array of NORMAL x y coordinates within the UNIT circle (magnitude no larger than 1)
        
        applies an already computed transform to xys, presumably computed by computetransform(fromxys, toxys, maxerr)
        
        Parameters:
            xys (list):
                a list of interleaved NORMAL x y coordinates within the UNIT circle, presumably the locations of the observed fiducials
        
        Returns:
            result (list):
                a list of interleaved NORMAL x y coordinates within the UNIT circle, hopefully the locations of the actual fiducials in real life
        """
        outxys = [a for a in xys]
        outxys[0::2] += numpy.polynomial.polynomial.polyval2d(xys[0::2],xys[1::2],self.tx_xypolycoefs)
        outxys[1::2] += numpy.polynomial.polynomial.polyval2d(xys[0::2],xys[1::2],self.ty_xypolycoefs)
        return outxys

def polypder2d(poly, wrtX_otherwise_wrtY):
    """Takes the derivative of a bivariate polynomial
    
    takes the derivative of poly with respect to x or y, x if wrtX_otherwise_wrtY is True, or y if wrtX_otherwise_wrtY is False
    
    coefficient order conforms to default numpy 2d polynomial evaluation:
        increasing row index    corresponds to increasing powers of x
        increasing column index corresponds to increasing powers of y
    
    Paramters:
        poly (ndarray):
            a 2d numpy array of coefficients, a polynomial in x and y
    
    Returns:
        result (ndarray):
            a 2d numpy array of coefficients, the derivative of poly
    """
    if wrtX_otherwise_wrtY:
        return numpy.multiply(numpy.repeat(numpy.arange(1,poly.shape[0]),poly.shape[1]).reshape((poly.shape[0]-1,poly.shape[1]  )),poly[1::, ::])
    else:
        return numpy.multiply(numpy.tile(  numpy.arange(1,poly.shape[1]),poly.shape[0]).reshape((poly.shape[0],  poly.shape[1]-1)),poly[::, 1::])

def unitize_xys(xys, radius):
    """Scales xys to be in the unit disk
    
    if radius is None, overshrinks by 1.5 (you shouldn't let radius be None, you should know the proper radius to scale)
    
    Paramters:
        xys (list):
            a list of interleaved x y coordinates
        radius (float):
            the radius that all the xys reside in, if None then an appropriate raidus that encloses all points comfortably is chosen
    
    Returns:
        result (list):
            a list of interleaved NORMAL x y coordinates within the UNIT circle
    """
    if radius is None:
        radius = 1.5*math.sqrt(numpy.max(numpy.square(xys[0::2])+numpy.square(xys[1::2])))
    return numpy.divide(xys, radius)

def transform_cart2hex_xy2ra(x, y):
    """Computes the transform from cartesian to hexagonal

    Paramters:
        x,y (float):
            a coordinate in cartesian

    Returns:
        r,a (float):
            a tuple containing
            the radius of the hexagon (always nonnegative) (the circumradius, not the apothem) and 
            the 'angle' is in the range [0,1) and increases linearly when the cartesian point is moved linearly along an oriented hexagon (with top and bottom aligned parallel to x-axis)
    """
    r = None
    a = None
    root3 = math.sqrt(3)
    root3x = root3*x
    yonroot3 = y/root3
    ratio = x/y
    root3ratio = root3*ratio
    if x==0 and y==0:
        r = 0
        a = 0
    elif abs(root3x)<=abs(y):
        r = 2*yonroot3
        a = abs(1-root3ratio)/2
        if y<0:
            a += 4
        else:
            a += 1
    else:
        r = abs(x)+abs(yonroot3)
        if x*y>0:
            a = 2/abs(1+root3ratio)
            if x<0:
                a += 3
        else:
            a = -2/abs(1-root3ratio)
            if x<0:
                a += 3
            else:
                a += 6
    return (r,a/6)

def transform_hex2cart_ra2xy(r, a):
    """Computes the transform from hexagonal to cartesian

    Paramters:
        r,a (float):
            a coordinate in hexagonal, r is made to be (should be) nonnegative and a is made to be (should be) in the range [0,1)

    Returns:
        x,y (float):
            a tuple containing the xy cartesian point on an oriented hexagon (top and bottom parallel to x-axis)
            with radius (circumradius not apothem) r, 
            a unitary distance a along the perimeter starting from the right going counterclockwise
    """
    a = a%1
    r = abs(r)
    if r==0:
        return (0,0)
    amodhalf = a%.5
    sixth = 1/6
    third = 1/3
    if sixth<=amodhalf and amodhalf<=third:
        x = 3*(1-4*amodhalf)/2
        y = 1/2
    elif modhalf<=sixth:
        x = 1-3*amodhalf
        y = 3*amodhalf
    else:
        x = 1/3/sqrt(3)-3*amodhalf
        y = 3*(1/2-amodhalf)
    sgn = 1-2*floor(2*a)
    return (x*r*sgn,y*r*sqrt(3)*sgn)

def unittest_plotzernikes(degree=10):
    """Plots all zernikes with n index from 0 to degree and all m indicies"""
    coefs = warpcoefs()
    import matplotlib.pyplot as plt
    import time
    x = numpy.arange(-1,1.0001,0.005)
    y = numpy.arange(-1,1.0001,0.005)
    x,y = numpy.meshgrid(x,y)
    dont = x**2+y**2>1
    for n in range(0,degree):
        for m in range(-n,n+1,2):
            starttime = time.time()
            c = coefs.get_zernike_xypolycoefs(n, m)
            print('compute time: ',time.time()-starttime)
            z = numpy.polynomial.polynomial.polyval2d(x, y, c)
            z[dont] = 0
            plt.pcolormesh(x, y, z)
            plt.show()

def unittest_plotphifields(degree=10):
    """Plots all phi grad/curl basis fields with n index from 1 to degree and all m indicies"""
    coefs = warpcoefs()
    import matplotlib.pyplot as plt
    import time
    x = numpy.arange(-1,1.0001,0.05)
    y = numpy.arange(-1,1.0001,0.05)
    x,y = numpy.meshgrid(x,y)
    do = x**2+y**2<=1
    x = x[do]
    y = y[do]
    length = .05
    for n in range(1,degree):
        starttime = time.time()
        tuplescontainingdxthendypolys = coefs.extendbasis(n)
        print('time to compute basis coefficients of degree ',n,': ',time.time()-starttime)
        for c in tuplescontainingdxthendypolys:
            zx = numpy.polynomial.polynomial.polyval2d(x,y,c[0])
            zy = numpy.polynomial.polynomial.polyval2d(x,y,c[1])
            r = max(numpy.sqrt(zx**2+zy**2))
            zx = length*zx/r
            zy = length*zy/r
            plt.quiver(x,y,zx,zy,scale=1)
            plt.show()
