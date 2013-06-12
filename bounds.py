# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 19:30:51 2013

@author: Vit Tucek

We rewrite all inequalities in form expr >= 0, solve for zero and take the next
integer value. If there are parameters involved we try to find the minimal root
possible.
"""

import numpy as np
import scipy as sp

from scipy.optimize import brentq
from numpy import ceil,floor,arange, log, sqrt

from prime_utils import get_primes,pi

ks = [3500000, 181000, 34429, 26024, 23283, 22949, 10719, 5000,4000,3000, 2000, 1000, 672]

primes = get_primes(max(ks))

# For lack of better options in python I've precomputed Moebius and totient 
# functions elsewhere. The file contains data from 1 to 500,000 and text
# files are also provided in case binary loading makes problems (which should
# not since these are only integers).
#mudata = np.fromfile('mobius.txt',dtype='int',sep='\n')
#eulerdata = np.fromfile('euler.txt',dtype='int',sep='\n')
mudata = np.fromfile('mobius.dat',dtype='int')
eulerdata = np.fromfile('euler.dat',dtype='int')

def mu(x):
    return mudata[int(x)-1]

def phi(x):
    return eulerdata[int(x)-1]
    
def muphi_sc(q):
    muq = float(mu(q))
    return muq*muq/phi(q)

muphi = np.vectorize(muphi_sc) # so we can apply muphi on arrays

def precompute_muomega():
    '''
    This routine precomputes the values that appear in second M-V inequality,
    i.e. it tabelates q -> mu(q)^2 * prod_{p divides q} 1/(p-1)
    '''
    res = []
    for q in range(len(mudata)):
        if q % 10000 == 0:
            print q
        piq = pi(primes,q)
        divcond = (q+1) % primes[:piq+1] == 0 # there is a shift due to the fact that Python's indices start at 0
        qprimes =  1.0/(primes[divcond]-1)
        p = qprimes.prod()
        muq = mudata[q]
        res.append(muq*muq*p)
    return np.array(res)

#muomegadata = precompute_muomega() # to comupte the values from mu by definition
#muomegadata = np.fromfile('muomega.txt',sep=',') # to load from text file
muomegadata = np.fromfile('muomega.dat') # to load from binary file which may not work depending on CPU arch

def muomega_sc(q):
    return muomegadata[int(q)-1]

muomega = np.vectorize(muomega_sc)

def approx_nr(k):
    return k*(log(k)+1)

#print 'Approximation of diameter H'
#for k in ks:
#    print approx_nr(k)

def brun_titchmarsh(k):
    def f(H):
        return 2*(H+1.0)/log(H+1.0) - k
    return int(ceil(brentq(f,1,k*k)))
    
def bombieri(k):
    def f(H):
        return 2*(H+1.0)/(log(H+1.0)-3) - k
    return int(ceil(brentq(f,1,k*k)))
#
#print 'Bounds from Brun-Titchmarsh'
#for k in ks:
#    print '|', '{0:,d}'.format(brun_titchmarsh(k))
  
def first_mv_f(k,Q):    
    '''
    This function returns the lowest H that satisfies the first M-V inequality
    The parameter Q > 1 is to be optimized over. Optimal Q is around sqrt(H).
    '''
    s = sum(muphi(arange(0,floor(Q+1))))
    return k*s - 1 - Q*Q
    
def first_mv(k):
    def obj_func(Q):
        return -first_mv_f(k,Q) # minus is here because scipy only provides minimizing routines
    optimalQ = sp.optimize.fminbound(obj_func,1,k*(log(k)+1))
    optimalH = first_mv_f(k,optimalQ)
    return int(ceil(optimalH))
    
#print 'Bounds from First M-V inequality for k <= 42860'
#for k in ks[2:]: # I have only first 500,000 of values of phi and mu
#    print '|', '{0:,d}'.format(first_mv(k))

def second_mv_f(H,k,z,c):
    '''
    We rewrite the second M-V inequality as -k*s(H,c,q,z) >= -1 and so we can
    seek for given c, k and Q such H that the value of this function is as 
    close to 1 as possible.
    .
    Parameter z is similar to Q in first_mv_f and c is a constant conjectured 
    to be 1 (c=1). The original choice of Montgomery-Vaughan was c = 1.5 which
    was later improved to c=sqrt(22)/pi and to c = 3.2/pi
    '''
    #    q = arange(1,floor(z+1),dtype='int')    
    q = arange(1,z,dtype='int')
    elements = muomegadata[q]/(H + 1.0 + c*q*z)
    return -k*elements.sum()

def second_mv_minH(k,z,c):
    '''
    We are looking for H that satisfies the M-V inequality as close as possible
    for given parameters k,z and c.
    '''
    def f(H):
        return second_mv_f(H,k,z,c)
    return sp.optimize.fminbound(f,k*(log(k)+1)*0.5,k*(log(k)+1)*1.5) #H should be somewhere in this interval right?

def second_mv(k,c):
    def obj_func(z):
        return -1.0*second_mv_minH(k,z,c) # minus is here because we want to maximize later and python provides only miniization
#    optimalz = sp.optimize.fminbound(obj_func,1,k*(log(k)+1)) # we maximize optimal H over parameter z
    az = sqrt(k*(log(k)+1)) # crude estimate of optimal z given by crude estimate on H
    print 'We are looking for optimal z around', az
    optimalz = sp.optimize.fminbound(obj_func,az*0.5,az*1.5) # we maximize optimal H over parameter z
    optimalH = second_mv_minH(k,optimalz,c)
    print 'optimal z:', optimalz
    print 'Are we close to -1?', second_mv_f(optimalH,k,optimalz,c)
    return int(ceil(optimalH))