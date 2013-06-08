import math
from bisect import bisect_left
import numpy as np
from matplotlib.pyplot import plot
#from pandas import unique

    
def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

def get_primes(max):
    return primesfrom2to(max+1000)

def primesfrom2to(n):
    '''
    Original author Robert Hanks
    http://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
     Input n>=6, Returns a array of primes, 2 <= p < n
    '''
    sieve = np.ones(n/3 + (n%6==2), dtype=np.bool)
    for i in xrange(1,int(n**0.5)/3+1):
        if sieve[i]:
            k=3*i+1|1
            sieve[       k*k/3     ::2*k] = False
            sieve[k*(k-2*(i&1)+4)/3::2*k] = False
    return np.r_[2,3,((3*np.nonzero(sieve)[0][1:]+1)|1)]

def is_prime(p,primes):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(primes, p)
    if i != len(primes) and primes[i] == p:
        return True
    else:
        return  False

def is_p_admissible(H,p):
    '''
    This is the biggest bottleneck.
    '''
    count = np.bincount(H % p, minlength = p)
    return p != np.count_nonzero(count) 
#    return p != len(unique(H % p)) # this unique is from Pandas
    
def is_admissible(H,primes):
    admissible = True
    bad_primes = []
    for p in primes:
        p_admissible = is_p_admissible(H,p)
        if not p_admissible:
            bad_primes.append(p)
        admissible = admissible and p_admissible
#        if p > H[-1]:
        if p > H[-1] or not admissible: #better performance but less info
            break
    return admissible, bad_primes

def pi(primes, x):
    'Find index of leftmost item greater than or equal to x in list primes'
    i = bisect_left(primes, x)
    if i != len(primes):
        return i
    raise ValueError

def sift(H,sieve):
    '''
    This accomplishes the actual sieving.
    '''
    for p,s in sieve:
        cond = H % p != s
        H = H[cond]
#        if max(abs(H)) < p: # this actually hurts performance
#            break
    return H

def pick_best(H,k):
    '''
    If our set H is larger than k then we just pick the subset of cardinality k with minimal diameter
    '''
    diameters = []
    for i in range(len(H)-k+1):
        diameters.append(H[i+k-1] - H[i])
    i = np.argmin(diameters)
    return H[i:i+k]

def results(H,primes):
    print 'diameter:', H[-1]-H[0], 'length:', len(H), 'density:', 100*float(len(H))/(H[-1]-H[0]), 'admissible:', is_admissible(H,primes)

def see_response(H,start=2):
    signal = np.diff(H)
    ft = np.fft.fft(signal*np.kaiser(len(signal),14))
    mgft = np.abs(ft)
    xVals = np.fft.fftfreq(len(signal), d=1.0)
    plot(xVals[start:len(mgft)], mgft[start:],'.')
    
def see_response_fromfile(filename,start=2,sep='\n'):
    H = np.fromfile(filename,dtype='int64',sep=sep)
    print H
    see_response(H,start)
    return H