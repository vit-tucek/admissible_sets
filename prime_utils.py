import math
from bisect import bisect_left
import numpy as np

def get_primes(max):
    return primesfrom2to(max)

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

def add_primes(N,primes):
    '''
    adds N primes
    '''
    while len(primes) < N:
        add_prime(primes)
    return primes

def is_p_admissible(H,p):
    '''
    This is the biggest bottleneck.
    The following code is copied from numpy.unique with a change from quicksort
    to mergesort since I found the latter to perform slightly better.
    '''
    r = H % p
    r.sort(kind='mergesort')
    return p != len(r[np.concatenate( ([True], r[1:] != r[:-1]) )])
#    return p != len(np.unique(H % p))

def is_admissible(H,primes):
    admissible = True
    bad_primes = []
    for p in primes:
        if not is_p_admissible(H,p):
            bad_primes.append(p)
        admissible = admissible and is_p_admissible(H,p)
#        if p > H[-1]:
        if p > H[-1] or not admissible: #better performance but less info
            break
    return admissible, bad_primes

def pi(a, x):
    'Find index of leftmost item greater than or equal to x in list a'
    i = bisect_left(a, x)
    if i != len(a):
        return i
    raise ValueError

def results(H,primes):
    print 'diameter:', H[-1]-H[0], 'length:', len(H), 'admissible:', is_admissible(H,primes)[0]