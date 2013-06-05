import math
from bisect import bisect_left
import numpy as np

def get_primes():
    return [2,3,5,7]

def add_prime(primes):
    prime = False
    n = primes[-1]
    while not prime:
        n = n+2;
        divisible = False
        p = primes[0]
        i = 0
        while not divisible and p <= math.sqrt(n):
            if n % p == 0:
                divisible = True
            i = i+1
            p = primes[i]
        if not divisible:
            prime = True
    primes.append(n)

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
#    return p != len(set(map(lambda x: x%p, H))) # std python
    return p != len(np.unique(H % p)) # numpy - much better performance
    

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