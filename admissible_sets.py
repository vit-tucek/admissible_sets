#import pdb
import numpy as np
from prime_utils import *

N = 100000
# we always assume that elements of H are in increasing order
H = np.arange(-N/2,N/2) 
# H = np.arange(N) # is this worse than symmetric interval around 0?

def process(H):
    primes = get_primes()
    maxH = max(abs(H))
    while primes[-1] < maxH:
        add_prime(primes)
#    sieve = [(p,1) for p in primes[:len(primes)/5]]
    sieve = eratosthenes(primes,math.log(N)*math.sqrt(N))
    H = sift(H,sieve)
    return H, primes
    
def sift(H,sieve):
    for p,s in sieve:
        cond = H % p != s
        H = H.compress(cond)
#        if max(abs(H)) < p: # this actually hurts performance
#            break
    return H
    

def eratosthenes(primes,cutoff):
    cutoff = pi(primes,cutoff) # returns the index of the first prime > N/2
    return [(p,0) for p in primes[:cutoff]]
    
def schinzel(primes,y,z):
    y = pi(primes,y)
    z = pi(primes,z+1)
    return [(p,1) for p in primes[:y]] + [(p,0) for p in primes[y:z+1]]
    
H,primes = process(H)
results(H,primes)