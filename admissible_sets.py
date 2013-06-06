#import pdb
import numpy as np
from prime_utils import * # this module contains prime-related stuff such as pi() as well as checkign for admissibility
from math import pow, log, sqrt

#N = 10000
k = 34429
N = 500000
# we always assume that elements of H are in increasing order
H = np.arange(-N/2,N/2)
#H = np.arange(N) # is this worse than symmetric interval around 0?

def process(H):
    '''
    In this function we populate the primes and pick a sieve.
    '''
    primes = get_primes(max(abs(H)))
#    sieve = [(p,1) for p in primes[:len(primes)/5]]
#    sieve = eratosthenes(primes,N/(math.log(N)*4))
#    y = 10
#    x = len(H)
#    z = x/(0.2*log(x)*pow(log(log(x)),pi(primes,y)))
#    sieve = schinzel(primes,y,z)
    sieve = schinzel(primes,20,2000)
    H = sift(H,sieve)
#    H = greedy(H,primes[pi(primes,z):pi(primes,z)+100])
    return H, primes

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

def eratosthenes(primes,cutoff):
    '''
    Eratosthenes sieve as described in [GR1998]
    '''
    cutoff = pi(primes,cutoff) # returns the index of the first prime > N/2
    return [(p,0) for p in primes[:cutoff]]

def schinzel(primes,y,z):
    '''
    Schinzel sieve as described in [GR1998]
    '''
    y = pi(primes,y)
#    print primes[:y]
    z = pi(primes,z+1)
    res =  [(p,1) for p in primes[:y]] + [(p,0) for p in primes[y:z+1]]
    return res

def greedy(H,primes_list):
    '''
    Performs a greedy siftin of H. For each prime in primelist it finds the one
    that saves as much of H as possible if it is sifted by previous primes.
    '''
    for p in primes_list:
        counts = count(H,p)
        s = np.argmin(counts)
        print 'greedy sifting with:',p,s,counts[s]
        H = sift(H,[ (p,s) ])
    return H

def count(H,p): # expensive!
    counts = []
    H = H % p
    for r in range(p):
        cond = H == r
        l = len(H[cond])
        counts.append(l)
        if l == 0: # there is a modulus that is not contained in H -> profit
            break
    return counts

def pick_best(H,k):
    '''
    If our set H is larger than we need we just pick the subset with minimal diameter
    '''
    diameters = []
    for i in range(len(H)-k+1):
        diameters.append(H[i+k-1] - H[i])
    i = np.argmin(diameters)
    return H[i:i+k]

def greedy_greedy(H,B):
    primes = get_primes(max(abs(H)))
    piB = pi(primes,sqrt(B))
    sieve = [(2,1)] + [(p,0) for p in primes[1:piB]]
    H = greedy(H,primes[piB:k])
    H = sift(H,sieve)
    return H,primes


#H,primes = process(H)
H = np.arange(-185662,202456)
H,primes = greedy_greedy(H,4739224)
results(H,primes) # this takes longest because is_admissible is expensive

if len(H) > k:
    H = pick_best(H,k)
    print 'best diameter of length', k, 'is:',H[-1]-H[0]