#import pdb
import numpy as np
from prime_utils import * # this module contains prime-related stuff such as pi() as well as checkign for admissibility
from math import pow, log, sqrt


#k = 34429
#N = 500000

k = 10719
N = 120000

# we always assume that elements of H are in increasing order
#H = np.arange(-N/2,N/2)
#H = np.arange(N) # is this worse than symmetric interval around 0?
H = np.arange(-45000,70000)

def optimal_schinzel(HH): 
    primes = get_primes(max(abs(HH)))
#    maxj = pi(primes,len(HH)/4)
    maxj = len(primes)
    results = []
    for i in range(6,15):
        H = HH
        for j in range(i+1,maxj,10):
            sieve =  [(p,1) for p in primes[:i]] + [(p,0) for p in primes[i:j+1]]
            sift(H,sieve)
            if is_admissible(H,primes)[0]:
                print i,j, 'diameter:', H[-1]-H[0], 'length:', len(H), 
                if len(H) > k:
                    print 'picking the best'
                    H = pick_best(H,k)
                res = (i,j,H[-1]-H[0])
                results.append( res )
    return results

def process(H):
    '''
    Just an example
    '''
    primes = get_primes(max(abs(H)))
#    sieve = [(p,1) for p in primes[:len(primes)/5]]
#    sieve = eratosthenes(primes,len(H)/(math.log(len(H))*4))
    y = 2
    x = len(H)
    z = x/(3.86*log(x)*pow(log(log(x)),pi(primes,y)))
#    print pi(primes,z)
    sieve = schinzel(primes,y,z)
#    sieve = schinzel(primes,20,2000)
    H = sift(H,sieve)
#    H = greedy(H,primes[pi(primes,z):pi(primes,z)+100])
    return H, primes

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
    z = pi(primes,z)+1
    res =  [(p,1) for p in primes[:y]] + [(p,0) for p in primes[y:z]]
    return res

def greedy(H,primes_list):
    '''
    Performs a greedy siftin of H. For each prime in primelist it finds the one
    that saves as much of H as possible if it is sifted by previous primes.
    '''
    for p in primes_list:
        counts = np.bincount(H % p, minlength = p)
        s = np.argmin(counts)
#        print 'greedy sifting with:',p,s,counts[s]
        H = sift(H,[ (p,s) ])
    return H

def greedy_greedy(H,B):
    '''
    Implementation of algorhitm by Andrew Sutherland
    See http://sbseminar.wordpress.com/2013/05/30/i-just-cant-resist-there-are-infinitely-many-pairs-of-primes-at-most-59470640-apart/#comment-23566
    '''
    primes = get_primes(max(abs(H)))
    piB = pi(primes,sqrt(B)) - 1
    sieve = [(2,1)] + [(p,0) for p in primes[1:piB+1]]
    print 'sifting mod 0 up to', primes[piB]
    H = sift(H,sieve)
    pik = pi(primes,k)
    print 'greedy sifting from', primes[piB], ' to', primes[pik]
    H = greedy(H,primes[piB:pik+1])
    return H,primes

#res = optimal_schinzel(H) # no gains

#H,primes = process(H)
#H = np.arange(-180568,207406)
H,primes = greedy_greedy(H,H[-1]-H[0])
results(H,primes) # this takes longest because is_admissible is expensive
#
if len(H) > k:
    H = pick_best(H,k)
    print 'best diameter of length', k, 'is:',H[-1]-H[0]
