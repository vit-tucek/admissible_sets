#import pdb
import numpy as np
from prime_utils import * # this module contains prime-related stuff such as pi() as well as checkign for admissibility
from math import pow, log, sqrt


#k = 34429
#N = 500000

k = 10719
N = 120000

# we always assume that elements of H are in increasing order
H = np.arange(-N/2,N/2)
#H = np.arange(N) # is this worse than symmetric interval around 0?
#H = np.arange(-N/3,2*N/3)
#H = np.arange(-58486,53786)

#def two_chunk_greedy(H,primes):
#    '''
#        Same as greedy except it works in chunks of size 2
#    '''
#    for i in range(0,len(primes),2):
#        if i > len(primes):
#            break
#        c1 = np.bincount(H % primes[i], minlength = primes[i])
#        c2 = np.bincount(H % primes[i+1], minlength = primes[i+1])
#        m = 100000
#        coord = (0,0)
#        for r1 in range(primes[i]):
#            for r2 in range(primes[i+1]):
#                s = c1[r1] + c2[r2]
#                if s < m:
#                    m = s
#                    coord = (r1,r2)
#        H = sift(H,[ (primes[i],coord[0]), (primes[i+1],coord[1]) ])
#	if is_admissible(H,primes[i+2:])[0]:
#         print 'admissible at prime number',i, ' which is', primes[i]
#         break
#    return H


def two_chunk_greedy(H,primes):
    '''
        Same as greedy except it works in chunks of size 2
    '''
    for i in range(0,len(primes),2):
        if i > len(primes):
            break
        m = 1000000000000000
        coord = (0,0)
        for r1 in range(primes[i]):
            for r2 in range(primes[i+1]):
                test = sift(H, [ (primes[i], r1), (primes[i+1], r2)])
                s = len(test)
                if s < m:
                    m = s
                    coord = (r1,r2)
        H = sift(H,[ (primes[i],coord[0]), (primes[i+1],coord[1]) ])
	if is_admissible(H,primes[i+2:])[0]:
         print 'admissible at prime number',i, ' which is', primes[i]
         break
    return H

def greedy(H,primes_list):
    '''
    Performs a greedy siftin of H. For each prime in primelist it finds the one
    that saves as much of H as possible if it is sifted by previous primes.
    '''
    for i in range(len(primes_list)):
        p = primes_list[i]
        counts = np.bincount(H % p, minlength = p)
        s = np.argmin(counts)
#        print 'greedy sifting with:',p,s,counts[s]
        H = sift(H,[ (p,s) ])
	if is_admissible(H,primes_list[i+1:])[0]:
         print 'admissible at prime number',i, ' which is', primes_list[i]
         break
    return H

print 'symmetric'
primes = get_primes(max(abs(H)))
H = two_chunk_greedy(H,primes[:len(primes)/16])
results(H,primes) # this takes longest because is_admissible is expensive

H = np.arange(N)
print 'positive'
primes = get_primes(max(abs(H)))
H = two_chunk_greedy(H,primes[:len(primes)/16])
results(H,primes) # this takes longest because is_admissible is expensive

if len(H) > k:
    H = pick_best(H,k)
    print 'best diameter of length', k, 'is:',H[-1]-H[0]
