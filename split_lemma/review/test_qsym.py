# test_qsym.py -- cross-check of the own Hilbert symbol code (qsym.hilbert) against
#  (i) PARI's hilbert(a,b,p) and (ii) exhaustive search for primitive zeros mod p^k
#  (theory-free).  Also checks the product formula.
import random, sys, time
from cypari import pari
from qsym import hilbert, hilbert_bruteforce, prime_factors, squarefree_part

random.seed(1)
t0 = time.time()
n_pari = n_bf = n_pf = 0
bad = []
# (i) PARI comparison, many pairs, all relevant primes + infinity
for trial in range(4000):
    a = random.choice([1, -1]) * random.randint(1, 10 ** random.randint(1, 8))
    b = random.choice([1, -1]) * random.randint(1, 10 ** random.randint(1, 6))
    ps = sorted(set([2, 3, 5, 7]) | set(prime_factors(a)) | set(prime_factors(b)))
    prod = hilbert(a, b, 0)
    for p in ps:
        h = hilbert(a, b, p)
        hp = int(pari.hilbert(a, b, p))
        n_pari += 1
        if h != hp:
            bad.append(('pari', a, b, p, h, hp))
        prod *= h
    if pari.hilbert(a, b, 0) != hilbert(a, b, 0):
        bad.append(('pari-inf', a, b))
    n_pf += 1
    if prod != 1:
        bad.append(('product formula', a, b))
# (ii) brute force, small primes, squarefree representatives
for trial in range(1500):
    a = random.choice([1, -1]) * random.randint(1, 400)
    b = random.choice([1, -1]) * random.randint(1, 400)
    for p in [2, 3, 5, 7, 11, 13]:
        h = hilbert(a, b, p)
        k = 5 if p == 2 else 2
        hb = hilbert_bruteforce(a, b, p, k)
        hb2 = hilbert_bruteforce(a, b, p, k + 1)
        n_bf += 1
        if not (h == hb == hb2):
            bad.append(('bruteforce', a, b, p, h, hb, hb2))
print(f"pari comparisons: {n_pari}, product-formula checks: {n_pf}, brute-force comparisons: {n_bf}")
print(f"disagreements: {len(bad)}")
for x in bad[:20]:
    print(x)
print(f"time {time.time()-t0:.1f}s")
