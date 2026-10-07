# forced_nontrivial.py -- reproduces forced_nontrivial.out from split_witnesses.json (output of verify_split.py).
# Counts the forced primes recorded by verify_split.py (p | 2 d d_F such that -d is a square in F_v for
# every v | p) and, among them, those where -d is NOT a square in Q_p, i.e. K_p = Q_p(sqrt(-d)) is a field
# contained in every completion F_v (the non-trivial case of the forced-prime argument).
# The count was first produced by an inline command; this file contains the same code.
import json, os
from cypari import pari

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'split_witnesses.json')))
Q = pari.nfinit(pari('x'))
nt = 0; tot = 0
for c in D:
    for f in c['forced']:
        if f['forced']:
            tot += 1
            sq = int(pari.nfislocalpower(Q, pari.idealprimedec(Q, f['p'])[0], -c['d'], 2))
            if not sq: nt += 1
print('forced primes', tot, ' non-trivial (K_p a field inside every F_v):', nt)
