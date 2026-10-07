# check_lemma2.py -- test the proof step, not just the final criterion:
# Lemma 6.1 of the note (printed as "Lemma 2" below):
#   det_K(H_K) == b*c*Nrd(T)  mod Nm(K^*),  where q^2 = -b and c = j^2 for some j in D
#   anticommuting with q (so D = (-b,c)).  Uses the Gram determinants stored by path1 (non-diagonal T
#   included: the proof diagonalizes T, the computation does not).
# Reads path1_out.json, writes check_lemma2_out.txt (both next to this script).
import json, os
from fractions import Fraction as F
from path1 import nrd, mul, primes_of, hilb

HERE = os.path.dirname(os.path.abspath(__file__))

rows = json.load(open(os.path.join(HERE, 'path1_out.json')))
ok = bad = 0
for r in rows:
    a, c0 = r['a'], r['c']
    q = [F(t) for t in r['q']]
    b = F(r['b']); detH = F(r['detH']); N = F(r['NrdT'])
    # pure x trace-orthogonal to q: -a q1 x1 - c0 q2 x2 + a c0 q3 x3 = 0
    g = [-a*q[1], -c0*q[2], a*c0*q[3]]
    ker = [v for v in ([g[1], -g[0], F(0)], [g[2], F(0), -g[0]], [F(0), g[2], -g[1]]) if any(v)]
    x = [F(0)] + ker[0]
    assert mul(q, x, a, c0) == [-t for t in mul(x, q, a, c0)]
    c = -nrd(x, a, c0)            # x^2 = -Nrd(x)
    ratio = detH / (b*c*N)
    P = primes_of(b, ratio)
    isnorm = ratio > 0 and all(hilb(-b, ratio, p) == 1 for p in P)
    if isnorm: ok += 1
    else: bad += 1
print(f'Lemma 2 (det_K H == b c Nrd(T) mod Nm K^*): holds in {ok}/{ok+bad} cases')
open(os.path.join(HERE, 'check_lemma2_out.txt'), 'w').write(f'Lemma 2 holds in {ok}/{ok+bad} cases (path1 Gram determinants)\n')
