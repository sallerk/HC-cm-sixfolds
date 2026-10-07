# nonsplit_examples.py -- for every test algebra D, exhibit a diagonal skew-Hermitian T on D^3 whose
# Delta = -Nrd(T) has Q(sqrt Delta) NOT embedded in D (predicted: NO K in D split), and check ~49 K's
# with path1 and the first 8 with path2 as well.
# Reads data.json, writes nonsplit_examples_out.txt (both next to this script).
# Note: importing examples.py runs examples.py as well (its output is printed first and
# examples_out.json is rewritten with identical content).
import json, itertools, os
from fractions import Fraction as F
from path1 import nrd, nrd_matrix
from exists import gram_split, is_sq_Qp, disc_primes, construct_K
from examples import hyperbolic_path2, diagT

HERE = os.path.dirname(os.path.abspath(__file__))

data = json.load(open(os.path.join(HERE, 'data.json')))
pures = [[0] + list(v) for v in itertools.product(range(-3, 4), repeat=3) if any(v)]
pures.sort(key=lambda v: sum(abs(t) for t in v))
report = []
for A in data['algebras']:
    a, c, disc = A['a'], A['c'], A['disc']
    dp = disc_primes(a, c)
    found = None
    for l1, l2, l3 in itertools.product(pures[:40], repeat=3):
        T = diagT(l1, l2, l3)
        TF = [[[F(t) for t in e] for e in row] for row in T]
        N = nrd(TF[0][0], a, c) * nrd(TF[1][1], a, c) * nrd(TF[2][2], a, c)
        Delta = -N
        if any(is_sq_Qp(Delta, p) for p in dp):
            found = (T, TF, Delta); break
    T, TF, Delta = found
    assert nrd_matrix(TF, a, c) == -Delta
    n1 = n2 = 0; s1 = s2 = 0
    for qi, qv in enumerate(A['qs']):
        b, dH, sp = gram_split(TF, [F(t) for t in qv], a, c)
        n1 += 1; s1 += sp
        if qi < 8:
            n2 += 1; s2 += hyperbolic_path2(a, c, T, qv)
    line = (f'disc {disc} (a,c)=({a},{c}) T=diag{tuple(tuple(T[i][i][1:]) for i in range(3))} Delta={Delta} '
            f'Q(sqrt Delta) in D: False | path1: {s1}/{n1} K split | path2: {s2}/{n2} K split')
    print(line); report.append(line)
open(os.path.join(HERE, 'nonsplit_examples_out.txt'), 'w').write('\n'.join(report) + '\n')
