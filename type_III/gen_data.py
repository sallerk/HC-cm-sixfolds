# gen_data.py -- generate test data for the type III discriminant criterion.
# D = (a,c)_Q definite (a,c<0), basis 1,i,j,k with i^2=a, j^2=c, k=ij.
# T = 3x3 skew-Hermitian matrix over D (T_ji = -conj(T_ij), T_ii pure).
# Q = list of pure quaternions q (K = Q(q), q^2 = -Nrd(q)).
# Output: data.json (next to this script). Shared input: path1.py (main computation) and
# path2.py (independent check) read only this file.
import json, random, itertools, os
from cypari import pari

HERE = os.path.dirname(os.path.abspath(__file__))

random.seed(20261005)

def ram_set(a, c, bound=200):
    ps = [p for p in pari(f'primes(primepi({bound}))')]
    S = [int(p) for p in ps if int(pari(f'hilbert({a},{c},{p})')) == -1]
    return S

# find (a,c) for target discriminants (product of an odd number of primes => definite)
targets = [2, 3, 5, 7, 11, 13, 17, 19, 23, 30, 42, 66, 70, 78, 102, 105, 110, 130]
found = {}
for a in range(-1, -80, -1):
    for c in range(-1, -400, -1):
        prs = set(int(x) for x in pari(f'factor({2*a*c})[,1]~'))
        S = [p for p in prs if int(pari(f'hilbert({a},{c},{p})')) == -1]
        d = 1
        for p in S: d *= p
        if d in targets and d not in found:
            found[d] = (a, c)
    if len(found) == len(targets):
        break
print('found', found)
assert set(found) == set(targets)

def rq(R):  # random quaternion
    return [random.randint(-R, R) for _ in range(4)]

def rpure(R):
    while True:
        v = [0] + [random.randint(-R, R) for _ in range(3)]
        if any(v): return v

def conj(x): return [x[0], -x[1], -x[2], -x[3]]
def neg(x): return [-t for t in x]

def randT(kind):
    T = [[None]*3 for _ in range(3)]
    if kind == 'diag':
        for r in range(3):
            T[r][r] = rpure(4)
            for s in range(3):
                if s != r: T[r][s] = [0, 0, 0, 0]
    else:
        R = 3 if kind == 'dense' else 6
        for r in range(3):
            T[r][r] = rpure(R)
            for s in range(r+1, 3):
                x = rq(R)
                T[r][s] = x
                T[s][r] = neg(conj(x))
    return T

qs = []
for v in itertools.product(range(-2, 3), repeat=3):
    if not any(v): continue
    # keep one of +-v, primitive
    import math
    g = math.gcd(math.gcd(abs(v[0]), abs(v[1])), abs(v[2]))
    if g != 1: continue
    first = next(t for t in v if t != 0)
    if first < 0: continue
    qs.append([0] + list(v))

data = {'algebras': []}
for d in targets:
    a, c = found[d]
    Ts = []
    for t in range(24):
        kind = ['diag', 'dense', 'wide'][t % 3]
        Ts.append({'kind': kind, 'T': randT(kind)})
    data['algebras'].append({'disc': d, 'a': a, 'c': c, 'Ts': Ts, 'qs': qs})
json.dump(data, open(os.path.join(HERE, 'data.json'), 'w'))
print('wrote data.json with', len(targets), 'algebras,', len(qs), 'q per T, 24 T each')
