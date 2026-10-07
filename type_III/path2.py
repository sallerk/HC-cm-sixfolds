# path2.py -- independent check of path1.py (no PARI; pure Python Fractions).
# D -> M_2(Q(r)), r^2 = c:  i -> [[0,a],[1,0]], j -> [[r,0],[0,-r]], k = ij -> [[0,-a r],[r,0]].
# Riemann form E(x,y) = Trd(sum x_m T_mn conj(y_n)) as a 12x12 rational matrix on the Q-basis
# {beta*eps_m}, beta in {1,i,j,k}. Q = matrix of left multiplication by q.
# Trace form of the K-Hermitian form: S(x,y) = -E(x, q y) (symmetric). The K-Hermitian form h
# with E = Tr(q^{-1} h) is split (hyperbolic) iff its trace form S is hyperbolic (Jacobson).
# Hyperbolicity: signature (6,6), det a square, Hasse invariants = those of <1,-1>^6.
# Delta = -Nrd(T) computed with the Q(sqrt c) embedding (6x6 det over Q(r), own elimination).
# Own Hilbert symbol implementation. Reads data.json; writes path2_out.json (both next to this script).
import json, sys, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- arithmetic in Q(r), r^2 = c : pairs (x,y) = x + y r
class QR:
    def __init__(s, c): s.c = c
    def add(s, u, v): return (u[0]+v[0], u[1]+v[1])
    def sub(s, u, v): return (u[0]-v[0], u[1]-v[1])
    def mul(s, u, v): return (u[0]*v[0] + s.c*u[1]*v[1], u[0]*v[1] + u[1]*v[0])
    def inv(s, u):
        n = u[0]*u[0] - s.c*u[1]*u[1]
        return (u[0]/n, -u[1]/n)
    def iszero(s, u): return u[0] == 0 and u[1] == 0

def mat2(x, a, c):
    # x = x0 + x1 i + x2 j + x3 k  -> 2x2 over Q(r)
    x0, x1, x2, x3 = x
    return [[(x0, x2), (x1*a, -x3*a)],
            [(x1, x3), (x0, -x2)]]

def mm(R, A, B):
    return [[R.add(R.mul(A[i][0], B[0][j]), R.mul(A[i][1], B[1][j])) for j in range(2)] for i in range(2)]
def madd(R, A, B): return [[R.add(A[i][j], B[i][j]) for j in range(2)] for i in range(2)]
def adj(A):  # conj = adjugate
    return [[A[1][1], (-A[0][1][0], -A[0][1][1])], [(-A[1][0][0], -A[1][0][1]), A[0][0]]]
def tr(R, A):
    t = R.add(A[0][0], A[1][1]); assert t[1] == 0; return t[0]

def coords(R, M, a, c):
    I = mat2([0,1,0,0], a, c); J = mat2([0,0,1,0], a, c); K = mat2([0,0,0,1], a, c)
    return [tr(R, M)/2, tr(R, mm(R, M, I))/(2*a), tr(R, mm(R, M, J))/(2*c), tr(R, mm(R, M, K))/(-2*a*c)]

def det_QR(R, M):
    n = len(M); M = [row[:] for row in M]; d = (F(1), F(0))
    for col in range(n):
        piv = next((r for r in range(col, n) if not R.iszero(M[r][col])), None)
        if piv is None: return (F(0), F(0))
        if piv != col:
            M[col], M[piv] = M[piv], M[col]; d = (-d[0], -d[1])
        d = R.mul(d, M[col][col]); iv = R.inv(M[col][col])
        for r in range(col+1, n):
            f = R.mul(M[r][col], iv)
            if R.iszero(f): continue
            M[r] = [R.sub(M[r][k], R.mul(f, M[col][k])) for k in range(n)]
    return d

# ---- Hilbert symbols (own implementation)
def sqfree_int(x):
    x = F(x); return x.numerator * x.denominator   # same square class
def vp(n, p):
    v = 0
    while n % p == 0: n //= p; v += 1
    return v, n
def leg(u, p):
    t = pow(u % p, (p-1)//2, p); return 1 if t == 1 else -1
def hilbert(x, y, p):
    x = sqfree_int(x); y = sqfree_int(y)
    if p == 0:  # infinity
        return -1 if (x < 0 and y < 0) else 1
    al, u = vp(x, p); be, v = vp(y, p)
    if p != 2:
        s = (-1) ** (al*be*((p-1)//2) % 2)
        return s * (leg(u, p) ** (be % 2)) * (leg(v, p) ** (al % 2))
    eps = lambda w: ((w-1)//2) % 2
    om = lambda w: ((w*w-1)//8) % 2
    e = (eps(u)*eps(v) + al*om(v) + be*om(u)) % 2
    return -1 if e else 1

def prime_factors(n):
    n = abs(n); S = set(); p = 2
    while p*p <= n:
        while n % p == 0: S.add(p); n //= p
        p += 1
    if n > 1: S.add(n)
    return S

def diagonalize(S):
    n = len(S); S = [row[:] for row in S]; d = []
    for k in range(n):
        if S[k][k] == 0:
            j = next((j for j in range(k+1, n) if S[k][j] != 0), None)
            if j is None:
                if all(S[k][t] == 0 for t in range(n)): raise ValueError('degenerate')
            # add row/col j to k (x_k -> x_k + x_j) -> S[k][k] += 2S[k][j] + S[j][j]
            if S[j][j] + 2*S[k][j] == 0:  # try difference
                sgn = -1
            else:
                sgn = 1
            for t in range(n): S[k][t] += sgn*S[j][t]
            for t in range(n): S[t][k] += sgn*S[t][j]
            assert S[k][k] != 0
        p = S[k][k]; d.append(p)
        for r in range(k+1, n):
            f = S[r][k] / p
            if f == 0: continue
            for t in range(k, n): S[r][t] -= f*S[k][t]
        for r in range(k+1, n): S[k][r] = F(0)
        for r in range(k+1, n): S[r][k] = F(0)
    return d

def is_square_rational(x):
    x = F(x)
    if x <= 0: return False
    import math
    n, m = x.numerator, x.denominator
    return math.isqrt(n)**2 == n and math.isqrt(m)**2 == m

def hasse(d, p):
    h = 1
    for i in range(len(d)):
        for j in range(i+1, len(d)):
            h *= hilbert(d[i], d[j], p)
    return h

HYP = [F(1), F(-1)] * 6

def run():
    data = json.load(open(os.path.join(HERE, 'data.json')))
    out = []
    for A in data['algebras']:
        a, c, disc = A['a'], A['c'], A['disc']
        R = QR(F(c))
        Ib = [[F(1),0,0,0], [0,F(1),0,0], [0,0,F(1),0], [0,0,0,F(1)]]
        Ib = [[F(t) for t in v] for v in Ib]
        mats = [mat2(v, a, c) for v in Ib]
        for ti, Trec in enumerate(A['Ts']):
            T = [[mat2([F(t) for t in e], a, c) for e in row] for row in Trec['T']]
            # Nrd(T): 6x6 over Q(r)
            big = [[None]*6 for _ in range(6)]
            for r_ in range(3):
                for s_ in range(3):
                    for u in range(2):
                        for v in range(2):
                            big[2*r_+u][2*s_+v] = T[r_][s_][u][v]
            dN = det_QR(R, big); assert dN[1] == 0
            N = dN[0]
            if N == 0: continue
            Delta = -N
            # E on Q-basis (m,beta)
            idx = [(m, be) for m in range(3) for be in range(4)]
            E = [[tr(R, mm(R, mm(R, mats[b1], T[m1][m2]), adj(mats[b2]))) for (m2, b2) in idx] for (m1, b1) in idx]
            for x in range(12):
                for y in range(12): assert E[x][y] == -E[y][x]
            for qv in A['qs']:
                q = [F(t) for t in qv]; Mq = mat2(q, a, c)
                # left mult by q: column for basis (m,beta) = coords of q*beta in slot m
                Qm = [[F(0)]*12 for _ in range(12)]
                for col, (m, be) in enumerate(idx):
                    co = coords(R, mm(R, Mq, mats[be]), a, c)
                    for t in range(4): Qm[4*m+t][col] = co[t]
                b = -(mm(R, Mq, Mq)[0][0][0])   # q^2 = -b
                S = [[-sum(E[x][z]*Qm[z][y] for z in range(12)) for y in range(12)] for x in range(12)]
                for x in range(12):
                    for y in range(12): assert S[x][y] == S[y][x]
                d = diagonalize(S)
                pos = sum(1 for t in d if t > 0)
                prod = F(1)
                for t in d: prod *= t
                # S has integer entries; at odd p not dividing det(S) the lattice Z^12 is unimodular, so the
                # Hasse invariant there is trivial for S and for <1,-1>^6. det(S) = prod(d) (unit-triangular
                # elimination with the pivot swaps used being determinant-preserving up to sign^2).
                assert all(x.denominator == 1 for row in S for x in row)
                hyp = (pos == 6) and is_square_rational(prod)
                if hyp:
                    import math
                    r = math.isqrt(prod.numerator)       # prod is a positive integer square here
                    P = {2} | prime_factors(r)
                    hyp = all(hasse(d, p) == hasse(HYP, p) for p in P)
                Pp = {2} | prime_factors(int(b.numerator*b.denominator)) | prime_factors(Delta.numerator) \
                     | prime_factors(Delta.denominator) | prime_factors(2*a*c)
                pred = all(hilbert(-b, Delta, p) == hilbert(a, c, p) for p in Pp)
                out.append({'disc': disc, 'ti': ti, 'q': [str(t) for t in q], 'b': str(b), 'Delta': str(Delta),
                            'split_hyperbolic': hyp, 'pred_formula': pred, 'pos': pos})
        print('done', disc, file=sys.stderr)
    json.dump(out, open(os.path.join(HERE, 'path2_out.json'), 'w'))
    n = len(out); ag = sum(r['split_hyperbolic'] == r['pred_formula'] for r in out)
    print(f'path2: {n} cases; split {sum(r["split_hyperbolic"] for r in out)}; formula agrees {ag}/{n}; sig(6,6) in {sum(r["pos"]==6 for r in out)}')

if __name__ == '__main__':
    run()
