# examples.py -- explicit test cases, each K checked by BOTH code paths:
#   path1.gram_split  (quaternion projection, det over K, Hilbert symbols, PARI)
#   path2 trace-form hyperbolicity (own arithmetic, own Hilbert symbols, no PARI)
# plus the criterion  split(K=Q(sqrt -b)) <=> (-b, Delta) ~ D,  Delta = -Nrd(T)  (note, Proposition 6.2).
# Cases include those relevant to Abdulali's list (Delta = -1, D containing sqrt-1 or sqrt-3).
# Prints the cases (saved as examples_out.txt) and writes examples_out.json next to this script.
import itertools, json, os
from fractions import Fraction as F
from path1 import nrd, nrd_matrix, hilb, primes_of
from exists import gram_split, is_sq_Qp, disc_primes
import path2 as P2

HERE = os.path.dirname(os.path.abspath(__file__))

def hyperbolic_path2(a, c, Tint, q):
    R = P2.QR(F(c))
    Ib = [[F(int(i == j)) for j in range(4)] for i in range(4)]
    mats = [P2.mat2(v, a, c) for v in Ib]
    T = [[P2.mat2([F(t) for t in e], a, c) for e in row] for row in Tint]
    idx = [(m, be) for m in range(3) for be in range(4)]
    E = [[P2.tr(R, P2.mm(R, P2.mm(R, mats[b1], T[m1][m2]), P2.adj(mats[b2]))) for (m2, b2) in idx] for (m1, b1) in idx]
    Mq = P2.mat2([F(t) for t in q], a, c)
    Qm = [[F(0)]*12 for _ in range(12)]
    for col, (m, be) in enumerate(idx):
        co = P2.coords(R, P2.mm(R, Mq, mats[be]), a, c)
        for t in range(4): Qm[4*m+t][col] = co[t]
    S = [[-sum(E[x][z]*Qm[z][y] for z in range(12)) for y in range(12)] for x in range(12)]
    d = P2.diagonalize(S)
    pos = sum(1 for t in d if t > 0); prod = F(1)
    for t in d: prod *= t
    Ps = {2}
    for t in d: Ps |= P2.prime_factors(t.numerator) | P2.prime_factors(t.denominator)
    return (pos == 6) and P2.is_square_rational(prod) and all(P2.hasse(d, p) == P2.hasse(P2.HYP, p) for p in Ps)

def diagT(l1, l2, l3):
    z = [0, 0, 0, 0]
    return [[l1, z, z], [z, l2, z], [z, z, l3]]

cases = [
    # (label, a, c, T, list of q)
    ('E1 D=(-2,-5) disc5, T=diag(i,j,k), Delta=-100~-1; contains sqrt-3, not sqrt-1', -2, -5,
     diagT([0,1,0,0], [0,0,1,0], [0,0,0,1]), [[0,1,0,1], [0,1,0,0], [0,0,1,0], [0,1,1,0], [0,1,1,1], [0,0,1,1]]),
    ('E2 D=(-1,-7) disc7, T=diag(i,i,i), Delta=-1; contains sqrt-1, not sqrt-3', -1, -7,
     diagT([0,1,0,0], [0,1,0,0], [0,1,0,0]), [[0,1,0,0], [0,0,1,0], [0,1,1,0], [0,0,0,1], [0,1,0,1]]),
    ('E3 D=(-1,-1) disc2, T=diag(i,i,i), Delta=-1', -1, -1,
     diagT([0,1,0,0], [0,1,0,0], [0,1,0,0]), [[0,1,0,0], [0,1,1,0], [0,1,1,1], [0,1,2,0], [0,1,1,2]]),
    ('E4 D=(-1,-3) disc3, T=diag(i,i,i), Delta=-1', -1, -3,
     diagT([0,1,0,0], [0,1,0,0], [0,1,0,0]), [[0,1,0,0], [0,0,1,0], [0,1,1,0], [0,0,0,1], [0,1,0,1]]),
    ('E5 D=(-1,-3) disc3, T=diag(i,i,j) (Abdulali Ex.3.1.3), Delta=-3', -1, -3,
     diagT([0,1,0,0], [0,1,0,0], [0,0,1,0]), [[0,1,0,0], [0,0,1,0], [0,1,1,0], [0,0,0,1], [0,1,0,1]]),
    ('E6 D=(-1,-1) disc2, T=diag(i,i,i+j+k), Delta=-3 (Q(sqrt-3) in D)', -1, -1,
     diagT([0,1,0,0], [0,1,0,0], [0,1,1,1]), [[0,1,0,0], [0,1,1,0], [0,1,1,1], [0,1,2,0], [0,1,1,2], [0,1,2,3]]),
    ('E7 D=(-1,-1) disc2, Delta=-7 (Q(sqrt-7) not in D): T=diag(i, i+j, 3i+2j+k)', -1, -1,
     diagT([0,1,0,0], [0,1,1,0], [0,3,2,1]), [[0,1,0,0], [0,1,1,0], [0,1,1,1], [0,1,2,0], [0,1,1,2], [0,1,2,3], [0,2,3,5]]),
]
out = []
for label, a, c, T, qs in cases:
    TF = [[[F(t) for t in e] for e in row] for row in T]
    N = nrd_matrix(TF, a, c); Delta = -N
    dp = disc_primes(a, c)
    emb = all(not is_sq_Qp(Delta, p) for p in dp)
    print(f'{label}\n   disc(D) primes {dp}, Nrd(T) = {N}, Delta = {Delta}, Q(sqrt Delta) in D: {emb}')
    for q in qs:
        qF = [F(t) for t in q]
        b, detH, sp1 = gram_split(TF, qF, a, c)
        sp2 = hyperbolic_path2(a, c, T, q)
        P = primes_of(b, Delta, a, c)
        pred = all(hilb(-b, Delta, p) == hilb(a, c, p) for p in P)
        print(f'   q={q}  b={b}  det_K H={detH}  split(path1)={sp1}  split(path2)={sp2}  criterion={pred}')
        out.append({'case': label, 'q': q, 'b': str(b), 'detH': str(detH), 'split1': sp1, 'split2': sp2, 'pred': pred,
                    'Delta': str(Delta), 'emb': emb})
        assert sp1 == sp2 == pred
json.dump(out, open(os.path.join(HERE, 'examples_out.json'), 'w'), indent=0)
print('all examples: path1 == path2 == criterion')
