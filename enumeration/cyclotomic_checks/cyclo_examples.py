# cyclo_examples.py -- independent re-check (no code shared with the enumeration scripts): defect of the Hodge group for
# explicit cyclotomic CM abelian varieties, using G = (Z/M)^x acting on embeddings by multiplication.
# Factors:
#   ('cyc', d, Phi)   : CM field Q(zeta_d), embeddings <-> a in (Z/d)^x, CM type Phi (list of a's)
#   ('iq', p, coset)  : K = Q(sqrt(-p)) inside Q(zeta_p) (p = 3 mod 4), embeddings <-> {squares, non-squares}
#                       of (Z/p)^x, CM type = the given coset ('sq' or 'nsq')
# Output: g, dim MT = rank span{t.mu}, defect d = g + 1 - dim MT, a rational basis of the odd defect
# space (characters of U_E vanishing on Hg), written per factor.
import sys
from math import gcd
from sympy import Matrix, Rational, lcm

def units(n):
    return [a for a in range(1, n) if gcd(a, n) == 1]

def build(factors):
    M = 1
    for f in factors:
        M = lcm(M, f[1])
    pts = []  # (factor index, label)
    for i, f in enumerate(factors):
        if f[0] == 'cyc':
            pts += [(i, a) for a in units(f[1])]
        else:
            pts += [(i, 'sq'), (i, 'nsq')]
    idx = {p: k for k, p in enumerate(pts)}
    def act(t, p):
        i, lab = p
        f = factors[i]
        if f[0] == 'cyc':
            return (i, (t * lab) % f[1])
        q = f[1]
        sq = set((x * x) % q for x in range(1, q))
        tsq = (t % q) in sq
        return (i, lab if tsq else ('nsq' if lab == 'sq' else 'sq'))
    mu = [0] * len(pts)
    for i, f in enumerate(factors):
        for lab in (f[2] if f[0] == 'cyc' else [f[2]]):
            mu[idx[(i, lab)]] = 1
    G = units(M)
    rows = []
    for t in G:
        v = [0] * len(pts)
        for k, p in enumerate(pts):
            if mu[k]:
                v[idx[act(t, p)]] = 1
        rows.append(v)
    # sanity: complex conjugation t = -1 maps mu to its complement on each factor (CM type)
    conj = rows[G.index(M - 1)]
    assert all(conj[k] + mu[k] == 1 for k in range(len(pts))), 'not a CM type'
    S = Matrix(rows)
    rk = S.rank()
    g = len(pts) // 2
    # odd defect: c with c[rho x] = -c[x], sum_x c_x (t.mu)_x = 0 for all t
    rho = {k: idx[act(M - 1, p)] for k, p in enumerate(pts)}
    reps = [k for k in range(len(pts)) if k < rho[k]]
    # variables c_r for r in reps; c[rho r] = -c_r
    A = []
    for r_ in rows:
        A.append([r_[k] - r_[rho[k]] for k in reps])
    ns = Matrix(A).nullspace()
    return pts, reps, g, rk, ns

def show(name, factors, expect_d=None):
    pts, reps, g, rk, ns = build(factors)
    d = g + 1 - rk
    print(f'== {name}: g={g} dimMT={rk} defect d={d} (odd nullspace dim {len(ns)})'
          + ('' if expect_d is None else f'  expected d={expect_d} -> {"OK" if d == expect_d else "MISMATCH"}'))
    assert len(ns) == d, 'defect bookkeeping mismatch'
    for v in ns:
        v = v / max(abs(x) for x in v if x != 0)
        print('   defect vector:', {str(pts[k]): str(v[j]) for j, k in enumerate(reps) if v[j] != 0})
    return d

if __name__ == '__main__':
    res = {}
    # controls
    res['J7'] = show('J(y^2=x^7+1): Q(zeta7), Phi={1,2,3}', [('cyc', 7, [1, 2, 3])], 0)
    res['Q21'] = show('Q(zeta21), Phi={1,2,4,5,8,10} (Shioda/APFV 12T2)', [('cyc', 21, [1, 2, 4, 5, 8, 10])], 1)
    # Milne Ex 1.12, n=3: Q(zeta7) with type {1,5,3} (signature (1,2) w.r.t. Q(sqrt-7)), times E_{sqrt-7} type 'sq'
    res['Milne3'] = show('Milne Ex1.12 n=3: Q(zeta7){1,3,5} x E_{sqrt-7}(sq)', [('cyc', 7, [1, 3, 5]), ('iq', 7, 'sq')], 1)
    # flagship (5,1): J_11 x E_{sqrt-11}; type {1..5}, K-signature (4,1): squares mod 11 = {1,3,4,5,9}
    res['J11xE'] = show('J(y^2=x^11+1) x E_{sqrt-11}(nsq)', [('cyc', 11, [1, 2, 3, 4, 5]), ('iq', 11, 'nsq')], 1)
    res['J11xE_sq'] = show('J(y^2=x^11+1) x E_{sqrt-11}(sq) [other E-type]', [('cyc', 11, [1, 2, 3, 4, 5]), ('iq', 11, 'sq')], None)
    # flagship (4,2): factor A_15 x S_5 of J(y^2=x^15+1)
    res['A15xS5'] = show('A_15{1,2,4,7} x S_5{1,2} (J_15 factor)', [('cyc', 15, [1, 2, 4, 7]), ('cyc', 5, [1, 2])], 2)
    res['J15'] = show('J(y^2=x^15+1) = A_15 x S_5 x E_3', [('cyc', 15, [1, 2, 4, 7]), ('cyc', 5, [1, 2]), ('cyc', 3, [1])], None)
    print(res)
