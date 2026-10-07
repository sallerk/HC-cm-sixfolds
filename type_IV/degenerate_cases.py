# degenerate_cases.py -- tests of the split-lemma PROCEDURES in degenerate situations.
#  (1) Product case F = Q x Q (L = K x K, A ~ A1 x A2 with K-signatures (2,1), (1,2)):
#      run proof 2 per factor (isotropic vectors by qfsolve on each 6-dim trace form), a = (1/d1, -1/d2),
#      check W = K e1 + K e2 + K (g1,g2) is h_a-isotropic, Witt index 6, Hilbert criterion.
#  (2) Forced primes by local type (inert/ramified in K and F, p = 2) over the (d,m) grid, with the
#      symbol values of the two ingredients (N Delta,-d)_p and (-d_F,-d)_p.
#  (3) Obstruction control: replacing c by c*q for a forced-type prime q makes the quaternary form
#      anisotropic at q (qfsolve returns q): the lemma depends on the shape c = -N(Delta) d_F.
import os, random, json
from fractions import Fraction as Fr
import numpy as np
from cypari import pari
from biq import Ctx, witt_index
import split_checks as S

HERE = os.path.dirname(os.path.abspath(__file__))   # JSON output is written next to this script
rng = random.Random(77)

# ---------- K-Hermitian helpers (K = Q(s), elements (x, y) = x + y s) ----------
def kmul(d, a, b): return (a[0]*b[0] - d*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def kconj(a): return (a[0], -a[1])
def kadd(a, b): return (a[0] + b[0], a[1] + b[1])
def kform(d, H, x, y):
    z = (Fr(0), Fr(0))
    for i in range(len(x)):
        for j in range(len(y)):
            z = kadd(z, kmul(d, kmul(d, x[i], H[i][j]), kconj(y[j])))
    return z
def ksig(d, H):
    A = np.array([[complex(float(H[i][j][0]), float(H[i][j][1]) * d ** 0.5) for j in range(3)] for i in range(3)])
    ev = np.linalg.eigvalsh(A)
    return (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()))
def rand_kherm(d):
    H = [[None]*3 for _ in range(3)]
    for i in range(3):
        H[i][i] = (Fr(rng.randint(-6, 6)), Fr(0))
        for j in range(i + 1, 3):
            z = (Fr(rng.randint(-3, 3)), Fr(rng.randint(-3, 3)))
            H[i][j] = z; H[j][i] = kconj(z)
    return H
def kdet(d, H):
    rows = [','.join('Mod(%s+(%s)*t,t^2+%d)' % (str(x), str(y), d) for x, y in r) for r in H]
    D = pari.lift(pari('matdet([' + ';'.join(rows) + '])'))
    assert int(pari.polcoef(D, 1, 't') == 0)
    c0 = pari.polcoef(D, 0, 't')
    return Fr(int(pari.numerator(c0)), int(pari.denominator(c0)))
def trace_gram(d, H):
    """rational Gram matrix of x -> H(x,x) on the Q-basis {e_i, s e_i}"""
    n = len(H)
    Q = [[Fr(0)]*(2*n) for _ in range(2*n)]
    for p in range(n):
        for q in range(n):
            re, im = H[p][q]
            for e in range(2):
                for f in range(2):
                    x, y = re, im
                    for _ in range(e): x, y = -d*y, x
                    for _ in range(f): x, y = d*y, -x
                    Q[2*p + e][2*q + f] = x
    return Q
def k_isotropic(d, H):
    Q = trace_gram(d, H)
    den = 1
    for r in Q:
        for x in r:
            den = den * x.denominator
    G = pari.matrix(6, 6, [int(x*den) for r in Q for x in r])
    v = pari.qfsolve(G)
    if str(pari.type(v)) == 't_MAT':
        v = v[0]
    assert str(pari.type(v)) == 't_COL', v
    v = [Fr(int(pari.numerator(t)), int(pari.denominator(t))) for t in v]
    x = [(v[2*i], v[2*i + 1]) for i in range(3)]
    assert kform(d, H, x, x) == (0, 0)
    return x
def kinv(d, a):
    n = a[0]*a[0] + d*a[1]*a[1]
    return (a[0]/n, -a[1]/n)
def cross3(d, a, b):
    def sub(x, y): return (x[0] - y[0], x[1] - y[1])
    return [sub(kmul(d, a[1], b[2]), kmul(d, a[2], b[1])), sub(kmul(d, a[2], b[0]), kmul(d, a[0], b[2])),
            sub(kmul(d, a[0], b[1]), kmul(d, a[1], b[0]))]
def factor_construction(d, H):
    e = k_isotropic(d, H)
    j = next(j for j in range(3) if kform(d, H, e, [(Fr(1), Fr(0)) if k == j else (Fr(0), Fr(0)) for k in range(3)]) != (0, 0))
    ej = [(Fr(1), Fr(0)) if k == j else (Fr(0), Fr(0)) for k in range(3)]
    lam = kconj(kinv(d, kform(d, H, e, ej)))
    f = [kmul(d, lam, z) for z in ej]
    hff = kform(d, H, f, f); half = (hff[0]/2, hff[1]/2)
    fp = [(fi[0] - kmul(d, half, ei)[0], fi[1] - kmul(d, half, ei)[1]) for fi, ei in zip(f, e)]
    assert kform(d, H, fp, fp) == (0, 0) and kform(d, H, e, fp) == (1, 0)
    def col(x):
        out = []
        for i in range(3):
            z = (Fr(0), Fr(0))
            for jj in range(3):
                z = kadd(z, kmul(d, H[i][jj], kconj(x[jj])))
            out.append(z)
        return out
    g = cross3(d, col(e), col(fp))
    assert kform(d, H, g, e) == (0, 0) and kform(d, H, g, fp) == (0, 0)
    delta = kform(d, H, g, g)
    assert delta[1] == 0 and delta[0] != 0
    return e, g, delta[0]

def product_case(d, n=15):
    out = dict(d=d, cases=0, ok=0)
    for _ in range(n):
        while True:
            H1, H2 = rand_kherm(d), rand_kherm(d)
            if ksig(d, H1) == (2, 1) and ksig(d, H2) == (1, 2):
                break
        out['cases'] += 1
        e1, g1, d1 = factor_construction(d, H1)
        e2, g2, d2 = factor_construction(d, H2)
        mixed = d1 > 0 > d2
        a1, a2 = 1/d1, -1/d2                        # a = mu/delta with mu = (1,-1), Tr(mu) = 0
        Z = [(Fr(0), Fr(0))]*3
        # h_a on K^6 = K^3 + K^3 (block diagonal a1 H1 + a2 H2)
        Hb = [[(Fr(0), Fr(0))]*6 for _ in range(6)]
        for i in range(3):
            for j in range(3):
                Hb[i][j] = (a1*H1[i][j][0], a1*H1[i][j][1])
                Hb[3+i][3+j] = (a2*H2[i][j][0], a2*H2[i][j][1])
        W = [e1 + Z, Z + e2, g1 + g2]
        iso = all(kform(d, Hb, x, y) == (0, 0) for x in W for y in W)
        wi, _ = witt_index(trace_gram(d, Hb))
        det = a1**3 * a2**3 * kdet(d, H1) * kdet(d, H2)
        C = Ctx(d, 2)
        hil = C.is_norm_K(-det)
        out['ok'] += int(mixed and a1 > 0 and a2 > 0 and iso and wi == 6 and hil)
    return out

def forced_table(per=3):
    pairs = [(d, m) for d in [1, 2, 3, 5, 6, 7, 11, 15, 19, 23, 31, 35]
             for m in [2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 21, 30]]
    tab = {}
    for d, m in pairs:
        C = Ctx(d, m)
        for _ in range(per):
            H = S.good_herm(C)
            if H is None:
                continue
            ND = C.NF(C.detL(H)); c = -ND * C.dF
            for p in C.forced_primes(c):
                ramK = int(pari('nfdisc(x^2+%d)' % d)) % p == 0
                ramF = C.dF % p == 0
                typ = ('p=2 ' if p == 2 else '') + ('ramified-both' if ramK and ramF else
                       ('inert-both' if not ramK and not ramF else 'MIXED?'))
                key = typ
                e1 = C.hilbert(ND, -d, p); e2 = C.hilbert(Fr(-C.dF), -d, p)
                rec = tab.setdefault(key, dict(count=0, sym_ND=set(), sym_dF=set(), examples=[]))
                rec['count'] += 1; rec['sym_ND'].add(e1); rec['sym_dF'].add(e2)
                if len(rec['examples']) < 6 and (d, m, p) not in rec['examples']:
                    rec['examples'].append((d, m, p))
    for k in tab:
        tab[k]['sym_ND'] = sorted(tab[k]['sym_ND']); tab[k]['sym_dF'] = sorted(tab[k]['sym_dF'])
    return tab

def obstruction_control(n=20):
    """c' = c q with q prime, -dm a square mod q, q inert in K (so q is a forced prime and (c q,-d)_q = -1)."""
    res = []
    cnt = 0
    for (d, m) in [(1, 2), (1, 7), (3, 5), (2, 3), (5, 3), (7, 21), (15, 10), (3, 15)]:
        C = Ctx(d, m)
        H = S.good_herm(C)
        ND = C.NF(C.detL(H)); c = -ND * C.dF
        q = 3
        while True:
            q = int(pari.nextprime(q + 1))
            if (2*d*m*c.numerator*c.denominator) % q == 0:
                continue
            if pari.kronecker(-d*m, q) == 1 and pari.kronecker(-d, q) == -1:
                break
        c2 = c * q
        den = c2.denominator
        v = pari.qfsolve(pari.matdiagonal([den*den, den*den*d, -c2.numerator*den, c2.numerator*den*m]))
        v0 = pari.qfsolve(pari.matdiagonal([c.denominator**2, c.denominator**2*d, -c.numerator*c.denominator,
                                             c.numerator*c.denominator*m]))
        res.append(dict(d=d, m=m, q=q, qfsolve_cq=str(v), qfsolve_c_is_solution=str(pari.type(v0)) in ('t_COL', 't_MAT')))
    return res

if __name__ == '__main__':
    prod = [product_case(d) for d in [1, 2, 3, 5, 7, 15, 23]]
    print('PRODUCT CASE (F = Q x Q):', prod)
    tab = forced_table()
    print('FORCED PRIMES BY TYPE:')
    for k, v in tab.items():
        print('  ', k, v)
    oc = obstruction_control()
    print('OBSTRUCTION CONTROL:')
    for r in oc:
        print('  ', r)
    json.dump(dict(product=prod, forced=tab, control=oc), open(os.path.join(HERE, 'degenerate_cases.json'), 'w'), indent=1, default=str)
