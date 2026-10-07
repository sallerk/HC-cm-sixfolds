# verify_indep.py -- INDEPENDENT VERIFIER (pure Python + sympy; no PARI, no Sage/GAP, no code from the enumeration folder).
# Rebuilds the Galois actions from an explicit algebraic model (not from PARI's nfgaloisconj):
#   L = Q(a, b, i),  a^2 = 3+sqrt2, b^2 = 3-sqrt2  (so a^2+b^2 = 6, a^2 b^2 = 7),  Gal(L/Q) = <s,t> x <c> = D4 x C2:
#     s: a->b, b->-a, i->i;   t: a->a, b->-b, i->i;   c: i->-i.
#   E8 = Q(a, i);  E4* = Q(w), w = i(a-b)/sqrt2  (w^2 = sqrt7 - 3, min poly x^4+6x^2+2);  K1 = Q(i);  K2 = Q(i*sqrt2).
#   X1: S5 x C2 on (Hom(E0,R) x {+-i}) u {+-i}.
# Embeddings are evaluated NUMERICALLY (g applied to the generator, matched to the roots) -- a different route from PARI.
# For every CM type it recomputes: d = g - rank(mu-orbit), the saturated defect lattice, and the minimal dimension D at
# which imaginary-quadratic Weil characters WITH MULTIPLICITIES generate Lambda_U as a lattice (own code, Dmax=16),
# plus their Q-rank and index.  Output: verify_out.txt
import itertools, cmath, os
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

HERE = os.path.dirname(os.path.abspath(__file__))

def rank(rows):
    return Matrix(rows).rank() if rows else 0

def sat_kernel(mu, g):
    ns = Matrix(mu).nullspace()
    if not ns:
        return []
    # integral basis of the saturated kernel: clear denominators, then saturate via Smith form trick
    from sympy import lcm, gcd
    B = []
    for v in ns:
        den = lcm([x.q for x in v])
        B.append([int(x * den) for x in v])
    # saturate: lattice L = Q-span(B) cap Z^g.  Use: L = { x in Z^g : mu x = 0 } directly by HNF of kernel of mu over Z.
    return B

def snf_index(rows, d):
    M = Matrix(rows)
    S = smith_normal_form(M, domain=ZZ)
    diag = [abs(S[i, i]) for i in range(min(S.shape)) if S[i, i] != 0]
    assert len(diag) == d
    p = 1
    for x in diag:
        p *= int(x)
    return p

def analyse(points, conj, gens, factors, iq_fields, Dmax=16, label='', subf=None):
    """points: list of hashable embedding ids (2g of them); conj: dict id->id; gens: list of dicts id->id (group gens);
    factors: list of lists of point ids (each a G-orbit = one simple factor's CM field embeddings);
    iq_fields: list of (name, {factor_index: set of point ids restricting to tau0}) ."""
    n = len(points); g = n // 2
    # representatives / U-coordinates
    reps = []
    seen = set()
    for p in points:
        if p not in seen:
            reps.append(p); seen.add(p); seen.add(conj[p])
    coord = {}
    for k, p in enumerate(reps):
        coord[p] = (k, 1); coord[conj[p]] = (k, -1)
    def uvec(S):
        v = [0] * g
        for p in S:
            k, s = coord[p]; v[k] += s
        return v
    # group closure
    idp = {p: p for p in points}
    G = {tuple(idp[p] for p in points)}
    frontier = [idp]
    elems = [idp]
    while frontier:
        nf = []
        for h in frontier:
            for s in gens:
                c = {p: s[h[p]] for p in points}
                key = tuple(c[p] for p in points)
                if key not in G:
                    G.add(key); elems.append(c); nf.append(c)
        frontier = nf
    # CM types
    pairs = [(p, conj[p]) for p in reps]
    out = []
    seen_types = set()
    for choice in itertools.product([0, 1], repeat=g):
        Phi = frozenset(pr[c] for pr, c in zip(pairs, choice))
        if Phi in seen_types:
            continue
        orb = set(frozenset(h[p] for p in Phi) for h in elems)
        seen_types |= orb
        mu = [uvec(T) for T in orb]          # +1/-1 vectors
        r = rank(mu)
        d = g - r
        # primitivity of each factor type: not induced from a proper CM subfield <=> no nontrivial block system
        # compatible with Phi -- checked below only against the CM subfields listed in subf (analyse_sage.py checks all block systems).
        # IQ Weil characters with multiplicities
        vecs = []
        for name, fib in iq_fields:
            fidx = sorted(fib)
            b = [uvec(fib[i]) for i in fidx]
            t = [2 * len(Phi & fib[i]) - len(fib[i]) for i in fidx]
            sz = [len(fib[i]) for i in fidx]
            for cv in itertools.product(*[range(-(Dmax // s), Dmax // s + 1) for s in sz]):
                D = sum(abs(x) * s for x, s in zip(cv, sz))
                if D == 0 or D > Dmax or sum(x * y for x, y in zip(cv, t)):
                    continue
                w = [sum(x * bb[j] for x, bb in zip(cv, b)) for j in range(g)]
                if any(w):
                    assert all(sum(a * m for a, m in zip(w, mm)) == 0 for mm in mu)
                    vecs.append((D, w))
        qrank = rank([w for D, w in vecs]) if vecs else 0
        Dmin, idx_at = (0, 1) if d == 0 else (None, None)
        best_index = None
        if d > 0:
            for D in sorted(set(v[0] for v in vecs)):
                W = [w for DD, w in vecs if DD <= D]
                if rank(W) == d:
                    ix = snf_index(W, d)
                    best_index = ix
                    if ix == 1:
                        Dmin = D; break
        prim = True
        for fi, fl in (subf or {}).items():
            for kf in fl:
                sizes = {}
                for p in factors[fi]:
                    sizes[kf(p)] = sizes.get(kf(p), 0) + 1
                assert len(set(sizes.values())) == 1 and 1 < len(sizes) < len(factors[fi]), sizes
        for fi, fl in (subf or {}).items():
            Pi = Phi & set(factors[fi])
            for kf in fl:
                fibs = {}
                for p in factors[fi]:
                    fibs.setdefault(kf(p), set()).add(p)
                if all(F <= Pi or not (F & Pi) for F in fibs.values()):
                    prim = False
        sig = {name: {i: 2 * len(Phi & fib[i]) - len(fib[i]) for i in fib} for name, fib in iq_fields}
        out.append(dict(Phi=sorted(map(str, Phi)), norb=len(orb), d=d, prim=prim, IQDmin=Dmin, IQqrank=qrank,
                        IQindex_at_full_rank=best_index, sig=sig))
    return len(elems), out

lines = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); lines.append(s)

# ---------------- D4 x C2 model ----------------
sq2 = 2 ** 0.5
A0 = (3 + sq2) ** 0.5; B0 = (3 - sq2) ** 0.5
# a Galois element = (image of a, image of b, image of i) as formal symbols in {+-a,+-b} x {+-i}; numeric evaluation:
VAL = {'a': A0, '-a': -A0, 'b': B0, '-b': -B0, 'i': 1j, '-i': -1j}
def neg(x): return x[1:] if x.startswith('-') else '-' + x
def compose(g1, g2):   # (g1 o g2): first g2 then g1;  g = dict on symbols a,b,i
    def ap(gg, x):
        return neg(gg[x[1:]]) if x.startswith('-') else gg[x]
    return {k: ap(g1, g2[k]) for k in 'abi'}
s = {'a': 'b', 'b': '-a', 'i': 'i'}
t = {'a': 'a', 'b': '-b', 'i': 'i'}
c = {'a': 'a', 'b': 'b', 'i': '-i'}
# closure
grp = []
fr = [{'a': 'a', 'b': 'b', 'i': 'i'}]
keys = set()
while fr:
    nf = []
    for h in fr:
        k = tuple(h[x] for x in 'abi')
        if k in keys: continue
        keys.add(k); grp.append(h)
        for gg in (s, t, c):
            nf.append(compose(gg, h))
    fr = nf
assert len(grp) == 16, len(grp)
def num(h):    # numeric images of a, b, i
    return VAL[h['a']], VAL[h['b']], VAL[h['i']]
# check relations a^2+b^2=6, a^2 b^2=7 preserved (automatic) and that sqrt2 = a^2-3 maps consistently
def elem_vals(h):
    a, b, i = num(h)
    r2 = a * a - 3
    return dict(a=a, b=b, i=i, r2=r2, w=i * (a - b) / r2, sqrtm2=i * r2, r7=a * b)
idh = grp[[tuple(h[x] for x in 'abi') for h in grp].index(('a', 'b', 'i'))]
base = elem_vals(idh)
# embeddings of E8 = Q(a,i): identified by (g(a), g(i)); of Q(w): by numeric value g(w); of Q(i): g(i); of Q(sqrt-2): g(i sqrt2)
def emb_E8(h): return ('E8', h['a'], h['i'])
def emb_w(h): v = elem_vals(h)['w']; return ('W', round(v.real, 6) + 0.0, round(v.imag, 6) + 0.0)
def emb_K1(h): return ('K1', h['i'])
def emb_K2(h): v = elem_vals(h)['sqrtm2']; return ('K2', round(v.imag, 6) + 0.0)

def build_D4(kinds):
    funcs = dict(E8=emb_E8, W=emb_w, K1=emb_K1, K2=emb_K2)
    pts, fac = [], []
    for kd in kinds:
        f = funcs[kd]
        S = []
        for h in grp:
            e = f(h)
            if e not in S: S.append(e)
        fac.append(S); pts += S
    # Galois action: g acts on embedding f(h) as f(g o h)
    def act(gg):
        m = {}
        for kd, S in zip(kinds, fac):
            f = funcs[kd]
            for h in grp:
                m[f(h)] = f(compose(gg, h))
        return m
    gens = [act(s), act(t), act(c)]
    conj = act(c)
    return pts, conj, gens, fac

def iq_fibres_D4(kinds, fac):
    # K1 = Q(i): restriction of an embedding to K1 determined by image of i; K2 = Q(i sqrt2): by image of i*sqrt2
    res = {}
    for name, test in (('Q(i)', lambda h: h['i'] == 'i'), ('Q(sqrt-2)', lambda h: elem_vals(h)['sqrtm2'].imag > 0)):
        fib = {}
        for fi, (kd, S) in enumerate(zip(kinds, fac)):
            if kd == 'W':
                continue          # Q(w) has no imaginary quadratic subfield (D4-primitive quartic)
            if kd == 'K1' and name != 'Q(i)': continue
            if kd == 'K2' and name != 'Q(sqrt-2)': continue
            f = dict(E8=emb_E8, K1=emb_K1, K2=emb_K2)[kd]
            fib[fi] = set(f(h) for h in grp if test(h))
        res[name] = fib
    return [(k, v) for k, v in res.items() if v]

def isb(x): return x in ('b', '-b')
E8SUB = [lambda p: p[2], lambda p: (p[2] == 'i') != isb(p[1]), lambda p: (p[2], isb(p[1])),
         lambda p: p[1] if p[2] == 'i' else neg(p[1])]   # Q(i), Q(sqrt-2), Q(zeta8), Q(i a)
for kinds, nm in ((['E8', 'W'], '(4,2): A(E8=Q(a,i)) x S(Q(w)), w^2 = sqrt7-3'),
                  (['E8', 'K1', 'K2'], '(4,1,1): A(E8) x E_{Q(i)} x E_{Q(sqrt-2)}')):
    pts, conj, gens, fac = build_D4(kinds)
    assert len(pts) == 12 and all(conj[conj[p]] == p and conj[p] != p for p in pts)
    iq = iq_fibres_D4(kinds, fac)
    order, res = analyse(pts, conj, gens, fac, iq, 16, nm, subf={0: E8SUB})
    P('==', nm, '| |G| image =', order, '| IQ fields:', [k for k, v in iq])
    prof = {}
    for r in res:
        key = (r['prim'], r['d'], r['IQDmin'], r['IQqrank'], r['IQindex_at_full_rank'], r['norb'])
        prof[key] = prof.get(key, 0) + 1
        if r['d'] >= 1 and r['prim']:
            P('   d=%d norb=%d IQDmin16=%s IQqrank=%d IQindex=%s sig=%s' % (r['d'], r['norb'], r['IQDmin'], r['IQqrank'],
                                                                         r['IQindex_at_full_rank'], r['sig']))
    P('   profile (d, IQDmin, IQqrank, IQindex, orbit size) -> #type-orbits:', sorted(prof.items(), key=str))

# ---------------- X1: S5 x C2 ----------------
pts = [('E', k, e) for k in range(5) for e in (1, -1)] + [('K', e) for e in (1, -1)]
def mk(sig, flip):
    m = {}
    for p in pts:
        if p[0] == 'E':
            m[p] = ('E', sig[p[1]], p[2] * flip)
        else:
            m[p] = ('K', p[1] * flip)
    return m
gens = [mk([1, 2, 3, 4, 0], 1), mk([1, 0, 2, 3, 4], 1), mk([0, 1, 2, 3, 4], -1)]
conj = gens[2]
fac = [[p for p in pts if p[0] == 'E'], [p for p in pts if p[0] == 'K']]
iq = [('Q(i)', {0: set(p for p in fac[0] if p[2] == 1), 1: {('K', 1)}})]
order, res = analyse(pts, conj, gens, fac, iq, 16, 'X1', subf={0: [lambda p: p[2]]})
P('== X1 (5,1): A(E0(i), E0 S5 quintic) x E_{Q(i)} | |G| image =', order)
prof = {}
for r in res:
    key = (r['prim'], r['d'], r['IQDmin'], r['IQqrank'], r['IQindex_at_full_rank'], r['norb'])
    prof[key] = prof.get(key, 0) + 1
    if r['d'] >= 1 and r['prim']:
        P('   d=%d norb=%d IQDmin16=%s IQqrank=%d IQindex=%s sig=%s' % (r['d'], r['norb'], r['IQDmin'], r['IQqrank'],
                                                                     r['IQindex_at_full_rank'], r['sig']))
P('   profile:', sorted(prof.items(), key=str))
open(os.path.join(HERE, 'verify_out.txt'), 'w').write('\n'.join(lines) + '\n')
