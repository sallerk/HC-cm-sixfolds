"""Explicit ordinary, geometrically simple abelian sixfolds over F_p (p = 43, 127) with NON-maximal
angle rank, i.e. with exceptional Tate classes, realised as Fermat-21 Jacobi-sum factors.
Main computation (PARI via cypari):
- pi = J(chi^a, chi^b) in Z[zeta21], with chi of order 21 on F_p^x;
- Weil polynomial = charpoly; checks pi*conj(pi) = p;
- ordinary: p does not divide the T^6 coefficient;
- geometrically simple: charpoly(pi^42) irreducible (all roots of unity in Q(zeta21) are 42nd roots);
- angle rank = rank(valuation matrix over sigma_s and primes above p) - 1.
Independent check (pure Python, no field arithmetic): Stickelberger exponent set
S = {t : <ta> + <tb> >= m}; rank of the G-circulant matrix [1_S(s t)] must equal the PARI rank."""
from math import gcd
from itertools import product
from fractions import Fraction
from cypari import pari

m = 21
units = [t for t in range(1, m) if gcd(t, m) == 1]
cyc = pari(f'polcyclo({m}, x)')
nf = pari.nfinit(cyc)

def jacobi(p, a, b):
    g = int(pari.znprimroot(p).lift())
    dlog, y = {}, 1
    for k in range(p - 1):
        dlog[y] = k; y = y * g % p
    v = [0] * m
    for t in range(2, p):
        v[(a * dlog[t] + b * dlog[(1 - t) % p]) % m] += 1
    return v

def elt(v, s=1):
    return pari('Mod(' + '+'.join(f'({c})*x^{(k*s) % m}' for k, c in enumerate(v) if c) + f', {cyc})')

def qrank(M):  # exact rank over Q, pure Python
    M = [[Fraction(x) for x in row] for row in M]; r = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]; M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
    return r

rows = []
for p in (43, 127):
    assert p % m == 1
    dec = pari.idealprimedec(nf, p)
    assert len(dec) == 12
    seen = set()
    for a, b in product(range(1, m), repeat=2):
        c = (-a - b) % m
        if c == 0 or any(gcd(z, m) != 1 for z in (a, b, c)):
            continue
        key = min(tuple(sorted(((a*t) % m, (b*t) % m, (c*t) % m))) for t in units)
        if key in seen:
            continue
        seen.add(key)
        v = jacobi(p, a, b)
        pi = elt(v)
        pibar = elt(v, m - 1)
        normok = (pi * pibar == p)
        P = pari.charpoly(pi)
        co = [int(P.polcoef(k)) for k in range(12, -1, -1)]
        ordinary = co[6] % p != 0
        irr = bool(pari.polisirreducible(P))
        geo = bool(pari.polisirreducible(pari.charpoly(pi**42)))
        V = [[int(pari('(n,a,q)->nfeltval(n,a,q)')(nf, elt(v, s).lift(), P_)) for P_ in dec] for s in units]
        rk = int(pari.matrank(pari.matrix(len(V), 12, [x for r_ in V for x in r_])))
        S = {t for t in units if (t*a % m) + (t*b % m) >= m}
        rk_st = qrank([[1 if (s*t) % m in S else 0 for t in units] for s in units])
        row = dict(p=p, abc=key, norm_ok=normok, ordinary=ordinary, irreducible=irr,
                   geom_simple=geo, rank_val=rk, rank_stickelberger=rk_st, angle_rank=rk - 1,
                   weil_poly=co)
        rows.append(row)
        print({k: row[k] for k in row if k != 'weil_poly'}, flush=True)
        if rk - 1 < 6:
            print('   Weil poly coeffs (T^12..T^0):', co, flush=True)
import json, os
json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'jacobi21.json'), 'w'), indent=1)
bad = [r for r in rows if not (r['norm_ok'] and r['irreducible'] and r['rank_val'] == r['rank_stickelberger'])]
print('rows', len(rows), 'consistency failures', len(bad))
