"""Nm_{E/K}(pi) = zeta * q^3 for the 15 rows of Table 3 of the note (Sec. 8).
Usage: python norm_zeta.py final_table.json > norm_zeta.out
Numerics (mpmath) only locate the factor P1 of the Weil polynomial P over K = Q(sqrt(-D)); the factorization
P = P1 * conj(P1) over K and the value zeta = P1(0)/q^3 (and its order) are then checked in exact arithmetic."""
import json, sys, itertools
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 120
rows = [r for r in json.load(open(sys.argv[1])) if r['name'].endswith('_deg')]
# K = Q(sqrt(-D)); element = (x, y) = x + y*sqrt(-D), x, y Fractions
def kmul(a, b, D): return (a[0]*b[0] - D*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def kconj(a): return (a[0], -a[1])
def polymul(P, Q, D):
    R = [(Fr(0), Fr(0))] * (len(P) + len(Q) - 1)
    for i, a in enumerate(P):
        for j, b in enumerate(Q):
            c = kmul(a, b, D); R[i+j] = (R[i+j][0] + c[0], R[i+j][1] + c[1])
    return R
for r in rows:
    D = r['K_d']; q = r['q']; P = r['weil_poly']  # coefficients T^12 .. T^0
    roots = mp.polyroots(P, maxsteps=400, extraprec=600)
    up = sorted([z for z in roots if mp.im(z) > 0], key=lambda z: mp.arg(z))
    assert len(up) == 6
    found = []
    for signs in itertools.product([0, 1], repeat=6):
        S = [z if s == 0 else mp.conj(z) for z, s in zip(up, signs)]
        c = [mp.mpc(1)]
        for z in S:  # multiply by (x - z)
            c = [c[0]] + [c[i] - z*c[i-1] for i in range(1, len(c))] + [-z*c[-1]]
        # round each coefficient to (x + y sqrt(-D)) with x, y in (1/2)Z
        P1 = []; ok = True
        for a in c:
            x = mp.re(a); y = mp.im(a) / mp.sqrt(D)
            xr, yr = mp.nint(2*x) / 2, mp.nint(2*y) / 2
            if abs(x - xr) > mp.mpf(10)**-40 or abs(y - yr) > mp.mpf(10)**-40: ok = False; break
            P1.append((Fr(int(2*xr), 2), Fr(int(2*yr), 2)))
        if ok: found.append(P1)
    assert len(found) == 2, (r['name'], len(found))
    P1 = found[0]
    prod = polymul(P1, [kconj(a) for a in P1], D)
    assert all(p[1] == 0 for p in prod) and [int(p[0]) for p in prod] == P and all(p[0].denominator == 1 for p in prod)
    N = P1[-1]  # constant term = product of the six roots (degree 6, sign +)
    zeta = (N[0] / q**3, N[1] / q**3)
    order = next(k for k in range(1, 13) if (lambda z: z == (Fr(1), Fr(0)))(
        (lambda k: [z := (Fr(1), Fr(0))] and [z := kmul(z, zeta, D) for _ in range(k)][-1])(k)))
    print(f"{r['row']:7s} K=Q(sqrt-{D}) q={q}: Nm(pi) = q^3 * ({zeta[0]} + {zeta[1]}*sqrt(-{D})), order {order}")
