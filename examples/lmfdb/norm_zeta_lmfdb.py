"""Nm_{E/K}(pi) = zeta * q^3 for the four LMFDB classes of Sec. 8 of the note (rows 12T23, 12T48, 12T77, 12T219).
Usage: python norm_zeta_lmfdb.py lmfdb_examples_in2.json lmfdb_examples_in.json > norm_zeta_lmfdb.out
Exact (PARI through cypari): for every imaginary quadratic subfield K of E = Q[T]/(P), factor P over K. When P is
the product of two sextic factors P1, P2 = conj(P1) over K, print zeta = P1(0)/q^3 and decide whether it is a root
of unity (and its order). P1(0) is the product of the six roots of P1, i.e. Nm_{E/K}(pi) for a K-embedding of E.
With ordinarity and geometric simplicity (verify_lmfdb*.json), a relation Nm_{E/K}(pi) = zeta q^3 with zeta a root
of unity gives angle rank at most 5, hence exactly 5 by Proposition 8.2 of the note.
K = Q[y]/(quadpoly(dK)), with y = (dK + sqrt(dK))/2 (PARI's quadpoly convention)."""
import json, sys
from cypari import pari

for fname in sys.argv[1:]:
    for r in json.load(open(fname)):
        q = r['q']
        P = pari.Pol(r['weil_poly'])  # coefficients T^12 .. T^0, variable x
        assert pari.polisirreducible(P) and pari.poldegree(P) == 12
        quads = [s[0] for s in pari.nfsubfields(pari.polredbest(P), 2)]
        imag = sorted({int(pari.nfdisc(g)) for g in quads if int(pari.nfdisc(g)) < 0}, reverse=True)
        print(f"{r['name']} (row {r['row']}), q = {q}; imaginary quadratic subfields of E: discriminants {imag}")
        for dK in imag:
            T = pari.quadpoly(dK, 'y')                       # y^2 - t*y + n
            t = -pari.polcoef(T, 1, 'y')
            fa = pari.nffactor(pari.nfinit(T), P)
            degs = [int(pari.poldegree(f)) for f in fa[0]]
            if degs != [6, 6]:
                print(f"  K = Q(sqrt({dK})): factor degrees {degs}; no relation of this form")
                continue
            P1, P2 = fa[0][0], fa[0][1]
            prod_ok = pari.liftall(P1 * P2) == P
            conj_ok = pari.liftall(P2) == pari.lift(pari.Mod(pari.subst(pari.liftall(P1), 'y', t - pari('y')), T))
            zeta = pari.Mod(pari.liftall(pari.polcoef(P1, 0)), T) / q**3
            one = pari.Mod(1, T)
            order = next((k for k in range(1, 25) if zeta**k == one), None)
            print(f"  K = Q(sqrt({dK})), y^2 - {t}*y + {pari.polcoef(T, 0, 'y')} = 0: P = P1 * P2: {bool(prod_ok)}; "
                  f"P2 = conj(P1): {bool(conj_ok)}; P1(0)/q^3 = {pari.lift(zeta)}; "
                  + (f"root of unity of order {order}" if order else "not a root of unity"))
