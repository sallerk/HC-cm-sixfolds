# Which Fermat-curve Jacobian factors (types H_{a,b,c}, a+b+c = 0 mod m) have reduced CM field Q(sqrt(-11))
# (i.e. are isogenous to a power of E_{sqrt-11})? Also locate J_11 = H_{1,1,9} mod 11.
from math import gcd
def H(a, b, c, m):
    return frozenset(t for t in range(1, m) if gcd(t, m) == 1 and ((t*a) % m + (t*b) % m + (t*c) % m) == m)
def qsqrtm11_subgroup(m):
    # elements t of (Z/m)^x fixing sqrt(-11) in Q(zeta_m) (11 | m): t mod 11 is a square mod 11
    sq = {(x*x) % 11 for x in range(1, 11)}
    return frozenset(t for t in range(1, m) if gcd(t, m) == 1 and t % 11 in sq)
for m in (11, 22, 33, 44, 55, 66, 77, 88, 99, 121):
    U = [t for t in range(1, m) if gcd(t, m) == 1]
    W = qsqrtm11_subgroup(m)
    hits = set()
    for a in range(1, m):
        for b in range(1, m):
            c = (-a - b) % m
            if c == 0: continue
            if gcd(gcd(a, b), gcd(c, m)) != 1: continue   # primitive part only (d = m)
            h = H(a, b, c, m)
            stab = frozenset(w for w in U if frozenset((w*t) % m for t in h) == h)
            if stab == W:
                hits.add(tuple(sorted((a, b, c))))
    print(f'm={m}: primitive types with reduced CM field Q(sqrt-11): {len(hits)}', sorted(hits)[:6])
print('J_11 type H_{1,1,9} mod 11 =', sorted(H(1, 1, 9, 11)))
