# inst.py -- random test instances for Theorem A: totally real algebras F = prod F_i, d, sign patterns.
import random
from cypari import pari


def coeffs_of(polstr_or_gen):
    """PARI polynomial (variable x or y) -> list of ints, low degree first."""
    P = pari(polstr_or_gen) if isinstance(polstr_or_gen, str) else polstr_or_gen
    P = pari.substpol(P, pari('x'), pari('y')) if P.variable() == pari('x') else P
    deg = int(pari.poldegree(P))
    return [int(pari.polcoef(P, k)) for k in range(deg + 1)]


def poly_y(coeffs):
    return " + ".join(f"({c})*y^{k}" for k, c in enumerate(coeffs))


def random_totally_real(m, rng, R=3):
    """Monic integral polynomial (coeff list, low first) of a random totally real field of degree m."""
    if m == 1:
        return [0, 1]
    while True:
        M = [[0] * m for _ in range(m)]
        for i in range(m):
            for j in range(i, m):
                v = rng.randint(-R, R)
                M[i][j] = M[j][i] = v
        rows = ";".join(",".join(str(x) for x in r) for r in M)
        P = pari(f"charpoly([{rows}], y)")
        if not pari.polisirreducible(P):
            continue
        P = pari.polredbest(P)
        P = pari.substpol(P, pari.variable(P), pari('y'))
        c = [int(pari.polcoef(P, k)) for k in range(m + 1)]
        # sanity: totally real
        nf = pari.nfinit(P)
        if int(nf[1][0]) != m:
            continue
        return c


def random_real_quadratic(rng, B=200):
    while True:
        k = rng.randint(2, B)
        if pari.issquarefree(k):
            return [-k, 0, 1]


def random_squarefree_d(rng, B=120):
    while True:
        d = rng.randint(1, B)
        if pari.issquarefree(d):
            return d


SHAPES = {
    'field2': [2], 'field4': [4], 'field6': [6], 'field8': [8],
    '1+1': [1, 1], '2+2': [2, 2], '1+3': [1, 3], '1+1+2': [1, 1, 2], '3+3': [3, 3], '2+4': [2, 4],
    '2+2+2': [2, 2, 2], '1+5': [1, 5], '1+1+4': [1, 1, 4], '1+1+1+1': [1, 1, 1, 1],
    '2+2+2+2': [2, 2, 2, 2], '4+4': [4, 4], '3+5': [3, 5], '1+1+3+3': [1, 1, 3, 3],
}


def random_instance(shape, rng):
    polys = []
    for m in SHAPES[shape]:
        if m == 2 and rng.random() < 0.5:
            polys.append(random_real_quadratic(rng))
        else:
            polys.append(random_totally_real(m, rng))
    d = random_squarefree_d(rng)
    return polys, d


def random_pattern(total, n, rng):
    """sign list of length total with exactly n entries -1."""
    idx = rng.sample(range(total), n)
    return [(-1 if i in idx else 1) for i in range(total)]
