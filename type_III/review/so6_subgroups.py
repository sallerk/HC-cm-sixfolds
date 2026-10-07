"""
so6_subgroups.py  (review of Lemma lem:typeIII-rho: Hg(A)_C = SO(beta_Y) without Flapan)

Argument: Hg(A)_C is connected, reductive, contained in SO(beta_Y) (it commutes with D and
preserves psi), acts irreducibly on Y (End_Hg(V) = End^0(A) = D forces End_Hg(Y) = C), and has
finite centre (centre acts by scalars on Y and lies in O(6)), so it is semisimple.  This script
enumerates all semisimple complex Lie algebras of rank <= 3 (rank of so(6)) and their irreducible
representations of dimension 6 that are non-trivial on every simple factor, and reports which
of them carry an invariant SYMMETRIC form.  Expected: only A3 = sl(4) on Lambda^2 C^4, i.e. so(6).
Weyl dimension formula via positive coroots generated from Cartan matrices; Frobenius-Schur type of
a self-dual irrep: symmetric iff <lambda, 2 rho^vee> is even.
"""
import itertools
from fractions import Fraction


CARTAN = {
    "A1": [[2]],
    "A2": [[2, -1], [-1, 2]],
    "A3": [[2, -1, 0], [-1, 2, -1], [0, -1, 2]],
    "B2": [[2, -2], [-1, 2]],          # alpha_1 long, alpha_2 short (convention irrelevant here)
    "B3": [[2, -1, 0], [-1, 2, -2], [0, -1, 2]],
    "C3": [[2, -1, 0], [-1, 2, -1], [0, -2, 2]],
    "G2": [[2, -1], [-3, 2]],
}


def positive_roots(C):
    """positive roots in the basis of simple roots, from Cartan matrix C (C[i][j] = <alpha_i^vee? ...>)
    uses: alpha + alpha_j is a root iff p - q > 0 in the alpha_j string (standard algorithm)"""
    n = len(C)
    simple = [tuple(int(i == j) for j in range(n)) for i in range(n)]
    roots = set(simple)
    layer = list(simple)
    while layer:
        new = []
        for r in layer:
            for j in range(n):
                # <r, alpha_j^vee> = sum_i r_i C[i][j]  (C[i][j] = <alpha_i, alpha_j^vee>)
                pairing = sum(r[i] * C[i][j] for i in range(n))
                # p = largest k with r - k alpha_j a root
                p = 0
                while True:
                    cand = tuple(r[i] - (p + 1) * (i == j) for i in range(n))
                    if cand in roots:
                        p += 1
                    else:
                        break
                q = p - pairing
                if q > 0:
                    s = tuple(r[i] + (i == j) for i in range(n))
                    if s not in roots:
                        roots.add(s); new.append(s)
        layer = new
    return sorted(roots)


def pos_coroots(C):
    # coroots of R are the roots of the dual system, whose Cartan matrix is the transpose
    n = len(C)
    CT = [[C[j][i] for j in range(n)] for i in range(n)]
    return positive_roots(CT)


def weyl_dim(typ, lam):
    cor = pos_coroots(CARTAN[typ])
    num = 1; den = 1
    for a in cor:
        num *= sum(a[i] * (lam[i] + 1) for i in range(len(lam)))
        den *= sum(a[i] for i in range(len(lam)))
    assert num % den == 0
    return num // den


def dual_weight(typ, lam):
    if typ in ("A2", "A3"):
        return tuple(reversed(lam))
    return tuple(lam)     # w0 = -1 for A1, B2, B3, C3, G2


def fs_sign(typ, lam):
    """+1 symmetric, -1 alternating (self-dual lam)"""
    cor = pos_coroots(CARTAN[typ])
    tot = sum(sum(a[i] * lam[i] for i in range(len(lam))) for a in cor)
    return 1 if tot % 2 == 0 else -1


def irreps_upto(typ, maxdim):
    n = len(CARTAN[typ])
    out = []
    for lam in itertools.product(range(7), repeat=n):
        if not any(lam):
            continue
        d = weyl_dim(typ, lam)
        if d <= maxdim:
            out.append((lam, d))
    return out


def sanity():
    assert len(pos_coroots(CARTAN["A3"])) == 6 and len(pos_coroots(CARTAN["B3"])) == 9
    assert len(pos_coroots(CARTAN["G2"])) == 6 and len(pos_coroots(CARTAN["C3"])) == 9
    assert weyl_dim("A3", (0, 1, 0)) == 6 and fs_sign("A3", (0, 1, 0)) == 1      # so(6)
    assert weyl_dim("C3", (1, 0, 0)) == 6 and fs_sign("C3", (1, 0, 0)) == -1     # sp(6)
    assert weyl_dim("G2", (1, 0)) in (7, 14) and weyl_dim("G2", (0, 1)) in (7, 14)
    assert weyl_dim("B3", (0, 0, 1)) in (7, 8) and weyl_dim("B3", (1, 0, 0)) in (7, 8)
    assert weyl_dim("A1", (5,)) == 6 and fs_sign("A1", (5,)) == -1


if __name__ == "__main__":
    sanity()
    rank = {t: len(C) for t, C in CARTAN.items()}
    types = list(CARTAN)
    found = []
    # semisimple algebras = multisets of simple types with total rank <= 3
    for r in range(1, 4):
        for combo in itertools.combinations_with_replacement(types, r):
            if sum(rank[t] for t in combo) > 3:
                continue
            reps = [irreps_upto(t, 6) for t in combo]
            for choice in itertools.product(*reps):
                dim = 1
                for (lam, d) in choice:
                    dim *= d
                if dim != 6:
                    continue
                selfdual = all(dual_weight(t, lam) == lam for t, (lam, d) in zip(combo, choice))
                sign = 1
                for t, (lam, d) in zip(combo, choice):
                    sign *= fs_sign(t, lam)
                kind = ("not self-dual" if not selfdual else ("ORTHOGONAL" if sign == 1 else "symplectic"))
                found.append((combo, choice, kind))
                print(f"{'x'.join(combo):10s} {[lam for lam, d in choice]}  dim 6  -> {kind}")
    orth = [f for f in found if f[2] == "ORTHOGONAL"]
    print(f"\n6-dimensional faithful irreducible representations with an invariant symmetric form: {len(orth)}")
    for f in orth:
        print("   ", f[0], [lam for lam, d in f[1]])
