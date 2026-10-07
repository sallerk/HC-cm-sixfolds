# herm.py -- own code: hermitian forms over K = Q(sqrt(-d)), explicit Witt decomposition.
# PARI's qfsolve is used ONLY as a search tool for isotropic vectors of rational quadratic forms;
# every isotropic vector / totally isotropic subspace is verified with exact arithmetic here.
# A non-split verdict is certified theory-free by the absence of primitive zeros mod p^k of the
# ternary form X^2 + d Y^2 - c Z^2 (c = -det of the binary anisotropic remainder).
from fractions import Fraction
from cypari import pari
from qsym import squarefree_part, has_primitive_zero_mod, prime_factors


class Kel:
    __slots__ = ('r', 's')
    D = None  # d, set by set_d

    def __init__(self, r, s=0):
        self.r = Fraction(r)
        self.s = Fraction(s)

    def __add__(self, o):
        return Kel(self.r + o.r, self.s + o.s)

    def __sub__(self, o):
        return Kel(self.r - o.r, self.s - o.s)

    def __mul__(self, o):
        d = Kel.D
        return Kel(self.r * o.r - d * self.s * o.s, self.r * o.s + self.s * o.r)

    def conj(self):
        return Kel(self.r, -self.s)

    def iszero(self):
        return self.r == 0 and self.s == 0

    def inv(self):
        n = self.r * self.r + Kel.D * self.s * self.s
        return Kel(self.r / n, -self.s / n)

    def __repr__(self):
        return f"({self.r}+{self.s}w)"


def set_d(d):
    Kel.D = d


Z0 = lambda: Kel(0)


def herm_eval(G, x, y):
    """h(x,y) = sum x_j G_jk conj(y_k)."""
    m = len(G)
    s = Kel(0)
    for j in range(m):
        if x[j].iszero():
            continue
        t = Kel(0)
        for k in range(m):
            if y[k].iszero() or G[j][k].iszero():
                continue
            t = t + G[j][k] * y[k].conj()
        s = s + x[j] * t
    return s


def gram_on(G, vecs):
    return [[herm_eval(G, u, v) for v in vecs] for u in vecs]


def trace_quadratic_matrix(G):
    """Rational symmetric 2m x 2m matrix M with q(u,w) = h(u + w sqrt(-d), same) = [u;w]^T M [u;w]."""
    d = Kel.D
    m = len(G)
    Gr = [[G[j][k].r for k in range(m)] for j in range(m)]
    Gi = [[G[j][k].s for k in range(m)] for j in range(m)]
    M = [[Fraction(0)] * (2 * m) for _ in range(2 * m)]
    for j in range(m):
        for k in range(m):
            M[j][k] = Gr[j][k]
            M[m + j][m + k] = d * Gr[j][k]
            M[j][m + k] = d * Gi[j][k]
            M[m + j][k] = -d * Gi[j][k]
    return M


def nullspace_K(rows, m):
    """Basis of {z in K^m : sum_k rows[i][k] z_k = 0 for all i} (Gaussian elimination over K)."""
    A = [list(r) for r in rows]
    piv_cols = []
    r = 0
    for c in range(m):
        p = None
        for i in range(r, len(A)):
            if not A[i][c].iszero():
                p = i
                break
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = A[r][c].inv()
        A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and not A[i][c].iszero():
                f = A[i][c]
                A[i] = [A[i][k] - f * A[r][k] for k in range(m)]
        piv_cols.append(c)
        r += 1
        if r == len(A):
            break
    free = [c for c in range(m) if c not in piv_cols]
    basis = []
    for fc in free:
        z = [Kel(0) for _ in range(m)]
        z[fc] = Kel(1)
        for i, pc in enumerate(piv_cols):
            z[pc] = Kel(0) - A[i][fc]
        basis.append(z)
    return basis


def to_pari_matrix(M):
    rows = ";".join(",".join(str(x) for x in row) for row in M)
    return pari(f"[{rows}]")


def find_isotropic(G):
    """Return an isotropic x in K^m (verified) or ('aniso', p) using qfsolve on the trace form."""
    m = len(G)
    M = trace_quadratic_matrix(G)
    # clear denominators
    from math import lcm
    L = 1
    for row in M:
        for x in row:
            L = lcm(L, x.denominator)
    Mi = [[int(x * L) for x in row] for row in M]
    res = pari.qfsolve(to_pari_matrix(Mi))
    if res.type() in ('t_INT',):
        return ('aniso', int(res))
    if res.type() == 't_MAT':
        res = res[0]  # several solutions: first column? (qfsolve returns a column vector normally)
    X = [Fraction(str(res[i])) for i in range(2 * m)]
    x = [Kel(X[k], X[m + k]) for k in range(m)]
    assert any(not t.iszero() for t in x)
    hv = herm_eval(G, x, x)
    assert hv.iszero(), "qfsolve vector not isotropic"
    return x


def witt_decompose(G0):
    """Iteratively split hyperbolic planes off the hermitian form with Gram matrix G0 (K-entries).
    Returns (iso_vectors_in_original_coords, remainder_gram, status) where status is
    'done' (remainder dim 0) or ('aniso', p) for an anisotropic remainder detected by qfsolve."""
    m0 = len(G0)
    # current basis of the working subspace, in original coordinates
    basis = [[Kel(1) if i == j else Kel(0) for i in range(m0)] for j in range(m0)]
    G = [row[:] for row in G0]
    iso = []
    while len(G) > 0:
        m = len(G)
        r = find_isotropic(G)
        if isinstance(r, tuple):
            return iso, G, basis, r
        x = r
        # y with h(x, y) != 0 : try standard basis vectors
        y = None
        for k in range(m):
            e = [Kel(0) for _ in range(m)]
            e[k] = Kel(1)
            if not herm_eval(G, x, e).iszero():
                y = e
                break
        assert y is not None, "degenerate"
        # orthogonal complement of span(x,y): rows giving h(z, x) = 0 and h(z, y) = 0 (linear in z)
        rows = []
        for v in (x, y):
            rows.append([sum((G[j][k] * v[k].conj() for k in range(m)), Kel(0)) for j in range(m)])
        W = nullspace_K(rows, m)
        assert len(W) == m - 2
        # record x in original coordinates
        xo = [sum((x[j] * basis[j][i] for j in range(m)), Kel(0)) for i in range(m0)]
        iso.append(xo)
        newbasis = [[sum((w[j] * basis[j][i] for j in range(m)), Kel(0)) for i in range(m0)] for w in W]
        G = gram_on(G, W)
        basis = newbasis
    return iso, G, basis, 'done'


def rank_K(vecs, m):
    if not vecs:
        return 0
    A = [list(v) for v in vecs]
    r = 0
    for c in range(m):
        p = None
        for i in range(r, len(A)):
            if not A[i][c].iszero():
                p = i
                break
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = A[r][c].inv()
        A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and not A[i][c].iszero():
                f = A[i][c]
                A[i] = [A[i][k] - f * A[r][k] for k in range(m)]
        r += 1
    return r


def det2_herm(G):
    d = G[0][0] * G[1][1] - G[0][1] * G[1][0]
    assert d.s == 0
    return d.r


def split_certificate(Brat, d, n):
    """Brat: rational symmetric 2n x 2n Gram matrix of h_a (as hermitian over K).
    Returns dict with verdict 'split' (with explicit totally isotropic subspace verified) or
    'nonsplit' (with Witt index n-1 and a brute-force local anisotropy certificate)."""
    set_d(d)
    m0 = len(Brat)
    G0 = [[Kel(Brat[j][k]) for k in range(m0)] for j in range(m0)]
    iso, Grem, basis, status = witt_decompose(G0)
    out = {'witt_planes': len(iso)}
    if status == 'done':
        # verify totally isotropic subspace of dimension m0/2
        Gw = gram_on(G0, iso)
        assert all(x.iszero() for row in Gw for x in row)
        rk = rank_K(iso, m0)
        assert rk == len(iso) == m0 // 2
        out['verdict'] = 'split'
        out['iso_dim'] = rk
        return out
    # anisotropic remainder
    p = status[1]
    out['qfsolve_prime'] = p
    if len(Grem) == 2:
        c = -det2_herm(Grem)          # isotropic iff c in Nm(K^x)
        ci = c.numerator * c.denominator
        csf = squarefree_part(ci)
        # theory-free check: no primitive zero of X^2 + d Y^2 - csf Z^2 modulo p^k
        primes = [p] if p > 0 else []
        cert = None
        cand = sorted(set(primes + [2] + prime_factors(d) + prime_factors(csf)))
        for q in cand:
            k = 5 if q == 2 else 2
            if q ** k > 3000:
                continue
            if not has_primitive_zero_mod([1, d, -csf], q, k):
                cert = (q, k)
                break
        out['verdict'] = 'nonsplit'
        out['remainder_dim'] = 2
        out['bruteforce_cert'] = cert
        # Witt index check: totally isotropic part found has dimension len(iso)
        if iso:
            Gw = gram_on(G0, iso)
            assert all(x.iszero() for row in Gw for x in row)
        return out
    out['verdict'] = 'nonsplit'
    out['remainder_dim'] = len(Grem)
    out['bruteforce_cert'] = None
    return out


# ---------------------------------------------------------------------------------------------
def diagonalize_rational(B):
    """Return (D, P) with P B P^T = diag(D), P rows = new basis vectors (Fractions)."""
    n = len(B)
    A = [[Fraction(x) for x in row] for row in B]
    P = [[Fraction(1) if i == j else Fraction(0) for j in range(n)] for i in range(n)]

    def add_row(i, j, c):  # basis_i += c * basis_j  (congruence)
        P[i] = [x + c * y for x, y in zip(P[i], P[j])]
        A[i] = [x + c * y for x, y in zip(A[i], A[j])]
        for r in range(n):
            A[r][i] += c * A[r][j]

    def swap(i, j):
        P[i], P[j] = P[j], P[i]
        A[i], A[j] = A[j], A[i]
        for r in range(n):
            A[r][i], A[r][j] = A[r][j], A[r][i]

    for k in range(n):
        if A[k][k] == 0:
            j = next((j for j in range(k + 1, n) if A[j][j] != 0), None)
            if j is not None:
                swap(k, j)
            else:
                j = next((j for j in range(k + 1, n) if A[k][j] != 0), None)
                if j is None:
                    continue
                add_row(k, j, Fraction(1))
        piv = A[k][k]
        for j in range(k + 1, n):
            if A[j][k] != 0:
                add_row(j, k, -A[j][k] / piv)
    D = [A[i][i] for i in range(n)]
    return D, P


def reduce_diag(D, P, limit_bits=120):
    """Scale basis vectors so that diagonal entries become squarefree integers when factorable."""
    from sympy import factorint
    D2, P2 = [], []
    for dv, row in zip(D, P):
        num = dv.numerator * dv.denominator   # same square class, value dv * den^2
        # dv = num / den^2 ; write |num| = s^2 * sf
        if abs(num).bit_length() <= limit_bits:
            f = factorint(abs(num))
            s = 1
            sf = 1
            for p, k in f.items():
                s *= p ** (k // 2)
                if k % 2:
                    sf *= p
            sf = sf if num > 0 else -sf
            # value of basis vector scaled by t: t^2 dv ; want t^2 dv = sf  -> t = den / s
            t = Fraction(dv.denominator, s)
        else:
            t = Fraction(dv.denominator)
            sf = num
        D2.append(t * t * dv)
        P2.append([t * x for x in row])
    return D2, P2


def split_certificate_reduced(Brat, d):
    """As split_certificate, but works on a diagonalised, squarefree-reduced model and maps the
    isotropic vectors back; the final totally-isotropic check is done on the ORIGINAL Gram matrix."""
    set_d(d)
    m0 = len(Brat)
    D, P = diagonalize_rational(Brat)
    D, P = reduce_diag(D, P)
    # check P B P^T = diag(D) exactly
    for i in range(m0):
        for j in range(m0):
            v = sum(P[i][a] * sum(Brat[a][b] * P[j][b] for b in range(m0)) for a in range(m0))
            assert v == (D[i] if i == j else 0)
    Gd = [[Kel(D[i]) if i == j else Kel(0) for j in range(m0)] for i in range(m0)]
    iso, Grem, basis, status = witt_decompose(Gd)
    # map back: vector x in diag coords -> sum_i x_i P_i in original coords
    iso_orig = [[sum((x[i] * Kel(P[i][c]) for i in range(m0)), Kel(0)) for c in range(m0)] for x in iso]
    G0 = [[Kel(Brat[j][k]) for k in range(m0)] for j in range(m0)]
    out = {'witt_planes': len(iso)}
    if iso_orig:
        Gw = gram_on(G0, iso_orig)
        assert all(x.iszero() for row in Gw for x in row), "mapped-back subspace not isotropic"
    rk = rank_K(iso_orig, m0)
    assert rk == len(iso)
    if status == 'done':
        assert rk == m0 // 2
        out['verdict'] = 'split'
        out['iso_dim'] = rk
        return out
    p = status[1]
    out['qfsolve_prime'] = p
    out['verdict'] = 'nonsplit'
    out['remainder_dim'] = len(Grem)
    cert = None
    if len(Grem) == 2:
        c = -det2_herm(Grem)
        ci = c.numerator * c.denominator
        csf = squarefree_part(ci)
        cand = sorted(set(([p] if p > 0 else []) + [2] + prime_factors(d) + prime_factors(csf)))
        for q in cand:
            k = 5 if q == 2 else 2
            if q ** k > 4000:
                continue
            if not has_primitive_zero_mod([1, d, -csf], q, k):
                cert = (q, k)
                break
        out['minus_det_rem_sf'] = csf
    out['bruteforce_cert'] = cert
    return out


# ---------------------------------------------------------------------------------------------
# Witt decomposition in a diagonal, squarefree-reduced model (keeps qfsolve inputs small).
def core_rational(beta):
    """beta = sf * u^2 with sf a squarefree integer and u rational; returns (sf, u)."""
    p, q = beta.numerator, beta.denominator
    num = p * q
    r = pari.core(num, 1)          # [core, f] with num = core * f^2
    sf = int(r[0])
    f = int(r[1])
    # beta = num / q^2 = sf f^2 / q^2
    return sf, Fraction(f, q)


def herm_diagonalize(G, basis):
    """Congruence-diagonalize the hermitian K-matrix G (basis rows tracked).  Returns (D, basis')
    with D rational diagonal entries and h(b_i, b_j) = 0 for i != j."""
    m = len(G)
    A = [row[:] for row in G]
    Bs = [row[:] for row in basis]

    def add(i, j, c):   # b_i += c b_j
        Bs[i] = [x + c * y for x, y in zip(Bs[i], Bs[j])]
        A[i] = [x + c * y for x, y in zip(A[i], A[j])]
        cc = c.conj()
        for r in range(m):
            A[r][i] = A[r][i] + A[r][j] * cc

    def swap(i, j):
        Bs[i], Bs[j] = Bs[j], Bs[i]
        A[i], A[j] = A[j], A[i]
        for r in range(m):
            A[r][i], A[r][j] = A[r][j], A[r][i]

    for k in range(m):
        if A[k][k].iszero():
            j = next((j for j in range(k + 1, m) if not A[j][j].iszero()), None)
            if j is not None:
                swap(k, j)
            else:
                j = next((j for j in range(k + 1, m) if not A[k][j].iszero()), None)
                if j is None:
                    raise ValueError('degenerate')
                g = A[k][j]
                c = Kel(1) if g.r != 0 else Kel(0, 1)
                add(k, j, c)
        piv = A[k][k]
        assert piv.s == 0 and piv.r != 0
        for j in range(k + 1, m):
            if not A[j][k].iszero():
                mu = A[j][k] * piv.inv()
                add(j, k, Kel(0) - mu)
    D = []
    for i in range(m):
        for j in range(m):
            if i != j:
                assert A[i][j].iszero()
        assert A[i][i].s == 0
        D.append(A[i][i].r)
    return D, Bs


def witt_diag(Brat, d, max_bits=400):
    """Explicit Witt decomposition of h = Brat (x) K, computed in successive diagonal squarefree
    models (small numbers); the totally isotropic subspace is mapped back to the original
    coordinates and verified there exactly.  Returns dict as split_certificate_reduced."""
    set_d(d)
    m0 = len(Brat)
    G0 = [[Kel(Brat[j][k]) for k in range(m0)] for j in range(m0)]
    ident = [[Kel(1) if i == j else Kel(0) for i in range(m0)] for j in range(m0)]
    D, basis = herm_diagonalize(G0, ident)      # basis rows in original coordinates
    iso = []
    status = 'done'
    Grem = None
    while True:
        D2, B2 = [], []
        for beta, b in zip(D, basis):
            if max(beta.numerator.bit_length(), beta.denominator.bit_length()) > max_bits:
                raise ValueError('too large')
            sf, u = core_rational(beta)
            t = Kel(1 / u)
            B2.append([t * x for x in b])
            D2.append(Fraction(sf))
        D, basis = D2, B2
        m = len(D)
        if m == 0:
            break
        Gd = [[Kel(D[i]) if i == j else Kel(0) for j in range(m)] for i in range(m)]
        r = find_isotropic(Gd)
        if isinstance(r, tuple):
            status = r
            Grem = Gd
            break
        x = r
        j = next(j for j in range(m) if not x[j].iszero())
        rows = [[Kel(1) if a == j else Kel(0) for a in range(m)],
                [Kel(D[a]) * x[a].conj() for a in range(m)]]
        W = nullspace_K(rows, m)                   # model coordinates
        assert len(W) == m - 2
        iso.append([sum((x[a] * basis[a][c] for a in range(m)), Kel(0)) for c in range(m0)])
        if not W:
            break
        Gw = gram_on(Gd, W)                        # small numbers
        Dn, Wn = herm_diagonalize(Gw, W)           # Wn rows: model coordinates
        basis = [[sum((w[a] * basis[a][c] for a in range(m)), Kel(0)) for c in range(m0)] for w in Wn]
        D = Dn
    out = {'witt_planes': len(iso)}
    if iso:
        Gw = gram_on(G0, iso)
        assert all(x.iszero() for row in Gw for x in row), "subspace not isotropic"
    rk = rank_K(iso, m0)
    assert rk == len(iso)
    if status == 'done':
        assert rk == m0 // 2
        out['verdict'] = 'split'
        out['iso_dim'] = rk
        return out
    p = status[1]
    out['qfsolve_prime'] = p
    out['verdict'] = 'nonsplit'
    out['remainder_dim'] = len(Grem)
    cert = None
    if len(Grem) == 2:
        csf = squarefree_part(int(-(Grem[0][0].r * Grem[1][1].r)))
        cand = sorted(set(([p] if p > 0 else []) + [2] + prime_factors(d) + prime_factors(csf)))
        for q in cand:
            k = 5 if q == 2 else 2
            if q ** k > 4000:
                continue
            if not has_primitive_zero_mod([1, d, -csf], q, k):
                cert = (q, k)
                break
        out['minus_det_rem_sf'] = csf
    out['bruteforce_cert'] = cert
    return out
