# alg.py -- own exact arithmetic for F = prod_i Q[y]/(f_i) (f_i monic integral, totally real),
# E = K (x) F with K = Q(sqrt(-d)).  Used for independent checks (no PARI here).
from fractions import Fraction
import mpmath

mpmath.mp.dps = 80


def poly_mulmod(a, b, f):
    """a, b: coefficient lists (low degree first) of length m; f monic of degree m (list len m+1)."""
    m = len(f) - 1
    prod = [Fraction(0)] * (2 * m - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0:
                continue
            prod[i + j] += ai * bj
    # reduce mod f
    for k in range(2 * m - 2, m - 1, -1):
        c = prod[k]
        if c != 0:
            prod[k] = Fraction(0)
            for t in range(m):
                prod[k - m + t] -= c * f[t]
    return prod[:m]


def det_frac(M):
    """Exact determinant of a square matrix of Fractions (Gaussian elimination)."""
    n = len(M)
    A = [[Fraction(x) for x in row] for row in M]
    det = Fraction(1)
    for c in range(n):
        piv = None
        for r in range(c, n):
            if A[r][c] != 0:
                piv = r
                break
        if piv is None:
            return Fraction(0)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = -det
        det *= A[c][c]
        inv = 1 / A[c][c]
        for r in range(c + 1, n):
            if A[r][c] != 0:
                fac = A[r][c] * inv
                for k in range(c, n):
                    A[r][k] -= fac * A[c][k]
    return det


class Factor:
    def __init__(self, f):
        # f: list of ints, low degree first, monic
        self.f = [Fraction(c) for c in f]
        self.m = len(f) - 1
        assert self.f[-1] == 1
        self.roots = self._real_roots()

    def _real_roots(self):
        m = self.m
        if m == 1:
            c = -self.f[0]
            return [mpmath.mpf(c.numerator) / c.denominator]
        coeffs = [mpmath.mpf(int(c)) if c.denominator == 1 else mpmath.mpf(c.numerator) / c.denominator
                  for c in reversed(self.f)]
        rts = mpmath.polyroots(coeffs, maxsteps=500, extraprec=400)
        out = []
        for r in rts:
            if abs(mpmath.im(r)) > mpmath.mpf(10) ** (-40):
                raise ValueError("not totally real")
            out.append(mpmath.re(r))
        return sorted(out)

    def one(self):
        return [Fraction(1)] + [Fraction(0)] * (self.m - 1)

    def elt(self, coeffs):
        c = [Fraction(x) for x in coeffs] + [Fraction(0)] * (self.m - len(coeffs))
        return c[:self.m]

    def mul(self, a, b):
        return poly_mulmod(a, b, self.f)

    def mulmat(self, a):
        """matrix of multiplication by a on the power basis (columns = images of y^k)."""
        m = self.m
        cols = []
        e = self.one()
        for k in range(m):
            basis = [Fraction(0)] * m
            basis[k] = Fraction(1)
            cols.append(self.mul(a, basis))
        return [[cols[k][r] for k in range(m)] for r in range(m)]

    def norm(self, a):
        return det_frac(self.mulmat(a))

    def trace(self, a):
        M = self.mulmat(a)
        return sum(M[i][i] for i in range(self.m))

    def evalf(self, a, r):
        s = mpmath.mpf(0)
        for k in range(self.m - 1, -1, -1):
            c = a[k]
            s = s * r + (mpmath.mpf(c.numerator) / c.denominator)
        return s

    def signs(self, a):
        out = []
        for r in self.roots:
            v = self.evalf(a, r)
            if abs(v) < mpmath.mpf(10) ** (-30):
                raise ValueError("sign undecided")
            out.append(1 if v > 0 else -1)
        return out

    def trace_form_disc(self):
        """det of (Tr(y^j y^k)) = discriminant of the power basis (exact)."""
        m = self.m
        pw = [self.one()]
        for k in range(1, 2 * m - 1):
            pw.append(self.mul(pw[-1], self.elt([0, 1]) if m > 1 else self.elt([0])))
        tr = [self.trace(p) for p in pw]
        return det_frac([[tr[j + k] for k in range(m)] for j in range(m)])


class Alg:
    """F = prod of factors; elements are lists of per-factor coefficient lists."""
    def __init__(self, polys):
        self.factors = [Factor(f) for f in polys]
        self.deg = sum(F.m for F in self.factors)
        assert self.deg % 2 == 0
        self.n = self.deg // 2

    def places(self):
        """list of (factor index, root index) -- the real places, in a fixed order."""
        return [(i, j) for i, F in enumerate(self.factors) for j in range(F.m)]

    def norm(self, a):
        r = Fraction(1)
        for F, ai in zip(self.factors, a):
            r *= F.norm(ai)
        return r

    def signs(self, a):
        out = []
        for F, ai in zip(self.factors, a):
            out += F.signs(ai)
        return out

    def mul(self, a, b):
        return [F.mul(x, y) for F, x, y in zip(self.factors, a, b)]

    def one(self):
        return [F.one() for F in self.factors]

    def disc(self):
        r = Fraction(1)
        for F in self.factors:
            r *= F.trace_form_disc()
        return r

    def gram_ha(self, a):
        """Gram matrix (rational, symmetric) of h_a(x,y) = Tr_{E/K}(a x ybar) on the K-basis given by
        the power bases of the factors: entries Tr_{F_i/Q}(a_i y^j y^k), block diagonal."""
        N = self.deg
        G = [[Fraction(0)] * N for _ in range(N)]
        off = 0
        for F, ai in zip(self.factors, a):
            m = F.m
            # a * y^t for t = 0..2m-2
            pw = [ai]
            ygen = F.elt([0, 1]) if m > 1 else None
            for t in range(1, 2 * m - 1):
                pw.append(F.mul(pw[-1], ygen))
            tr = [F.trace(p) for p in pw]
            for j in range(m):
                for k in range(m):
                    G[off + j][off + k] = tr[j + k]
            off += m
        return G
