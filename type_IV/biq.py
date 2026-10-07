# biq.py -- shared exact arithmetic for the biquadratic type-IV split lemma (note, Sec. 7).
#
# Setting: K = Q(s), s^2 = -d (d > 0 squarefree); F = Q(w), w^2 = m (m > 1 squarefree); L = K F.
# Elements of L: 4-tuples of Fractions (a0, a1, a2, a3) = a0 + a1 s + a2 w + a3 s w.
# Embeddings over the same tau: K -> C (s -> i sqrt d):  sigma1: w -> +sqrt m (real place f1),
#                                                       sigma2: w -> -sqrt m (real place f2).
# H: 3x3 L-Hermitian matrix (H[j][i] = conj(H[i][j])), H(x,y) = sum x_i H_ij conj(y_j) (linear in x).
# h_a(x,y) = Tr_{L/K}(a H(x,y)), a in F; K-basis of V = L^3: b_(i,j) = w^j e_i (i=0..2, j=0..1).
from fractions import Fraction as Fr
import math, cmath
import numpy as np
from cypari import pari
try:
    pari.allocatemem(2 * 10**9, silent=True)
except TypeError:
    pari.allocatemem(2 * 10**9)

class Ctx:
    def __init__(self, d, m):
        self.d, self.m = d, m
        self.dF = int(pari.nfdisc(pari('x^2-%d' % m)))
        self._bnf = None; self._T = None

    # ---------- L arithmetic ----------
    def mul(self, x, y):
        d, m = self.d, self.m
        a0, a1, a2, a3 = x; b0, b1, b2, b3 = y
        return (a0*b0 - d*a1*b1 + m*a2*b2 - d*m*a3*b3,
                a0*b1 + a1*b0 + m*a2*b3 + m*a3*b2,
                a0*b2 + a2*b0 - d*a1*b3 - d*a3*b1,
                a0*b3 + a3*b0 + a1*b2 + a2*b1)
    @staticmethod
    def add(x, y): return tuple(a + b for a, b in zip(x, y))
    @staticmethod
    def sub(x, y): return tuple(a - b for a, b in zip(x, y))
    @staticmethod
    def neg(x): return tuple(-a for a in x)
    @staticmethod
    def conj(x): return (x[0], -x[1], x[2], -x[3])          # complex conjugation (s -> -s)
    @staticmethod
    def rho(x): return (x[0], x[1], -x[2], -x[3])           # Gal(L/K) generator (w -> -w)
    @staticmethod
    def L(a0=0, a1=0, a2=0, a3=0): return (Fr(a0), Fr(a1), Fr(a2), Fr(a3))
    @staticmethod
    def Fel(x, y=0): return (Fr(x), Fr(0), Fr(y), Fr(0))    # x + y w in F
    def normQ(self, x):
        p = self.mul(self.mul(x, self.conj(x)), self.mul(self.rho(x), self.rho(self.conj(x))))
        assert p[1] == p[2] == p[3] == 0
        return p[0]
    def inv(self, x):
        n = self.normQ(x)
        assert n != 0
        y = self.mul(self.conj(x), self.mul(self.rho(x), self.rho(self.conj(x))))
        return tuple(c / n for c in y)
    def div(self, x, y): return self.mul(x, self.inv(y))
    def NF(self, f):                                         # N_{F/Q} of f in F
        assert f[1] == 0 and f[3] == 0
        return f[0]*f[0] - self.m*f[2]*f[2]
    def trLK(self, x): return (2*x[0], 2*x[1])               # Tr_{L/K}: K element (re, im-coeff of s)
    def emb(self, x, k):                                     # complex embedding sigma_k, k = 1, 2
        sd, sm = math.sqrt(self.d), math.sqrt(self.m) * (1 if k == 1 else -1)
        return complex(float(x[0]) + float(x[2])*sm, float(x[1])*sd + float(x[3])*sd*sm)
    def realF(self, f, k):
        return float(f[0]) + float(f[2]) * math.sqrt(self.m) * (1 if k == 1 else -1)
    def totpos(self, f):
        return self.realF(f, 1) > 0 and self.realF(f, 2) > 0 and self.NF(f) > 0

    # ---------- Hermitian matrices over L ----------
    def Hform(self, H, x, y):
        z = self.L()
        for i in range(3):
            for j in range(3):
                z = self.add(z, self.mul(self.mul(x[i], H[i][j]), self.conj(y[j])))
        return z
    def is_herm(self, H):
        return all(H[j][i] == self.conj(H[i][j]) for i in range(3) for j in range(3))
    def signature(self, H, k):
        A = np.array([[self.emb(H[i][j], k) for j in range(3)] for i in range(3)])
        ev = np.linalg.eigvalsh(A)
        return (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()))
    def detL(self, H):
        M = self.mul; S = self.sub; A = self.add
        t1 = M(H[0][0], S(M(H[1][1], H[2][2]), M(H[1][2], H[2][1])))
        t2 = M(H[0][1], S(M(H[1][0], H[2][2]), M(H[1][2], H[2][0])))
        t3 = M(H[0][2], S(M(H[1][0], H[2][1]), M(H[1][1], H[2][0])))
        D = A(S(t1, t2), t3)
        assert D[1] == 0 and D[3] == 0, D                  # det of Hermitian matrix lies in F
        return D

    # ---------- K-Hermitian form h_a and its invariants ----------
    def gram_K(self, H, a):
        """6x6 Gram matrix of h_a on the K-basis w^j e_i, as list of K-elements (re, s-coeff)."""
        W = [self.L(1), self.L(0, 0, 1)]
        G = [[None]*6 for _ in range(6)]
        for i in range(3):
            for j in range(2):
                for k in range(3):
                    for l in range(2):
                        z = self.mul(a, self.mul(self.mul(W[j], W[l]), H[i][k]))   # conj(w^l) = w^l
                        G[2*i + j][2*k + l] = self.trLK(z)
        return G
    def detK(self, G):
        d = self.d
        rows = []
        for r in G:
            rows.append(','.join('Mod(%s+(%s)*t,t^2+%d)' % (str(x), str(y), d) for x, y in r))
        D = pari('matdet([' + ';'.join(rows) + '])')
        D = pari.lift(D)
        c1 = pari.polcoef(D, 1, 't'); c0 = pari.polcoef(D, 0, 't')
        assert int(c1 == 0), D
        return Fr(int(pari.numerator(c0)), int(pari.denominator(c0)))
    def gram_Q(self, G):
        """12x12 rational Gram matrix of the quadratic form x -> h(x,x) on the Q-basis {b, s b}."""
        d = self.d
        n = len(G)
        Q = [[Fr(0)]*(2*n) for _ in range(2*n)]
        for p in range(n):
            for q in range(n):
                re, im = G[p][q]                       # h(b_p, b_q) = re + im s
                # h(s^e b_p, s^f b_q) = s^e conj(s)^f h(b_p,b_q);  rational part of that:
                for e in range(2):
                    for f in range(2):
                        # multiply (re + im s) by s^e (-s)^f
                        x, y = re, im
                        for _ in range(e): x, y = -d*y, x          # *s
                        for _ in range(f): x, y = d*y, -x          # *(-s)
                        Q[2*p + e][2*q + f] = x                    # Re_K part = rational part
        return Q

    # ---------- local arithmetic ----------
    def hilbert(self, a, b, p):
        a = pari(a.numerator) / pari(a.denominator) if isinstance(a, Fr) else pari(a)
        b = pari(b.numerator) / pari(b.denominator) if isinstance(b, Fr) else pari(b)
        return int(pari.hilbert(a, b, p))
    def bad_primes(self, *xs):
        N = 2 * self.d * self.m
        for x in xs:
            x = Fr(x); N *= x.numerator * x.denominator
        return sorted(set(abs(int(p)) for p in pari.factor(abs(N))[0]) - {1})
    def is_norm_K(self, x):
        x = Fr(x)
        if x <= 0: return False
        return all(self.hilbert(x, -self.d, p) == 1 for p in self.bad_primes(x))
    def is_square_Qp(self, a, p):
        a = Fr(a)
        return bool(pari('issquare(%s+O(%d^40))' % (str(a), p)))
    def forced_primes(self, x):
        """primes p (among those where a symbol can be nontrivial) with -dm in Q_p^2 and -d not in Q_p^2."""
        out = []
        for p in self.bad_primes(x):
            if self.is_square_Qp(Fr(-self.d*self.m), p) and not self.is_square_Qp(Fr(-self.d), p):
                out.append(p)
        return out

    # ---------- norm equations in L/F ----------
    def rnf(self):
        if self._T is None:
            self._bnf = pari('bnfinit(y^2-%d,1)' % self.m)
            self._T = pari.rnfisnorminit(self._bnf, pari('x^2+%d' % self.d))
        return self._T
    def norm_LF_solve(self, t):
        """t in F (tuple). Return x in L with x conj(x) = t, or None."""
        tp = pari('%s+(%s)*y' % (str(t[0]), str(t[2])))
        res = pari.rnfisnorm(self.rnf(), tp)
        if int(res[1] != 1):
            return None
        sol = pari.lift(pari.lift(res[0]))            # polynomial in x (outer) with coeffs poly in y
        def coef(pol, var, k):
            c = pari.polcoef(pol, k, var)
            return c
        out = []
        for kx in range(2):
            cx = coef(sol, 'x', kx)
            for ky in range(2):
                c = coef(cx, 'y', ky)
                out.append(Fr(int(pari.numerator(c)), int(pari.denominator(c))))
        # out = [c(x^0 y^0), c(x^0 y^1), c(x^1 y^0), c(x^1 y^1)];  x = s, y = w
        z = (out[0], out[2], out[1], out[3])
        n = self.mul(z, self.conj(z))
        assert n == (t[0], Fr(0), t[2], Fr(0)), (n, t)
        return z

# ---------- rational quadratic forms: Witt index via qfsolve ----------
def _to_int_matrix(Q):
    den = 1
    for r in Q:
        for x in r:
            den = den * Fr(x).denominator // math.gcd(den, Fr(x).denominator)
    return [[int(Fr(x) * den) for x in r] for r in Q]

def witt_index(Q):
    """Witt index of the nondegenerate rational quadratic form with Gram matrix Q, computed by
    repeatedly finding an isotropic vector (PARI qfsolve, Simon's algorithm) and splitting off a
    hyperbolic plane.  Returns (index, reason) where reason is the qfsolve verdict on the anisotropic
    kernel (-1: real obstruction, p: p-adic obstruction, 'dim0')."""
    Qi = _to_int_matrix(Q)
    g = 0
    for r in Qi:
        for x in r:
            g = math.gcd(g, x)
    G = pari.matrix(len(Q), len(Q), [x // g for r in Qi for x in r])
    idx = 0
    while True:
        n = int(pari.matsize(G)[0])
        if n == 0:
            return idx, 'dim0'
        if n == 1:
            return idx, 'dim1'
        v = pari.qfsolve(G)
        if str(pari.type(v)) == 't_MAT':          # columns span a totally isotropic subspace
            v = v[0]
        if str(pari.type(v)) != 't_COL':
            return idx, int(v)
        assert int(pari.mattranspose(v) * G * v == 0)
        if n == 2:
            return idx + 1, 'dim0'
        v = v / pari.content(v)                    # primitive integral isotropic vector
        Gv = G * v
        # pick basis vector e_j with B(v, e_j) != 0 and |B(v, e_j)| minimal (keeps discriminants small)
        cand = [i for i in range(n) if int(Gv[i] != 0)]
        j = min(cand, key=lambda i: abs(int(Gv[i])))
        e = pari('vectorv(%d,k,k==%d)' % (n, j + 1))
        P = pari.matconcat([v, e])
        # orthogonal complement of span(v, e) (a hyperbolic plane) as a SATURATED sublattice of Z^n
        # (integral kernel), so that disc(complement) divides B(v,e)^2 * disc(G) and stays factorable
        Kb = pari.matkerint(pari.mattranspose(P) * G)
        Gn = pari.mattranspose(Kb) * G * Kb
        # clear denominators, keep integral
        den = pari.denominator(Gn)
        G = Gn * den
        G = G / pari.content(G)
        idx += 1
