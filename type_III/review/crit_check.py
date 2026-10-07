"""
crit_check.py  (review of Theorem C(i), Lemma lem:typeIII-det, Prop. prop:typeIII-crit; fresh code)

For random definite D = (a,c)_Q, random non-degenerate D-valued skew-hermitian T on V = D^3 (not
necessarily diagonal) and random imaginary quadratic K = Q(q), q pure in D:
  * psi = Trd o T on V (12x12 over Q); the K-hermitian form phi(x,y) = psi(x,qy) + q psi(x,y)
    (Milne's normalisation, as in the note's section 2.1), built directly from psi;
  * det phi in Q (Gaussian elimination over K = Q(q), exact); signature of psi(x,qy) on V_R;
  * split  <=>  -det(phi) in Nm(K^x)  (Landherr; signature (3,3)):
      decided (A) locally: (-b, -det phi)_v = 1 for all places v (own Hilbert symbols), and
              (B) globally: an explicit rational solution of X^2 + b Y^2 = N Z^2, N = -det phi
                  (sympy's Legendre solver, solution verified), or a local obstruction otherwise;
  * Nrd(T) = sqrt(det_Q(x -> xT)) computed separately; Delta = -Nrd(T);
  * criterion of Prop. typeIII-crit(i):  (-b, Delta)_Q = D  (equal local invariants at all places);
  * Lemma typeIII-det: det phi == -b c Delta mod Nm(K^x), c = j^2 for some pure j anticommuting with q;
  * the note's identity phi = -2 conj(H_K), H_K = q pi_K(T(x,y))  (checked entrywise);
  * existence (Prop. typeIII-crit(ii)): local test "Delta not a square in Q_p for all p | disc D"
    versus a constructive search (pure y with y^2 in Delta Q^x2, x anticommuting with y, K = Q(x),
    split checked by the direct computation) and versus random K.
Controls (must FAIL to match): the R1 formula det = (-b)^3 Delta; the sign-flipped Delta' = +Nrd(T);
a wrong algebra D' in the criterion.  Hilbert symbols are cross-checked against the product formula
and against sympy's global solver.
"""
from fractions import Fraction as Fr
import random, sys, json, math, itertools
from sympy import factorint
from sympy.solvers.diophantine.diophantine import diop_ternary_quadratic_normal
from sympy import symbols

X_, Y_, Z_ = symbols("X Y Z", integer=True)
SYMPY_BAD = [0]

# ----------------------------------------------------------------------------- number theory


def sqfree_part(n):
    """squarefree integer with n / sqf a rational square (n nonzero rational)"""
    n = Fr(n)
    num, den = n.numerator, n.denominator
    m = num * den  # same square class
    s = -1 if m < 0 else 1
    m = abs(m)
    out = 1
    for p, e in factorint(m).items():
        if e % 2:
            out *= p
    return s * out


def vp(n, p):
    n = abs(n)
    v = 0
    while n % p == 0:
        n //= p; v += 1
    return v


def legendre(u, p):
    u %= p
    if u == 0:
        return 0
    return 1 if pow(u, (p - 1) // 2, p) == 1 else -1


def hilbert(a, b, p):
    """Hilbert symbol (a,b)_p for nonzero integers a,b (squarefree not required); p prime or 0 (= infinity)"""
    if p == 0:
        return -1 if (a < 0 and b < 0) else 1
    al, be = vp(a, p), vp(b, p)
    u, v = a // p ** al, b // p ** be
    if p != 2:
        eps = (p - 1) // 2
        r = (-1) ** (al * be * eps)
        r *= legendre(u, p) ** be * legendre(v, p) ** al
        return r
    eu, ev = ((u - 1) // 2) % 2, ((v - 1) // 2) % 2
    wu, wv = ((u * u - 1) // 8) % 2, ((v * v - 1) // 8) % 2
    e = (eu * ev + al * wv + be * wu) % 2
    return -1 if e else 1


def places(*nums):
    ps = {2}
    for n in nums:
        for p in factorint(abs(int(n))):
            ps.add(p)
    return sorted(ps) + [0]


def is_local_norm_everywhere(b, N):
    """N in Nm(Q(sqrt(-b))^x)  via Hasse norm theorem"""
    b, N = sqfree_part(b), sqfree_part(N)
    return all(hilbert(-b, N, p) == 1 for p in places(b, N))


def global_norm_solution(b, N):
    """explicit (X,Y,Z) != 0 with X^2 + b Y^2 = N Z^2, or None (sympy Legendre solver)"""
    b, N = sqfree_part(b), sqfree_part(N)
    if N == 1:
        return (1, 0, 1)
    try:
        sol = diop_ternary_quadratic_normal(X_ ** 2 + b * Y_ ** 2 - N * Z_ ** 2)
    except Exception:
        sol = None
    if not sol or sol == (None, None, None):
        return None
    x, y, z = [int(t) for t in sol]
    if (x, y, z) == (0, 0, 0):
        return None
    assert x * x + b * y * y == N * z * z
    return (x, y, z)


def ram_set(a, c):
    return frozenset(p for p in places(a, c) if hilbert(a, c, p) == -1)


def is_square_Qp(n, p):
    n = sqfree_part(n)
    if p == 0:
        return n > 0
    if vp(n, p) % 2:
        return False
    u = n // p ** vp(n, p)
    if p == 2:
        return u % 8 == 1
    return legendre(u, p) == 1


# ----------------------------------------------------------------------------- quaternions


class Quat:
    def __init__(self, a, c):
        self.a, self.c = a, c

    def mul(self, x, y):
        a, c = self.a, self.c
        x0, x1, x2, x3 = x; y0, y1, y2, y3 = y
        # i^2=a, j^2=c, ij=k=-ji, k^2=-ac, ik = a j, ki = -a j, jk = -c i, kj = c i
        return (x0 * y0 + a * x1 * y1 + c * x2 * y2 - a * c * x3 * y3,
                x0 * y1 + x1 * y0 - c * x2 * y3 + c * x3 * y2,
                x0 * y2 + x2 * y0 + a * x1 * y3 - a * x3 * y1,
                x0 * y3 + x3 * y0 + x1 * y2 - x2 * y1)

    @staticmethod
    def conj(x):
        return (x[0], -x[1], -x[2], -x[3])

    def nrd(self, x):
        a, c = self.a, self.c
        return x[0] ** 2 - a * x[1] ** 2 - c * x[2] ** 2 + a * c * x[3] ** 2

    @staticmethod
    def trd(x):
        return 2 * x[0]

    @staticmethod
    def add(x, y):
        return tuple(s + t for s, t in zip(x, y))

    @staticmethod
    def sub(x, y):
        return tuple(s - t for s, t in zip(x, y))

    @staticmethod
    def scal(r, x):
        return tuple(r * t for t in x)


def check_quat_table(Qa):
    i, j, k, one = (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (1, 0, 0, 0)
    a, c = Qa.a, Qa.c
    assert Qa.mul(i, i) == (a, 0, 0, 0) and Qa.mul(j, j) == (c, 0, 0, 0)
    assert Qa.mul(i, j) == k and Qa.mul(j, i) == (0, 0, 0, -1)
    assert Qa.mul(k, k) == (-a * c, 0, 0, 0)
    rnd = random.Random(1)
    for _ in range(50):
        x, y, z = [tuple(rnd.randint(-5, 5) for _ in range(4)) for _ in range(3)]
        assert Qa.mul(Qa.mul(x, y), z) == Qa.mul(x, Qa.mul(y, z))
        assert Qa.nrd(Qa.mul(x, y)) == Qa.nrd(x) * Qa.nrd(y)
        assert Qa.mul(x, Qa.conj(x)) == (Qa.nrd(x), 0, 0, 0)


# ----------------------------------------------------------------------------- linear algebra over Q / K


def det_Q(M):
    M = [[Fr(x) for x in row] for row in M]
    n = len(M); d = Fr(1)
    for col in range(n):
        piv = next((r for r in range(col, n) if M[r][col] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != col:
            M[col], M[piv] = M[piv], M[col]; d = -d
        d *= M[col][col]
        for r in range(col + 1, n):
            f = M[r][col] / M[col][col]
            if f:
                M[r] = [x - f * y for x, y in zip(M[r], M[col])]
    return d


def rank_Q(rows):
    M = [[Fr(x) for x in row] for row in rows]
    rk = 0; ncols = len(M[0]) if M else 0
    for col in range(ncols):
        piv = next((r for r in range(rk, len(M)) if M[r][col] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        for r in range(len(M)):
            if r != rk and M[r][col] != 0:
                f = M[r][col] / M[rk][col]
                M[r] = [x - f * y for x, y in zip(M[r], M[rk])]
        rk += 1
    return rk


def inertia(S):
    """(n_pos, n_neg, n_zero) of a rational symmetric matrix, exact (symmetric elimination)"""
    S = [[Fr(x) for x in row] for row in S]
    n = len(S); pos = neg = 0
    idx = list(range(n))
    while idx:
        # find a nonzero diagonal entry
        d = next((i for i in idx if S[i][i] != 0), None)
        if d is None:
            # find nonzero off-diagonal; replace e_i by e_i + e_j
            pair = next(((i, j) for i in idx for j in idx if i != j and S[i][j] != 0), None)
            if pair is None:
                break
            i, j = pair
            for t in range(n):
                S[i][t] += S[j][t]
            for t in range(n):
                S[t][i] += S[t][j]
            continue
        piv = S[d][d]
        if piv > 0:
            pos += 1
        else:
            neg += 1
        idx.remove(d)
        for i in idx:
            f = S[i][d] / piv
            if f:
                for t in idx:
                    S[i][t] -= f * S[d][t]
        for i in idx:
            S[i][d] = S[d][i] = Fr(0)
    return pos, neg, n - pos - neg


class KField:
    """K = Q(q), q^2 = -b; elements (u, v) = u + v q"""
    def __init__(self, b):
        self.b = Fr(b)

    def mul(self, x, y):
        return (x[0] * y[0] - self.b * x[1] * y[1], x[0] * y[1] + x[1] * y[0])

    def inv(self, x):
        n = x[0] ** 2 + self.b * x[1] ** 2
        return (x[0] / n, -x[1] / n)

    @staticmethod
    def sub(x, y):
        return (x[0] - y[0], x[1] - y[1])

    @staticmethod
    def conj(x):
        return (x[0], -x[1])


def det_K(G, Kf):
    G = [list(row) for row in G]
    n = len(G); d = (Fr(1), Fr(0))
    for col in range(n):
        piv = next((r for r in range(col, n) if G[r][col] != (0, 0)), None)
        if piv is None:
            return (Fr(0), Fr(0))
        if piv != col:
            G[col], G[piv] = G[piv], G[col]; d = (-d[0], -d[1])
        d = Kf.mul(d, G[col][col])
        ip = Kf.inv(G[col][col])
        for r in range(col + 1, n):
            f = Kf.mul(G[r][col], ip)
            if f != (0, 0):
                G[r] = [Kf.sub(x, Kf.mul(f, y)) for x, y in zip(G[r], G[col])]
    return d


# ----------------------------------------------------------------------------- the type III data


def basis_V():
    return [(r, beta) for r in range(3) for beta in range(4)]


def vec(r, beta):
    v = [(0, 0, 0, 0)] * 3
    e = [0, 0, 0, 0]; e[beta] = 1
    v[r] = tuple(e)
    return v


def T_form(Qa, T, x, y):
    tot = (0, 0, 0, 0)
    for r in range(3):
        for s in range(3):
            tot = Qa.add(tot, Qa.mul(Qa.mul(x[r], T[r][s]), Qa.conj(y[s])))
    return tot


def left(Qa, d, x):
    return [Qa.mul(d, xr) for xr in x]


def coords(x):
    return [x[r][beta] for r in range(3) for beta in range(4)]


def from_coords(cv):
    return [tuple(cv[4 * r:4 * r + 4]) for r in range(3)]


def random_pure(rnd, B=3):
    while True:
        p = (0, rnd.randint(-B, B), rnd.randint(-B, B), rnd.randint(-B, B))
        if any(p):
            return p


def random_T(rnd, Qa, kind):
    T = [[(0, 0, 0, 0)] * 3 for _ in range(3)]
    for r in range(3):
        T[r][r] = random_pure(rnd, 3)
    if kind != "diag":
        for r in range(3):
            for s in range(r + 1, 3):
                B = 2 if kind == "sparse" else 4
                t = tuple(rnd.randint(-B, B) for _ in range(4))
                if kind == "sparse" and rnd.random() < 0.5:
                    t = (0, 0, 0, 0)
                T[r][s] = t
                T[s][r] = Qa.scal(-1, Qa.conj(t))
    return T


def nrd_T(Qa, T):
    """Nrd of T in M_3(D) as sqrt(det_Q(x -> x T)) (left-D-linear map on D^3)"""
    rows = []
    for (r, beta) in basis_V():
        x = vec(r, beta)
        xT = [(0, 0, 0, 0)] * 3
        for s in range(3):
            acc = (0, 0, 0, 0)
            for t in range(3):
                acc = Qa.add(acc, Qa.mul(x[t], T[t][s]))
            xT[s] = acc
        rows.append(coords(xT))
    dq = det_Q(rows)
    if dq == 0:
        return Fr(0)
    num, den = dq.numerator, dq.denominator
    sn, sd = math.isqrt(num), math.isqrt(den)
    assert sn * sn == num and sd * sd == den, "det of right multiplication is not a square"
    return Fr(sn, sd)


def psi_matrix(Qa, T):
    B = basis_V()
    return [[Fr(Qa.trd(T_form(Qa, T, vec(*e1), vec(*e2)))) for e2 in B] for e1 in B]


def bil(Psi, x, y):
    cx, cy = coords(x), coords(y)
    return sum(cx[i] * Psi[i][j] * cy[j] for i in range(12) for j in range(12) if cx[i] and cy[j])


def K_basis(Qa, q):
    """6 vectors v_r from the Q-basis with {v_r, q v_r} a Q-basis of V"""
    chosen, rows = [], []
    for (r, beta) in basis_V():
        v = vec(r, beta)
        cand = rows + [coords(v), coords(left(Qa, q, v))]
        if rank_Q(cand) == len(cand):
            chosen.append(v); rows = cand
        if len(chosen) == 6:
            break
    assert len(chosen) == 6
    return chosen


def phi_gram(Qa, Psi, q, Kb):
    G = []
    for x in Kb:
        row = []
        for y in Kb:
            u = bil(Psi, x, left(Qa, q, y))       # psi(x, q y)
            v = bil(Psi, x, y)                    # psi(x, y)
            row.append((Fr(u), Fr(v)))            # u + v q
        G.append(row)
    return G


def pure_orthogonal(Qa, q):
    """a nonzero pure quaternion j anticommuting with q"""
    a, c = Qa.a, Qa.c
    # anticommute <=> orthogonal for the norm form on pure part: -a q1 x1 - c q2 x2 + a c q3 x3 = 0
    w = (-a * q[1], -c * q[2], a * c * q[3])
    for cand in [(w[1], -w[0], 0), (w[2], 0, -w[0]), (0, w[2], -w[1])]:
        if any(cand):
            j = (0,) + tuple(cand)
            assert Qa.mul(j, q) == Qa.scal(-1, Qa.mul(q, j))
            return j
    raise ValueError


def decompose_K_Kj(Qa, q, j, t):
    """t = alpha + beta j with alpha, beta in K = Q(q); returns (alpha, beta) as K-pairs"""
    qj = Qa.mul(q, j)
    basis = [(1, 0, 0, 0), q, j, qj]
    # solve t = x0*1 + x1*q + x2*j + x3*qj over Q
    M = [[Fr(basis[s][r]) for s in range(4)] + [Fr(t[r])] for r in range(4)]
    n = 4
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col] / M[col][col]
                M[r] = [x - f * y for x, y in zip(M[r], M[col])]
    sol = [M[r][n] / M[r][r] for r in range(n)]
    return (sol[0], sol[1]), (sol[2], sol[3])


def crit_holds(b, Delta, a, c):
    """(-b, Delta)_Q iso (a, c)_Q"""
    return all(hilbert(-sqfree_part(b), sqfree_part(Delta), p) == hilbert(a, c, p)
               for p in places(b, Delta, a, c))


def analyse(Qa, T, q, Psi=None, check_HK=False):
    a, c = Qa.a, Qa.c
    if Psi is None:
        Psi = psi_matrix(Qa, T)
    b = Fr(Qa.nrd(q))
    assert Qa.mul(q, q) == (-b, 0, 0, 0)
    Kb = K_basis(Qa, q)
    G = phi_gram(Qa, Psi, q, Kb)
    Kf = KField(b)
    # hermitian check
    for r in range(6):
        for s in range(6):
            assert G[s][r] == Kf.conj(G[r][s]), "phi not hermitian"
    dG = det_K(G, Kf)
    assert dG[1] == 0, "det phi not rational"
    detphi = dG[0]
    N = Nrd = nrd_T(Qa, T)
    Delta = -Nrd
    res = dict(b=int(sqfree_part(b)), detphi_sqf=int(sqfree_part(detphi)), Delta=int(sqfree_part(Delta)))
    res["detphi_neg"] = detphi < 0
    split_local = is_local_norm_everywhere(b, -detphi)
    res["split"] = split_local
    res["crit"] = crit_holds(b, Delta, a, c)
    # Lemma typeIII-det: det phi == -b c' Delta mod norms, c' = j^2
    j = pure_orthogonal(Qa, q)
    cj = Qa.mul(j, j)[0]
    res["lemma_det"] = is_local_norm_everywhere(b, detphi * (-b * cj * Delta))
    # controls
    res["ctrl_R1"] = is_local_norm_everywhere(b, -(-(b ** 3) * Delta))      # R1: det == (-b)^3 Delta
    res["ctrl_sign"] = crit_holds(b, Nrd, a, c)                              # Delta' = +Nrd(T)
    # sign-consistent variant of R1 (drops the factor c, i.e. replaces c by -1): det == b*Delta
    res["ctrl_noc"] = is_local_norm_everywhere(b, -(b * Delta))
    # perturbed discriminant: Delta * 3 (a wrong square class)
    res["ctrl_pert"] = crit_holds(b, 3 * Delta, a, c)
    if check_HK:
        # phi = -2 conj(H_K), H_K(x,y) = q * pi_K(T(x,y))
        for x in Kb[:3]:
            for y in Kb[:3]:
                t = T_form(Qa, T, x, y)
                alpha, beta = decompose_K_Kj(Qa, q, j, t)
                HK = Kf.mul((Fr(0), Fr(1)), alpha)
                u = bil(Psi, x, left(Qa, q, y)); v = bil(Psi, x, y)
                lhs = (Fr(u), Fr(v))
                rhs = (-2 * HK[0], 2 * HK[1])     # -2 conj(HK)
                assert lhs == rhs, ("phi != -2 conj H_K", lhs, rhs)
        res["phi_eq_m2conjHK"] = True
    return res, Psi, detphi, Delta


def embeds(Delta, a, c):
    R = ram_set(a, c)
    return all(not is_square_Qp(Delta, p) for p in R if p != 0)   # at infinity Delta < 0 always


def find_split_K_constructive_planes(Qa, Delta, rnd, tries=60):
    """pure y with y^2 in Delta Q^x2 found on random planes <u,v> of pure quaternions:
    solve N0(s u + t v) + Delta z^2 = 0 (sympy general ternary solver); returns (y, x) or None"""
    from sympy.solvers.diophantine.diophantine import diop_ternary_quadratic
    from sympy import symbols
    s_, t_, z_ = symbols("s t z", integer=True)
    a, c = Qa.a, Qa.c
    D0 = int(sqfree_part(Delta))
    for _ in range(tries):
        u, v = random_pure(rnd, 3), random_pure(rnd, 3)
        # Gram of N0 on (u, v):  N0(s u + t v) = A s^2 + 2 B s t + C t^2
        A = Qa.nrd(u); C = Qa.nrd(v)
        B2 = Qa.nrd(Qa.add(u, v)) - A - C          # = 2B
        if A * C * 4 - B2 * B2 == 0:
            continue
        sol = diop_ternary_quadratic(A * s_ ** 2 + B2 * s_ * t_ + C * t_ ** 2 + D0 * z_ ** 2)
        if sol is None or sol == (None, None, None) or not any(sol):
            continue
        ss, tt, zz = [int(w) for w in sol]
        y = Qa.add(Qa.scal(ss, u), Qa.scal(tt, v))
        if not any(y):
            continue
        n = Qa.nrd(y)
        if n != -D0 * zz * zz:
            # sympy's general ternary solver occasionally returns a wrong tuple (seen ~1/300 in a
            # separate test); such outputs are discarded, never used
            SYMPY_BAD[0] += 1
            continue
        x = pure_orthogonal(Qa, y)
        return y, x
    return None


def find_split_K_constructive(Qa, Delta, B=12):
    """pure y with y^2 in Delta Q^x2, then x anticommuting with y; returns x or None"""
    a, c = Qa.a, Qa.c
    target = sqfree_part(-Delta)
    best = None
    for y1 in range(-B, B + 1):
        for y2 in range(-B, B + 1):
            for y3 in range(-B, B + 1):
                y = (0, y1, y2, y3)
                if not any(y):
                    continue
                n = Qa.nrd(y)       # y^2 = -n ; need -n in Delta Q^x2  <=>  n in -Delta Q^x2
                if sqfree_part(n) == target:
                    x = pure_orthogonal(Qa, y)
                    return y, x
    return None


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 2026
    n_alg = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    n_T = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    n_K = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    rnd = random.Random(seed)
    # Hilbert symbol sanity: product formula, and agreement with global solvability (sympy)
    hs_ok = 0; hs_glob = 0
    for _ in range(300):
        u = rnd.choice([-1, 1]) * rnd.randint(1, 400); v = rnd.choice([-1, 1]) * rnd.randint(1, 400)
        prod = 1
        for p in places(u, v):
            prod *= hilbert(u, v, p)
        assert prod == 1, ("product formula fails", u, v)
        hs_ok += 1
        # u x^2 + v y^2 = z^2 solvable  <=>  all symbols 1
        loc = all(hilbert(u, v, p) == 1 for p in places(u, v))
        us, vs = sqfree_part(u), sqfree_part(v)
        try:
            sol = diop_ternary_quadratic_normal(us * X_ ** 2 + vs * Y_ ** 2 - Z_ ** 2)
            glob = sol is not None and sol != (None, None, None) and any(sol)
        except Exception:
            glob = False
        assert loc == glob, ("local/global mismatch", u, v, loc, glob)
        hs_glob += 1
    print(f"[hilbert] product formula ok {hs_ok}/300; local-vs-sympy-global agreement {hs_glob}/300")

    stats = dict(cases=0, split=0, mism_crit=0, lemma_fail=0, detphi_pos=0, sig_bad=0, sig_checked=0,
                 ctrl_R1_mism=0, ctrl_sign_mism=0, ctrl_wrongD_mism=0, ctrl_noc_mism=0, ctrl_pert_mism=0,
                 glob_checked=0, glob_fail=0, HK_checked=0)
    exist = dict(pairs=0, embeds=0, constructive_split=0, constructive_notfound=0, constructive_fail=0,
                 nonembed=0, nonembed_split_found=0, random_K_split_in_embed=0, embed_with_some_split=0,
                 split_but_not_embed=0)
    algs = []
    while len(algs) < n_alg:
        a, c = -rnd.randint(1, 40), -rnd.randint(1, 40)
        if sqfree_part(a) == 1 or sqfree_part(c) == 1:
            continue
        algs.append((a, c))
    algs = [(-1, -1), (-2, -5), (-1, -3), (-3, -5), (-1, -7)] + algs
    for (a, c) in algs:
        Qa = Quat(a, c)
        check_quat_table(Qa)
        R = ram_set(a, c)
        for tI in range(n_T):
            kind = ["diag", "sparse", "dense"][tI % 3]
            T = random_T(rnd, Qa, kind)
            if nrd_T(Qa, T) == 0:
                continue
            Psi = psi_matrix(Qa, T)
            # psi alternating
            assert all(Psi[i][j] == -Psi[j][i] for i in range(12) for j in range(12))
            Delta = None
            anysplit = False
            for kI in range(n_K):
                q = random_pure(rnd, 3)
                res, _, detphi, Delta = analyse(Qa, T, q, Psi, check_HK=(kI == 0 and tI < 2))
                stats["cases"] += 1
                stats["split"] += res["split"]
                stats["mism_crit"] += (res["split"] != res["crit"])
                stats["lemma_fail"] += (not res["lemma_det"])
                stats["detphi_pos"] += (not res["detphi_neg"])
                stats["ctrl_R1_mism"] += (res["split"] != res["ctrl_R1"])
                stats["ctrl_sign_mism"] += (res["split"] != res["ctrl_sign"])
                stats["ctrl_noc_mism"] += (res["split"] != res["ctrl_noc"])
                stats["ctrl_pert_mism"] += (res["split"] != res["ctrl_pert"])
                # wrong algebra in the criterion
                a2, c2 = -rnd.randint(1, 40), -rnd.randint(1, 40)
                stats["ctrl_wrongD_mism"] += (res["split"] != crit_holds(Fr(Qa.nrd(q)), Delta, a2, c2))
                stats["HK_checked"] += ("phi_eq_m2conjHK" in res)
                anysplit |= res["split"]
                if kI == 0:
                    # signature of the real form psi(x, q y) on V_R must be (6,6)
                    S = [[bil(Psi, from_coords([Fr(int(i == t)) for t in range(12)]),
                              left(Qa, q, from_coords([Fr(int(j_ == t)) for t in range(12)])))
                          for j_ in range(12)] for i in range(12)]
                    S = [[(S[i][j_] + S[j_][i]) / 2 for j_ in range(12)] for i in range(12)]
                    stats["sig_checked"] += 1
                    stats["sig_bad"] += (inertia(S) != (6, 6, 0))
                if kI < 2:
                    # global (explicit) norm check of the split decision
                    sol = global_norm_solution(Qa.nrd(q), -detphi)
                    stats["glob_checked"] += 1
                    stats["glob_fail"] += ((sol is not None) != res["split"])
            # existence
            exist["pairs"] += 1
            emb = embeds(Delta, a, c)
            more = anysplit
            for _ in range(6):
                q = random_pure(rnd, 4)
                r3, _, _, _ = analyse(Qa, T, q, Psi)
                more |= r3["split"]
            if more and not emb:
                exist["split_but_not_embed"] += 1          # would contradict Prop. typeIII-crit(ii) '=>'
            if emb:
                exist["embeds"] += 1
                exist["random_K_split_in_embed"] += more
                found = find_split_K_constructive_planes(Qa, Delta, rnd)
                if found is None:
                    exist["constructive_notfound"] += 1
                else:
                    y, x = found
                    r2, _, _, _ = analyse(Qa, T, x, Psi)
                    if r2["split"]:
                        exist["constructive_split"] += 1
                    else:
                        exist["constructive_fail"] += 1
                    more |= r2["split"]
                exist["embed_with_some_split"] += more
            else:
                exist["nonembed"] += 1
                exist["nonembed_split_found"] += more
        print(f"D=({a},{c}) ram={sorted(R)}: cumulative {stats}"); sys.stdout.flush()
    print("FINAL stats:", json.dumps(stats))
    print("FINAL existence:", json.dumps(exist))
    print("discarded invalid sympy ternary solutions:", SYMPY_BAD[0])


if __name__ == "__main__":
    main()
