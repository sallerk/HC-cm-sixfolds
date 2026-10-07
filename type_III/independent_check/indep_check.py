# indep_check.py -- independent re-implementation checking Lemma 6.1 and Proposition 6.2 of the note
# (printed below as "Lemma2", "Thm3" and, for the second statement of Prop. 6.2, "Cor4").
# Pure Python, exact rationals (fractions), own Hilbert symbols, own factorisation.
# No code or data shared with the main computation (path1.py / path2.py in the parent folder, which use
# PARI / M_2(Q(sqrt c)) / structure constants).
# Usage: python indep_check.py SEED NCASES NK OUT.json   (runs used: SEED = 101, 202, 303, 404, NCASES = 250,
# NK = 8; summary line on stdout saved as run_SEED.out.txt; OUT.json is written relative to the current directory)
#
# Route: D = (al,be)_Q, al,be < 0 (definite). V = D^3 (row vectors, D acting on the LEFT).
# T(x,y) = sum_{r,s} x_r M_rs conj(y_s), M skew-Hermitian (M_sr = -conj M_rs, M_rr pure).
# E = Trd o T (12x12 alternating, rational). K = Q(q), q pure random, b = Nrd(q), q^2 = -b.
# Hermitian form computed DIRECTLY from E and the 12x12 action of q:
#     H(x,y) = E(qx,y) + q E(x,y)  in K       (K-linear in x, Hermitian; Tr_{K/Q}(q^-1 H) = 2E)
# det_K of H on a RANDOM sparse K-basis of V (6 vectors of Q^12), Gaussian elimination over K.
# Nrd(T) = + sqrt(det_Q(x -> xM on Q^12))  (simple module of M_3(D): det = Nrd^2; Nrd >= 0 over H).
# Delta = -Nrd(T).
# Tests:  (i) det_K H in Q, < 0; signature of E(qx,y) is (6,6)
#         (ii) SPLIT := -det H in Nm K^x  (Hasse norm thm: (-b,-detH)_v = 1 all v)
#         (iii) THM3 := (-b,Delta)_v == (al,be)_v all v                     (Prop. 6.2, first statement)
#         (iv) LEMMA2 := detH * (-b c Delta) in Nm K^x, c = j'^2, j' pure anticommuting with q   (Lemma 6.1)
#         (v) COR4 := for each (D,T): [Q(sqrt Delta) embeds in D] vs [some K among many split]
#                    (Prop. 6.2, second statement)
import random, sys, json, math
from fractions import Fraction as Fr

# ---------------- integer factorisation (own) ----------------
def is_probable_prime(n):
    if n < 2: return False
    small = [2,3,5,7,11,13,17,19,23,29,31,37]
    for p in small:
        if n % p == 0: return n == p
    d, s = n-1, 0
    while d % 2 == 0: d//=2; s+=1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

def pollard_brent(n):
    if n % 2 == 0: return 2
    while True:
        y, c, m = random.randrange(1,n), random.randrange(1,n), 128
        g = r = q = 1
        while g == 1:
            x = y
            for _ in range(r): y = (y*y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r-k)):
                    y = (y*y + c) % n; q = q*abs(x-y) % n
                g = math.gcd(q, n); k += m
            r *= 2
        if g == n:
            g = 1
            while g == 1:
                ys = (ys*ys + c) % n; g = math.gcd(abs(x-ys), n)
        if g != n: return g

def factor(n, out=None):
    if out is None: out = {}
    n = abs(n)
    if n == 1: return out
    for p in range(2, 2000):
        while n % p == 0:
            out[p] = out.get(p,0)+1; n//=p
        if n == 1: return out
    stack = [n]
    while stack:
        m = stack.pop()
        if m == 1: continue
        if is_probable_prime(m): out[m] = out.get(m,0)+1; continue
        d = pollard_brent(m); stack += [d, m//d]
    return out

def primes_of(x):
    x = Fr(x)
    return set(factor(x.numerator).keys()) | set(factor(x.denominator).keys())

# ---------------- Hilbert symbol (own; Serre, Cours d'arithmetique III.1.2) ----------------
def sqfree_int(x):
    """integer in the square class of the nonzero rational x"""
    x = Fr(x); return x.numerator * x.denominator

def hilbert(a, b, p):
    a, b = sqfree_int(a), sqfree_int(b)
    if p == 0:  # real place
        return -1 if (a < 0 and b < 0) else 1
    al = 0
    while a % p == 0: a//=p; al+=1
    be = 0
    while b % p == 0: b//=p; be+=1
    if p == 2:
        eps = lambda u: ((u-1)//2) % 2
        om = lambda u: ((u*u-1)//8) % 2
        e = eps(a)*eps(b) + al*om(b) + be*om(a)
        return -1 if e % 2 else 1
    leg = lambda u: 1 if pow(u % p, (p-1)//2, p) == 1 else -1
    s = (-1) ** ((al*be*((p-1)//2)) % 2)
    if be % 2: s *= leg(a)
    if al % 2: s *= leg(b)
    return s

def places(*xs):
    S = {2}
    for x in xs: S |= primes_of(x)
    return [0] + sorted(S)

def is_norm_from(minus_b, x):
    """x in Nm(Q(sqrt(minus_b))^x) ?  (Hasse norm theorem, cyclic)"""
    return all(hilbert(minus_b, x, v) == 1 for v in places(minus_b, x))

# ---------------- quaternion algebra (al,be) ----------------
def qmul(x, y, al, be):
    a1,b1,c1,d1 = x; a2,b2,c2,d2 = y
    return (a1*a2 + al*b1*b2 + be*c1*c2 - al*be*d1*d2,
            a1*b2 + b1*a2 - be*c1*d2 + be*d1*c2,
            a1*c2 + c1*a2 + al*b1*d2 - al*d1*b2,
            a1*d2 + d1*a2 + b1*c2 - c1*b2)
def qconj(x): return (x[0], -x[1], -x[2], -x[3])
def qadd(x, y): return tuple(u+v for u,v in zip(x,y))
def qnrd(x, al, be): return x[0]**2 - al*x[1]**2 - be*x[2]**2 + al*be*x[3]**2
Z4 = (0,0,0,0)

# ---------------- exact linear algebra ----------------
def det_Q(A):
    A = [[Fr(v) for v in row] for row in A]; n = len(A); d = Fr(1)
    for i in range(n):
        piv = next((r for r in range(i,n) if A[r][i] != 0), None)
        if piv is None: return Fr(0)
        if piv != i: A[i],A[piv] = A[piv],A[i]; d = -d
        d *= A[i][i]
        for r in range(i+1,n):
            f = A[r][i]/A[i][i]
            if f: A[r] = [A[r][k]-f*A[i][k] for k in range(n)]
    return d

def rank_Q(A):
    A = [[Fr(v) for v in row] for row in A]; m = len(A); n = len(A[0]); rk = 0
    for c in range(n):
        piv = next((r for r in range(rk,m) if A[r][c] != 0), None)
        if piv is None: continue
        A[rk],A[piv] = A[piv],A[rk]
        for r in range(m):
            if r != rk and A[r][c] != 0:
                f = A[r][c]/A[rk][c]; A[r] = [A[r][k]-f*A[rk][k] for k in range(n)]
        rk += 1
    return rk

def inertia_sym(S):
    """(n_plus, n_minus) of a rational symmetric matrix, exact (random congruence + LDL)."""
    n = len(S)
    for _ in range(20):
        P = [[random.randint(-3,3) for _ in range(n)] for _ in range(n)]
        if det_Q(P) == 0: continue
        PtS = [[sum(P[k][i]*S[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        A = [[sum(PtS[i][k]*P[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        A = [[Fr(v) for v in row] for row in A]
        pos = neg = 0; ok = True
        for i in range(n):
            if A[i][i] == 0: ok = False; break
            if A[i][i] > 0: pos += 1
            else: neg += 1
            for r in range(i+1,n):
                f = A[r][i]/A[i][i]
                if f: A[r] = [A[r][k]-f*A[i][k] for k in range(n)]
        if ok: return pos, neg
    return None

# K = Q(q), q^2 = -b ; element (u,v) = u + v q
def kmul(x, y, b): return (x[0]*y[0] - b*x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def kinv(x, b):
    n = x[0]**2 + b*x[1]**2; return (x[0]/n, -x[1]/n)
def ksub(x, y): return (x[0]-y[0], x[1]-y[1])
def det_K(A, b):
    A = [[(Fr(u),Fr(v)) for (u,v) in row] for row in A]; n = len(A); d = (Fr(1),Fr(0))
    for i in range(n):
        piv = next((r for r in range(i,n) if A[r][i] != (0,0)), None)
        if piv is None: return (Fr(0),Fr(0))
        if piv != i: A[i],A[piv] = A[piv],A[i]; d = (-d[0],-d[1])
        d = kmul(d, A[i][i], b); inv = kinv(A[i][i], b)
        for r in range(i+1,n):
            f = kmul(A[r][i], inv, b)
            if f != (0,0): A[r] = [ksub(A[r][k], kmul(f, A[i][k], b)) for k in range(n)]
    return d

# ---------------- the abelian-variety data ----------------
def vec_to_quats(x): return [tuple(x[4*r:4*r+4]) for r in range(3)]
def quats_to_vec(qs): return [c for qq in qs for c in qq]

class Setup:
    def __init__(self, al, be, M):
        self.al, self.be, self.M = al, be, M
        self.basis = [[1 if i == k else 0 for i in range(12)] for k in range(12)]
        self.G = [[self.E(ei, ej) for ej in self.basis] for ei in self.basis]
    def T(self, x, y):
        X, Y = vec_to_quats(x), vec_to_quats(y); s = Z4
        for r in range(3):
            for t in range(3):
                s = qadd(s, qmul(qmul(X[r], self.M[r][t], self.al, self.be), qconj(Y[t]), self.al, self.be))
        return s
    def E(self, x, y): return 2*self.T(x, y)[0]
    def Eg(self, x, y):  # via Gram matrix
        return sum(x[i]*self.G[i][j]*y[j] for i in range(12) for j in range(12) if x[i] and y[j])
    def leftmul(self, q, x):
        return quats_to_vec([qmul(q, xx, self.al, self.be) for xx in vec_to_quats(x)])
    def nrdT(self):
        # matrix of x -> xM on Q^12 ; det = Nrd(M)^2
        R = []
        for e in self.basis:
            X = vec_to_quats(e); img = []
            for t in range(3):
                s = Z4
                for r in range(3): s = qadd(s, qmul(X[r], self.M[r][t], self.al, self.be))
                img.append(s)
            R.append(quats_to_vec(img))
        d = det_Q(R)
        assert d > 0, d
        num, den = d.numerator, d.denominator
        rn, rd = math.isqrt(num), math.isqrt(den)
        assert rn*rn == num and rd*rd == den, "det not a square"
        return Fr(rn, rd)

def rand_pure(lo=-2, hi=2):
    while True:
        x = (0, random.randint(lo,hi), random.randint(lo,hi), random.randint(lo,hi))
        if x != Z4: return x
def rand_quat(lo=-2, hi=2, pzero=0.3):
    if random.random() < pzero: return Z4
    return tuple(random.randint(lo,hi) for _ in range(4))

def rand_M(kind):
    M = [[Z4]*3 for _ in range(3)]
    for r in range(3): M[r][r] = rand_pure()
    if kind != 'diag':
        for r in range(3):
            for s in range(r+1,3):
                m = rand_quat(pzero=0.2 if kind == 'dense' else 0.6)
                M[r][s] = m; M[s][r] = tuple(-c for c in qconj(m))
    return M

def anticommuting_pure(q, al, be):
    q1,q2,q3 = q[1],q[2],q[3]
    # B(x,y) = -al x1y1 - be x2y2 + al be x3y3 ; want B(q,j)=0
    cands = [(0, be*q2, -al*q1, 0), (0, al*be*q3, 0, al*q1), (0, 0, al*q3, be*q2)]
    for j in cands:
        if j != Z4: return j

def hermitian_det(S, q, b, rng_basis=True):
    # random sparse K-basis
    for _ in range(200):
        vs = []
        for a in range(6):
            v = [0]*12
            for _k in range(random.randint(1,3)): v[random.randrange(12)] = random.choice([-2,-1,1,1,2])
            vs.append(v)
        qvs = [S.leftmul(q, v) for v in vs]
        if rank_Q(vs + qvs) == 12: break
    else:
        raise RuntimeError("no K-basis")
    Hm = [[(S.Eg(qvs[a], vs[c]), S.Eg(vs[a], vs[c])) for c in range(6)] for a in range(6)]
    # Hermitian check: H(c,a) = conj H(a,c)
    for a in range(6):
        for c in range(6):
            assert Hm[c][a] == (Hm[a][c][0], -Hm[a][c][1])
    return det_K(Hm, b), vs, qvs

def embeds(Delta, al, be):
    """Q(sqrt Delta) embeds in D=(al,be)  <=> Delta non-square in Q_p for every p ramified in D."""
    ram = [p for p in places(al, be) if p != 0 and hilbert(al, be, p) == -1]
    def is_sq_p(x, p):
        x = sqfree_int(x); v = 0
        while x % p == 0: x//=p; v+=1
        if v % 2: return False
        if p == 2: return x % 8 == 1
        return pow(x % p, (p-1)//2, p) == 1
    return all(not is_sq_p(Delta, p) for p in ram), ram

ALGS = [(-1,-1),(-1,-3),(-2,-5),(-1,-7),(-3,-5),(-2,-13),(-1,-11),(-5,-7),(-3,-11),(-7,-19),(-2,-3),(-6,-10)]

def one_case(al, be, M, qs, do_sig):
    S = Setup(al, be, M)
    # sanity: E alternating, nondegenerate
    for i in range(12):
        assert S.G[i][i] == 0
        for j in range(12): assert S.G[i][j] == -S.G[j][i]
    if det_Q(S.G) == 0: return None
    N = S.nrdT(); Delta = -N
    rows = []
    for q in qs:
        b = qnrd(q, al, be)
        qq = qmul(q, q, al, be); assert qq == (-b,0,0,0)
        # compatibility E(qx,y) = E(x, qbar y) on basis
        for (i,j) in [(0,5),(3,7),(1,11),(6,10)]:
            ei, ej = S.basis[i], S.basis[j]
            assert S.Eg(S.leftmul(q,ei), ej) == S.Eg(ei, S.leftmul(qconj(q),ej))
        dH, vs, qvs = hermitian_det(S, q, b)
        assert dH[1] == 0, "det not rational"
        dH = dH[0]
        sig = None
        if do_sig:
            Sym = [[S.Eg(S.leftmul(q, S.basis[i]), S.basis[j]) for j in range(12)] for i in range(12)]
            for i in range(12):
                for j in range(12): assert Sym[i][j] == Sym[j][i]
            sig = inertia_sym(Sym)
        split = is_norm_from(-b, -dH)
        thm3 = all(hilbert(-b, Delta, v) == hilbert(al, be, v) for v in places(-b, Delta, al, be))
        j = anticommuting_pure(q, al, be)
        assert qmul(j, q, al, be) == tuple(-c for c in qmul(q, j, al, be))
        c = -qnrd(j, al, be)
        lemma2 = is_norm_from(-b, dH * (-b*c*Delta))
        rows.append(dict(q=q, b=b, c=c, detH=str(dH), split=split, thm3=thm3, lemma2=lemma2, sig=sig))
    return dict(al=al, be=be, M=[[list(m) for m in row] for row in M], Nrd=str(N), Delta=str(Delta), rows=rows)

def main(seed, ncases, nK, out):
    random.seed(seed)
    res = []
    for n in range(ncases):
        if n % 3 == 0:
            al, be = random.choice(ALGS)
        else:
            al, be = -random.randint(1,30), -random.randint(1,30)
        kind = random.choice(['diag','sparse','dense'])
        M = rand_M(kind)
        qs = []
        while len(qs) < nK:
            q = rand_pure(-3,3)
            g = math.gcd(math.gcd(q[1],q[2]),q[3])
            if g == 1 and q not in qs: qs.append(q)
        r = one_case(al, be, M, qs, do_sig=(n % 5 == 0))
        if r is None: continue
        r['kind'] = kind
        r['embeds'], r['ram'] = embeds(Fr(r['Delta']), al, be)
        res.append(r)
    json.dump(res, open(out,'w'))
    tot = sum(len(r['rows']) for r in res)
    nsplit = sum(x['split'] for r in res for x in r['rows'])
    bad3 = sum(x['split'] != x['thm3'] for r in res for x in r['rows'])
    bad2 = sum(not x['lemma2'] for r in res for x in r['rows'])
    sigs = [tuple(x['sig']) for r in res for x in r['rows'] if x['sig'] is not None]
    negdet = sum(Fr(x['detH']) < 0 for r in res for x in r['rows'])
    cor_bad = sum((any(x['split'] for x in r['rows']) and not r['embeds']) for r in res)
    cor_emb_nosplit = sum((r['embeds'] and not any(x['split'] for x in r['rows'])) for r in res)
    print(f"seed {seed}: (D,T) {len(res)}, (D,T,K) {tot}, split {nsplit}, split!=Thm3 {bad3}, Lemma2 failures {bad2}, "
          f"detH<0 {negdet}/{tot}, signatures {set(sigs)} ({len(sigs)} computed), "
          f"split-but-not-embedded {cor_bad}, embedded-but-no-sampled-K-split {cor_emb_nosplit}, "
          f"embedded {sum(r['embeds'] for r in res)}")

if __name__ == '__main__':
    seed, ncases, nK, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    main(seed, ncases, nK, out)
