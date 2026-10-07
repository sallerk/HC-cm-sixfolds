"""INDEPENDENT VERIFIER (no PARI/cypari, no code shared with produce.py / produce_lib.gp).
Input: only the Weil polynomial P (integer coefficients), p, q and |G| (order of the Galois group of P,
taken from GAP GaloisType in verify_gap.py). Libraries: python-flint (FLINT), mpmath, fractions.
 V1  exact: monic, q-symmetry a_{12-i} = q^{6-i} a_i, q a power of p.
 V2  exact Weil property: real Weil polynomial h(y) (P(T) = T^6 h(T + q/T)) has 6 real roots in
     [-floor(2 sqrt q), floor(2 sqrt q)] (Sturm sequence over Q, own implementation) => all |roots| = sqrt q.
 V3  irreducible over Q (FLINT fmpz_poly.factor).
 V4  ordinary: p does not divide a_6.
 V5  geometric simplicity: Q_N = charpoly(x^N in F_l[x]/P) squarefree for a large prime l, N = lcm{n: phi(n) | |G|}
     (any root of unity in the splitting field L has order n with phi(n) | [L:Q] = |G|).
 V6  angle rank (numerical, mpmath): LLL (FLINT) on (theta_1/2pi, ..., theta_6/2pi, 1); count relations r;
     angle rank = 6 - r. Gram-Schmidt lower bound B for any further relation.
 V7  rigorous certificate of each found relation n: A'^2 = q^m with A' = prod_{n_k>0} pi_k^{n_k} prod_{n_k<0}
     conj(pi_k)^{|n_k|}, m = sum|n_k|; x = A'^2 - q^m is an algebraic integer of degree <= |G| with all
     conjugates <= 2 q^m in absolute value, so |x| < (2q^m)^{-(|G|-1)} forces x = 0 (precision chosen for this).
Usage: python verify.py <input.json> <output.json>; input: list of {name, weil_poly (T^12..T^0), p, q, G_order}."""
import sys, json, math
from fractions import Fraction
import mpmath
import flint

def is_prime_small(n):
    if n < 2: return False
    i = 2
    while i * i <= n:
        if n % i == 0: return False
        i += 1
    return True

def phi(n):
    r, m, d = n, n, 2
    while d * d <= m:
        if m % d == 0:
            while m % d == 0: m //= d
            r -= r // d
        d += 1
    if m > 1: r -= r // m
    return r

def roots_of_unity_bound(G):
    N = 1
    for n in range(1, 2 * G * G + 10):
        if G % phi(n) == 0:
            N = N * n // math.gcd(N, n)
    return N

def real_weil_poly(a, q):
    # s_k(y) = T^k + q^k T^-k: s_0 = 2, s_1 = y, s_{k+1} = y s_k - q s_{k-1}; polys as coefficient lists (low->high)
    s = [[2], [0, 1]]
    for k in range(1, 6):
        nxt = [0] + s[k]
        prev = s[k - 1] + [0] * (len(nxt) - len(s[k - 1]))
        s.append([nxt[i] - q * prev[i] for i in range(len(nxt))])
    h = [0] * 7
    h[0] += a[6]
    for k in range(1, 7):
        for i, c in enumerate(s[k]):
            h[i] += a[6 - k] * c
    return h  # low -> high

def poly_eval(c, x):  # c low->high, x Fraction
    r = Fraction(0)
    for co in reversed(c): r = r * x + co
    return r

def poly_rem(A, B):  # Fractions, low->high
    A = A[:]
    while len(A) >= len(B) and any(A):
        if A[-1] == 0: A.pop(); continue
        f = A[-1] / B[-1]; d = len(A) - len(B)
        for i in range(len(B)): A[d + i] -= f * B[i]
        A.pop()
    while A and A[-1] == 0: A.pop()
    return A

def sturm_count(c, lo, hi):
    P0 = [Fraction(x) for x in c]
    P1 = [Fraction(i * c[i]) for i in range(1, len(c))]
    seq = [P0, P1]
    while True:
        r = poly_rem(seq[-2], seq[-1])
        if not r: break
        seq.append([-x for x in r])
    def V(x):
        vals = [poly_eval(p, Fraction(x)) for p in seq]
        vals = [v for v in vals if v != 0]
        return sum(1 for i in range(len(vals) - 1) if (vals[i] > 0) != (vals[i + 1] > 0))
    assert poly_eval(c, Fraction(lo)) != 0 and poly_eval(c, Fraction(hi)) != 0
    return V(lo) - V(hi)

def geom_simple(a, N, ell):
    P = flint.nmod_poly(list(reversed(a)), ell)
    z = flint.nmod_poly([0, 1], ell).pow_mod(N, P)
    rows = []
    cur = flint.nmod_poly([1], ell)
    for i in range(12):
        v = (cur * z) % P
        co = [int(c) for c in v.coeffs()] + [0] * 12
        rows.append(co[:12])
        cur = (cur * flint.nmod_poly([0, 1], ell)) % P
    Mz = flint.nmod_mat([[rows[j][i] for j in range(12)] for i in range(12)], ell)
    Q = Mz.charpoly()
    g = Q.gcd(Q.derivative())
    return g.degree() == 0

def refine_roots(a, dps):
    mpmath.mp.dps = 60
    co = [mpmath.mpf(x) for x in a]
    rts = mpmath.polyroots(co, maxsteps=500, extraprec=200)
    mpmath.mp.dps = dps + 20
    co = [mpmath.mpf(x) for x in a]
    dco = [co[i] * (12 - i) for i in range(12)]
    out = []
    for r in rts:
        z = mpmath.mpc(r)
        for _ in range(int(math.log2(dps / 40.0 + 2)) + 8):
            z = z - mpmath.polyval(co, z) / mpmath.polyval(dco, z)
        out.append(z)
    return out

def gram_schmidt_norms(rows):
    mpmath.mp.dps = 60
    B = [[mpmath.mpf(int(x)) for x in r] for r in rows]
    Bs = []
    for b in B:
        v = b[:]
        for u in Bs:
            uu = sum(x * x for x in u)
            mu = sum(x * y for x, y in zip(b, u)) / uu
            v = [x - mu * y for x, y in zip(v, u)]
        Bs.append(v)
    return [mpmath.sqrt(sum(x * x for x in v)) for v in Bs]

def verify(ex):
    a = [int(x) for x in ex['weil_poly']]; p = int(ex['p']); q = int(ex['q']); G = int(ex['G_order'])
    out = dict(name=ex['name'])
    e = round(math.log(q, p)); out['q_is_power_of_p'] = (p ** e == q and is_prime_small(p))
    out['monic_symmetric'] = (len(a) == 13 and a[0] == 1 and all(a[12 - i] == q ** (6 - i) * a[i] for i in range(7)))
    h = real_weil_poly(a, q)
    # (i) h has 6 real roots: Sturm count on [-B, B], B = Cauchy bound
    B = 1 + max(abs(c) for c in h[:-1])
    # (ii) no root with y^2 > 4q: g(z) = h(y) h(-y), z = y^2 (degree-6 polynomial in z), count roots in (4q, B^2]
    hm = [c * (-1) ** i for i, c in enumerate(h)]
    prod = [0] * 13
    for i, c in enumerate(h):
        for j, d in enumerate(hm): prod[i + j] += c * d
    assert all(prod[k] == 0 for k in range(1, 13, 2))
    g = [prod[2 * k] for k in range(7)]
    out['weil_sturm'] = (sturm_count(h, -B, B) == 6 and sturm_count(g, 4 * q, B * B) == 0)
    fac = flint.fmpz_poly(list(reversed(a))).factor()
    out['irreducible'] = (len(fac[1]) == 1 and fac[1][0][1] == 1 and fac[1][0][0].degree() == 12)
    out['ordinary'] = (a[6] % p != 0)
    N = roots_of_unity_bound(G); out['N'] = N
    ells = [2**61 - 1, 2**62 - 57]  # primes (checked below)
    out['geom_simple'] = any(geom_simple(a, N, l) for l in ells if flint.fmpz(l).is_prime())
    # angle rank via LLL
    D = 420
    rts = refine_roots(a, D)
    mpmath.mp.dps = D
    up = sorted([r for r in rts if mpmath.im(r) > 0], key=lambda z: mpmath.arg(z))
    assert len(up) == 6
    out['max_abs_dev'] = float(max(abs(abs(r) / mpmath.sqrt(q) - 1) for r in rts))
    xs = [mpmath.arg(z) / (2 * mpmath.pi) for z in up] + [mpmath.mpf(1)]
    C = mpmath.mpf(10) ** 300
    rows = []
    for i in range(7):
        rows.append([1 if j == i else 0 for j in range(7)] + [int(mpmath.nint(C * xs[i]))])
    L = flint.fmpz_mat(rows).lll()
    Lr = [[int(L[i, j]) for j in range(8)] for i in range(7)]
    rels = [r[:7] for r in Lr if abs(r[7]) < 10**6 and max(abs(x) for x in r[:7]) < 10**6]
    gs = gram_schmidt_norms(Lr)
    nrel = len(rels)
    out['relations'] = rels
    out['angle_rank_numerical'] = 6 - nrel
    out['bound_no_further_relation'] = float(min(gs[nrel:])) if nrel < 7 else None
    # rigorous certificate for each relation
    certs = []
    for n in rels:
        nk = n[:6]
        m = sum(abs(x) for x in nk)
        need = (m + (G - 1) * (m + math.log(2, q)) + 2) * math.log10(q) + 30
        Dc = int(need) + 50
        rr = refine_roots(a, Dc)
        mpmath.mp.dps = Dc
        upc = sorted([r for r in rr if mpmath.im(r) > 0], key=lambda z: mpmath.arg(z))
        Ap = mpmath.mpc(1)
        for z, c in zip(upc, nk):
            Ap *= (z if c > 0 else mpmath.conj(z)) ** abs(c)
        x = Ap ** 2 - mpmath.mpf(q) ** m
        thr = mpmath.power(2 * mpmath.mpf(q) ** m, -(G - 1))
        certs.append(dict(relation=n, m=m, digits=Dc, log10_abs_x=float(mpmath.log10(abs(x))) if x != 0 else None,
                          log10_threshold=float(mpmath.log10(thr)), certified=bool(abs(x) < thr)))
    out['relation_certificates'] = certs
    out['ALL_OK'] = (out['q_is_power_of_p'] and out['monic_symmetric'] and out['weil_sturm'] and out['irreducible']
                     and out['ordinary'] and out['geom_simple'])
    return out

if __name__ == '__main__':
    exs = json.load(open(sys.argv[1]))
    res = []
    for ex in exs:
        r = verify(ex)
        res.append(r)
        print(json.dumps({k: v for k, v in r.items() if k not in ('relation_certificates',)}), flush=True)
        for c in r['relation_certificates']:
            print('   cert:', c, flush=True)
    json.dump(res, open(sys.argv[2], 'w'), indent=1)
