# witness.py -- search for split E-compatible polarizations psi_a (own code; PARI only for bnfinit
# units and qfsolve as search tools).  delta(a) = (-1)^n N(a) disc(F) must lie in Nm(K^x).
import random
from fractions import Fraction
from cypari import pari
try:
    pari.allocatemem(600 * 10 ** 6)
except Exception:
    pass
from alg import Alg
from qsym import hilbert, factor_over, small_primes, squarefree_part, prime_factors
from inst import poly_y

PRIMES = small_primes(400)


def coeff_list_from_pari(e, m):
    e = pari.lift(e) if e.type() == 't_POLMOD' else e
    if e.type() in ('t_INT', 't_FRAC'):
        return [Fraction(str(e))] + [Fraction(0)] * (m - 1)
    return [Fraction(str(pari.polcoef(e, k, pari('y')))) for k in range(m)]


class Instance:
    def __init__(self, polys, d):
        self.polys = polys
        self.d = d
        self.Fa = Alg(polys)
        self.n = self.Fa.n
        self.N = 2 * self.n
        self.discF = self.Fa.disc()                  # power-basis discriminant (own code)
        self.discF_int = self.discF.numerator * self.discF.denominator
        self.base_primes = sorted(set([2] + prime_factors(d) + prime_factors(self.discF_int)))
        # t_p = ((-1)^n d_F, -d)_p
        self.t = {p: hilbert((-1) ** self.n * self.discF_int, -d, p) for p in self.base_primes}
        self._bnf = {}
        self._units_cache = {}
        self._sunit_cache = {}

    def bnf(self, i):
        if i not in self._bnf:
            self._bnf[i] = pari.bnfinit(pari(poly_y(self.polys[i])), 1)
        return self._bnf[i]

    # ---- PARI helpers (units only) ----
    def units(self, i):
        f = self.polys[i]
        m = len(f) - 1
        if m == 1:
            return []
        if i in self._units_cache:
            return self._units_cache[i]
        bnf = self.bnf(i)
        fu = pari('(b)->lift(b.fu)')(bnf)
        res = [coeff_list_from_pari(u, m) for u in fu]
        self._units_cache[i] = res
        return res

    # ---- delta and the norm test ----
    def delta(self, a):
        Na = self.Fa.norm(a)
        return (-1) ** self.n * Na * self.discF, Na

    def is_split(self, a):
        """True iff delta(a) in Nm(K^x) (own Hilbert symbols).  Returns (bool, bad primes, delta)."""
        dl, Na = self.delta(a)
        x = dl.numerator * dl.denominator
        if x < 0:
            return False, ('inf',), dl
        f, cof = factor_over(x, PRIMES)
        ps = set(f.keys()) | set(self.base_primes)
        if cof != 1:
            ps |= set(prime_factors(cof))
        bad = tuple(p for p in sorted(ps) if hilbert(x, -self.d, p) == -1)
        return len(bad) == 0, bad, dl

    def norm_certificate(self, dl):
        """qfsolve on X^2 + d Y^2 - delta_sf Z^2; returns verified (X,Y,Z,delta_sf) or None."""
        x = dl.numerator * dl.denominator
        f, cof = factor_over(x, PRIMES)
        sf = squarefree_part(x, PRIMES if cof == 1 else None)
        M = pari(f"matdiagonal([1,{self.d},{-sf}])")
        r = pari.qfsolve(M)
        if r.type() == 't_INT':
            return None
        X, Y, Z = (int(r[0]), int(r[1]), int(r[2]))
        assert X * X + self.d * Y * Y == sf * Z * Z and Z != 0
        return (X, Y, Z, sf)


def random_elt_with_signs(inst, eps, rng, B=4, maxtries=4000):
    """random small a in F^x with sign v(a) = eps(v) at every real place (factor by factor)."""
    Fa = inst.Fa
    a = []
    off = 0
    for F in Fa.factors:
        target = eps[off:off + F.m]
        off += F.m
        for _ in range(maxtries):
            c = [Fraction(rng.randint(-B, B)) for _ in range(F.m)]
            if all(x == 0 for x in c):
                continue
            try:
                if F.signs(c) == target:
                    break
            except ValueError:
                continue
        else:
            return None
        a.append(c)
    return a


# ---------------------------------------------------------------------------------------------
def gf2_solve(rows, target, ncols):
    """rows: list of int bitmasks (bit j = coefficient of unknown j); target: list of bits.
    Returns a solution bitmask or None."""
    R = [(r, t) for r, t in zip(rows, target)]
    piv = []
    sol_rows = []
    rank = 0
    pivcols = []
    R = list(R)
    for col in range(ncols):
        pr = None
        for i in range(rank, len(R)):
            if (R[i][0] >> col) & 1:
                pr = i
                break
        if pr is None:
            continue
        R[rank], R[pr] = R[pr], R[rank]
        for i in range(len(R)):
            if i != rank and (R[i][0] >> col) & 1:
                R[i] = (R[i][0] ^ R[rank][0], R[i][1] ^ R[rank][1])
        pivcols.append(col)
        rank += 1
    for i in range(rank, len(R)):
        if R[i][0] == 0 and R[i][1] == 1:
            return None
    x = 0
    for i, col in enumerate(pivcols):
        if R[i][1]:
            x |= (1 << col)
    return x


def sunits(inst, i, Pset):
    """S-unit generators of factor i for S = all primes above the rational primes in Pset (PARI)."""
    f = inst.polys[i]
    m = len(f) - 1
    if m == 1:
        return [[Fraction(p)] for p in sorted(Pset)]
    key = (i, tuple(sorted(Pset)))
    if key in inst._sunit_cache:
        return inst._sunit_cache[key]
    bnf = inst.bnf(i)
    S = pari('(b,L)->my(v=List()); for(k=1,#L, my(D=idealprimedec(b,L[k])); for(j=1,#D, listput(v,D[j]))); Vec(v)')(bnf, pari(str(sorted(Pset))))
    res = pari.bnfsunit(bnf, S)
    out = []
    for e in res[0]:
        if e.type() == 't_MAT':
            e = pari.nffactorback(bnf, e)
        e = pari.nfbasistoalg(bnf, e) if e.type() == 't_COL' else e
        out.append(coeff_list_from_pari(e, m))
    inst._sunit_cache[key] = out
    return out


def norm_primes(Nc):
    x = abs(Nc.numerator * Nc.denominator)
    f, cof = factor_over(x, PRIMES)
    ps = set(f.keys())
    if cof != 1:
        ps |= set(prime_factors(cof))
    return ps


def generators(inst, rng, Pset, n_small=8, B=3, smooth_bound=60, max_norm=10 ** 8):
    """list of (factor index, element coeffs, norm): -1, units, S-units (S above Pset), and a few
    random small elements with smooth norms (these keep witnesses small)."""
    gens = []
    Fa = inst.Fa
    for i, F in enumerate(Fa.factors):
        m = F.m
        gens.append((i, [Fraction(-1)] + [Fraction(0)] * (m - 1), Fraction((-1) ** m)))
        for u in inst.units(i):
            gens.append((i, u, F.norm(u)))
        for u in sunits(inst, i, Pset):
            gens.append((i, u, F.norm(u)))
        cnt = tries = 0
        seen = set()
        while m > 1 and cnt < n_small and tries < 5000:
            tries += 1
            c = tuple(Fraction(rng.randint(-B, B)) for _ in range(m))
            if all(x == 0 for x in c) or c in seen:
                continue
            seen.add(c)
            Nc = F.norm(list(c))
            if Nc == 0:
                continue
            x = abs(Nc.numerator * Nc.denominator)
            if x > max_norm:
                continue
            f, cof = factor_over(x, [p for p in PRIMES if p <= smooth_bound] + sorted(Pset))
            if cof != 1:
                continue
            gens.append((i, list(c), Nc))
            cnt += 1
    return gens


def linear_witness(inst, eps, rng, extra_rounds=4):
    """Solve for a = prod g^{e_g} with sign(a) = eps and (N(a),-d)_p = t_p for all p."""
    Fa = inst.Fa
    places = Fa.places()
    for rnd in range(extra_rounds):
        Pset = set(inst.base_primes) | set(p for p in PRIMES if p <= 10 + 30 * rnd)
        gens = generators(inst, rng, Pset)
        P = set(inst.base_primes)
        for (_, _, Nc) in gens:
            P |= norm_primes(Nc)
        P = sorted(P)
        ncols = len(gens)
        # rows: places then primes
        rows = []
        target = []
        sign_bits = []
        for (i, c, Nc) in gens:
            sg = Fa.factors[i].signs(c)
            full = []
            for (fi, ri) in places:
                full.append(1 if (fi == i and sg[ri] == -1) else 0)
            sign_bits.append(full)
        for r in range(len(places)):
            mask = 0
            for j in range(ncols):
                if sign_bits[j][r]:
                    mask |= (1 << j)
            rows.append(mask)
            target.append(1 if eps[r] == -1 else 0)
        for p in P:
            mask = 0
            for j, (i, c, Nc) in enumerate(gens):
                if hilbert(Nc, -inst.d, p) == -1:
                    mask |= (1 << j)
            rows.append(mask)
            tp = inst.t.get(p, 1)
            target.append(1 if tp == -1 else 0)
        sol = gf2_solve(rows, target, ncols)
        if sol is None:
            continue
        a = Fa.one()
        for j, (i, c, Nc) in enumerate(gens):
            if (sol >> j) & 1:
                a[i] = Fa.factors[i].mul(a[i], c)
        return a, rnd, ncols, len(rows)
    return None


# ---------------------------------------------------------------------------------------------
class SignPool:
    """Pools of small random elements per factor, indexed by sign vector (fast sampling of
    'random polarizations' with a prescribed sign pattern)."""
    def __init__(self, inst, rng, size=3000, B=3):
        import numpy as np
        self.inst = inst
        self.rng = rng
        self.pools = []
        for F in inst.Fa.factors:
            m = F.m
            roots = np.array([float(r) for r in F.roots])
            V = np.vander(roots, m, increasing=True)      # V[j,k] = r_j^k
            pool = {}
            seen = set()
            for _ in range(size):
                c = tuple(rng.randint(-B, B) for _ in range(m))
                if all(x == 0 for x in c) or c in seen:
                    continue
                seen.add(c)
                vals = V @ np.array(c, dtype=float)
                if np.min(np.abs(vals)) < 1e-6:
                    continue
                key = tuple(1 if v > 0 else -1 for v in vals)
                pool.setdefault(key, []).append(c)
            self.pools.append(pool)

    def sample(self, eps):
        """random a with sign pattern eps (verified with mpmath); products of two pool elements are
        used when a pattern is missing from the pool.  Returns None if impossible."""
        a = []
        off = 0
        for F, pool in zip(self.inst.Fa.factors, self.pools):
            tgt = tuple(eps[off:off + F.m])
            off += F.m
            if tgt in pool:
                c = [Fraction(x) for x in self.rng.choice(pool[tgt])]
            else:
                keys = list(pool.keys())
                self.rng.shuffle(keys)
                c = None
                for k1 in keys:
                    k2 = tuple(s * t for s, t in zip(k1, tgt))
                    if k2 in pool:
                        c1 = [Fraction(x) for x in self.rng.choice(pool[k1])]
                        c2 = [Fraction(x) for x in self.rng.choice(pool[k2])]
                        c = F.mul(c1, c2)
                        break
                if c is None:
                    return None
            if F.signs(c) != list(tgt):
                return None
            a.append(c)
        return a
