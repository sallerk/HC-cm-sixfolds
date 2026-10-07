"""
gen_check.py  (review of Theorem C / Proposition prop:typeIII-gen, fresh code)

Claim checked: for U' = C^{2k}, the SO(6)-invariants of Lambda^*(Y (x) U'), Y = C^6 with the
standard form sum y_i^2, are spanned by products of
   theta_{st} = sum_i (y_i (x) u_s) ^ (y_i (x) u_t)                    (s < t; degree 2, O(6)-invariant)
   Omega_nu   = sum_{(a_1..a_6) of content nu} (y_1(x)u_{a_1}) ^ ... ^ (y_6(x)u_{a_6})   (|nu| = 6; degree 6)
(the Omega_nu span Lambda^6 Y (x) Sym^6 U', i.e. the span of (y_1(x)u)^...^(y_6(x)u), u in U').
Over C the Hodge classes of A^k are exactly these invariants (Hg(A)_C = SO(beta_Y)); the theta_{st}
span the divisor classes and Lambda^6 Y (x) Sym^6 U' is the span of the pull-backs f^*w
(pullback_span.py checks that last identification with an actual quaternion algebra).

Method, per dominant GL(U')-weight mu in degree d:
  * upper bound  inv(mu) = dim of SO(6)-invariants in the mu-weight space  (Weyl character formula,
    inv_dims.py; no linear algebra);
  * lower bound  r(mu) = rank over F_p of the images of products of generators of weight mu under the
    projection v -> (sum_m v_m h(m)^j)_{j<R}, h a hash of the basis monomial m (Vandermonde-type
    projection).  rank_Fp(projection) <= rank_Fp(vectors) <= rank_Q(vectors) <= inv(mu), because the
    vectors are integral and SO(6)-invariant (checked below).  So r(mu) == inv(mu) PROVES that the
    products span the invariants of weight mu.  Monomials are processed in a pseudo-random order and
    the computation stops as soon as r == inv.  r < inv after all monomials is retried with another
    prime/hash before being reported as SHORT.
  Both spaces are GL(U')-stable, so equality on all dominant weights of degree d gives equality in
  degree d.

Controls: so(6)-invariance of the generators (and its failure for a planted non-invariant element);
mode theta_only (must fall SHORT from degree 6 on: divisor classes alone do not suffice);
mode bad (a planted non-invariant degree-2 generator must make the rank EXCEED inv somewhere).
"""
import itertools, sys, json, time, os, random
import numpy as np
from multiprocessing import Pool

from inv_dims import inv_dim, dominant_weights

NY = 6
PRIMES = [(11, 1000003), (23, 1000033), (37, 999983)]
VERIFY_PRODUCTS = os.environ.get("VERIFY_PRODUCTS", "1") == "1"


def idx(i, s):
    return s * NY + i


def wedge(a, b, p):
    """a, b: dict mask -> coef (mod p). returns a ^ b (a on the left)"""
    out = {}
    get = out.get
    for mb, cb in b.items():
        highs = []
        m = mb
        while m:
            low = m & -m
            highs.append(low.bit_length())      # shift amount j+1
            m ^= low
        for ma, ca in a.items():
            if ma & mb:
                continue
            c = 0
            for sh in highs:
                c += (ma >> sh).bit_count()
            v = ca * cb
            if c & 1:
                v = -v
            key = ma | mb
            out[key] = (get(key, 0) + v) % p
    return {k_: v for k_, v in out.items() if v}


_GEN_CACHE = {}


def make_generators(k, p, mode="normal"):
    key = (k, p, mode)
    if key in _GEN_CACHE:
        return _GEN_CACHE[key]
    k2 = 2 * k
    gens = []  # (name, weight tuple, dict)
    for s, t in itertools.combinations(range(k2), 2):
        d = {}
        for i in range(NY):
            a, b = idx(i, s), idx(i, t)          # a < b always
            d[(1 << a) | (1 << b)] = 1
        w = [0] * k2; w[s] += 1; w[t] += 1
        gens.append((f"th{s}{t}", tuple(w), d))
    if mode != "theta_only":
        for nu in itertools.product(range(7), repeat=k2):
            if sum(nu) != 6:
                continue
            content = []
            for s, c in enumerate(nu):
                content += [s] * c
            d = {}
            for seq in set(itertools.permutations(content)):
                elem = {0: 1}
                for i in range(NY):
                    elem = wedge(elem, {1 << idx(i, seq[i]): 1}, p)
                for m, c in elem.items():
                    d[m] = (d.get(m, 0) + c) % p
            d = {m: c for m, c in d.items() if c}
            gens.append((f"Om{''.join(map(str, nu))}", tuple(nu), d))
    if mode == "bad":
        w = [0] * k2; w[0] += 1; w[1] += 1
        gens.append(("BAD01", tuple(w), {(1 << idx(0, 0)) | (1 << idx(1, 1)): 1}))
    _GEN_CACHE[key] = gens
    return gens


def so6_action(elem, i, j, p):
    """derivation action of X_ij = E_ij - E_ji of so(6) on Lambda(Y (x) U')"""
    out = {}
    for m, c in elem.items():
        mm = m
        while mm:
            low = mm & -mm
            pos = low.bit_length() - 1
            mm ^= low
            s, a = divmod(pos, NY)
            if a == j:
                tgt, coef = idx(i, s), 1
            elif a == i:
                tgt, coef = idx(j, s), -1
            else:
                continue
            rest = m ^ (1 << pos)
            if rest & (1 << tgt):
                continue
            s1 = (m & ((1 << pos) - 1)).bit_count()
            s2 = (rest & ((1 << tgt) - 1)).bit_count()
            v = c * coef * (-1 if (s1 + s2) & 1 else 1)
            new = rest | (1 << tgt)
            out[new] = (out.get(new, 0) + v) % p
    return {m: c for m, c in out.items() if c}


def check_invariance(gens, p):
    bad = []
    for name, w, d in gens:
        for i, j in itertools.combinations(range(NY), 2):
            if so6_action(d, i, j, p):
                bad.append(name)
                break
    return bad


def monomials_of_weight(gens, mu):
    k2 = len(mu)
    res = []

    def rec(start, rem, cur):
        if not any(rem):
            res.append(tuple(cur))
            return
        for g in range(start, len(gens)):
            w = gens[g][1]
            ok = True
            for s in range(k2):
                if w[s] > rem[s]:
                    ok = False; break
            if not ok:
                continue
            rec(g, tuple(rem[s] - w[s] for s in range(k2)), cur + [g])
    rec(0, tuple(mu), [])
    return res


def hash_vals(masks, seed, p):
    with np.errstate(over='ignore'):
        z = masks + np.uint64((0x9E3779B97F4A7C15 + seed * 0xD1B54A32D192ED03) & 0xFFFFFFFFFFFFFFFF)
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        z = z ^ (z >> np.uint64(31))
    return (z % np.uint64(p - 1)).astype(np.int64) + 1     # nonzero values in [1, p-1]


def project(vec, R, seed, p):
    """(sum_m c_m h(m)^j mod p)_{j<R}"""
    if not vec:
        return [0] * R
    masks = np.fromiter(vec.keys(), dtype=np.uint64, count=len(vec))
    coefs = np.fromiter(vec.values(), dtype=np.int64, count=len(vec))
    h = hash_vals(masks, seed, p)
    out = []
    cur = coefs % p
    for _ in range(R):
        out.append(int(cur.sum() % p))      # entries < 2^20, so no int64 overflow below 2^43 terms
        cur = (cur * h) % p
    return out


class Echelon:
    """incremental row echelon basis over F_p"""
    def __init__(self, R, p):
        self.R, self.p = R, p
        self.rows = {}   # pivot col -> row (list), normalised pivot 1

    def add(self, v):
        p = self.p
        v = [x % p for x in v]
        for col in sorted(self.rows):
            if v[col]:
                f = v[col]
                r = self.rows[col]
                v = [(a - f * b) % p for a, b in zip(v, r)]
        for col, x in enumerate(v):
            if x:
                inv = pow(x, p - 2, p)
                v = [(a * inv) % p for a in v]
                # eliminate this col from existing rows (keep reduced form)
                for c2, r2 in self.rows.items():
                    if r2[col]:
                        f = r2[col]
                        self.rows[c2] = [(a - f * b) % p for a, b in zip(r2, v)]
                self.rows[col] = v
                return True
        return False

    def rank(self):
        return len(self.rows)


def rank_of_weight(args):
    k, mu, R_extra, mode = args
    t0 = time.time()
    target = inv_dim(mu)
    tries = []
    best = 0
    nmons = 0
    for seed, p in PRIMES:
        gens = make_generators(k, p, mode)
        mons = monomials_of_weight(gens, mu)
        nmons = len(mons)
        R = target + R_extra
        ech = Echelon(R, p)
        rnd = random.Random(1000 * seed + sum((i + 1) * x for i, x in enumerate(mu)))
        order = list(range(len(mons)))
        rnd.shuffle(order)
        # cache theta-only sub-products (up to 3 factors) to save work
        cache = {(): {0: 1}}

        def prod(key):
            if key in cache:
                return cache[key]
            v = wedge(prod(key[:-1]), gens[key[-1]][2], p)
            if len(key) <= 3:
                cache[key] = v
            return v
        used = 0
        n_inv_checked = 0
        for o in order:
            mon = mons[o]
            used += 1
            vec = prod(mon)
            added = ech.add(project(vec, R, seed, p))
            if added and VERIFY_PRODUCTS and mode == "normal":
                # every product that enlarges the span is re-checked to be so(6)-invariant
                for i, j in itertools.combinations(range(NY), 2):
                    if so6_action(vec, i, j, p):
                        raise RuntimeError(f"non-invariant product {mon} at weight {mu}")
                n_inv_checked += 1
            if mode == "normal" and ech.rank() >= target:
                break
            if mode == "bad" and ech.rank() > target:
                break
        tries.append([seed, p, used, ech.rank()] + ([n_inv_checked] if VERIFY_PRODUCTS and mode == "normal" else []))
        best = max(best, ech.rank())
        if mode != "normal" or best >= target:
            break
    return dict(mu=list(mu), inv=target, rank=best, n_monomials=nmons, tries=tries,
                secs=round(time.time() - t0, 2))


def main():
    k = int(sys.argv[1]); dmax = int(sys.argv[2])
    mode = sys.argv[3] if len(sys.argv) > 3 else "normal"
    dmin = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    nproc = min(7, int(os.environ.get("NPROC", "7")))
    p = PRIMES[0][1]
    gens = make_generators(k, p, "bad")
    bad = check_invariance(gens, p)
    print(f"[control] k={k}: {len(gens)} generators tested; failing so(6)-invariance: {bad}  (expected: only BAD01)")
    sys.stdout.flush()
    jobs = []
    for d in range(dmin, dmax + 1, 2):
        for mu in dominant_weights(2 * k, d):
            jobs.append((k, mu, 6, mode))
    # heaviest first for load balance
    jobs.sort(key=lambda j: -inv_dim(j[1]))
    out = []
    with Pool(nproc) as pool:
        for res in pool.imap_unordered(rank_of_weight, jobs):
            out.append(res)
            flag = "OK" if res["rank"] == res["inv"] else ("SHORT" if res["rank"] < res["inv"] else "EXCEEDS")
            print(f"k={k} mu={tuple(res['mu'])} deg={sum(res['mu'])} inv={res['inv']} rank={res['rank']} "
                  f"monomials={res['n_monomials']} tries={res['tries']} {flag} ({res['secs']}s)")
            sys.stdout.flush()
    by_deg = {}
    for r in out:
        d = sum(r["mu"])
        a = by_deg.setdefault(d, [0, 0, 0, 0])
        a[0] += r["inv"]; a[1] += r["rank"]; a[2] += (r["rank"] == r["inv"]); a[3] += 1
    print(f"summary k={k} mode={mode} (degree: sum inv over dominant weights, sum rank, #weights with rank==inv / #weights):")
    for d in sorted(by_deg):
        print(f"  {d}: {by_deg[d][0]} {by_deg[d][1]} {by_deg[d][2]}/{by_deg[d][3]}")
    out.sort(key=lambda r: (sum(r["mu"]), [-x for x in r["mu"]]))
    with open(f"gen_check_k{k}_d{dmin}-{dmax}_{mode}.json", "w") as f:
        json.dump(out, f, indent=0)


if __name__ == "__main__":
    main()
