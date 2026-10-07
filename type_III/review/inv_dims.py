"""
inv_dims.py  (review of Theorem C, fresh code)

Dimension of SO(6)-invariants in the GL(U')-weight spaces of
    Lambda^*(Y (x) U'),   Y = C^6 (standard rep of SO(6)),  U' = C^{2k}.
The mu-weight space (mu = (mu_1..mu_{2k}), mu_s = number of tensor slots carrying u_s)
is  (x)_s Lambda^{mu_s}(Y)  as an SO(6)-module.  The multiplicity of the trivial
representation is computed with the Weyl character formula for D_3:
    mult_0(M) = sum_{w in W(D3)} det(w) * m_M(rho - w rho),   rho = (2,1,0).
No linear algebra is involved, so this is independent of gen_check.py.

Also computes (control) the multiplicity of the trivial rep of SL(6) restricted ... not needed.
"""
import itertools, functools, sys, json

N = 3  # rank of D_3
WEIGHTS_Y = []
for i in range(N):
    e = [0] * N; e[i] = 1; WEIGHTS_Y.append(tuple(e))
    e = [0] * N; e[i] = -1; WEIGHTS_Y.append(tuple(e))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


@functools.lru_cache(maxsize=None)
def char_wedge(a):
    """character of Lambda^a(C^6) as dict weight->mult"""
    ch = {}
    for S in itertools.combinations(range(6), a):
        w = (0,) * N
        for s in S:
            w = add(w, WEIGHTS_Y[s])
        ch[w] = ch.get(w, 0) + 1
    return ch


def mult(ch1, ch2):
    out = {}
    for w1, m1 in ch1.items():
        for w2, m2 in ch2.items():
            w = add(w1, w2)
            out[w] = out.get(w, 0) + m1 * m2
    return out


def perm_sign(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


# Weyl group of D_3: signed permutations with an even number of sign changes; det = sign(perm)
WEYL = []
for p in itertools.permutations(range(N)):
    for signs in itertools.product([1, -1], repeat=N):
        if signs.count(-1) % 2 == 0:
            WEYL.append((p, signs, perm_sign(p)))
assert len(WEYL) == 24
RHO = (2, 1, 0)


def w_act(w, v):
    p, signs, _ = w
    return tuple(signs[i] * v[p[i]] for i in range(N))


def mult_trivial(ch):
    tot = 0
    for w in WEYL:
        wr = w_act(w, RHO)
        key = tuple(RHO[i] - wr[i] for i in range(N))
        tot += w[2] * ch.get(key, 0)
    return tot


@functools.lru_cache(maxsize=None)
def inv_dim_sorted(mu_sorted):
    """dim of SO(6)-invariants in (x)_s Lambda^{mu_s}(C^6); mu given as sorted tuple"""
    ch = {(0,) * N: 1}
    for a in mu_sorted:
        if a == 0:
            continue
        ch = mult(ch, char_wedge(a))
    return mult_trivial(ch)


def inv_dim(mu):
    return inv_dim_sorted(tuple(sorted(mu, reverse=True)))


def dominant_weights(k2, d):
    """partitions of d into at most k2 parts, each <= 6, padded to length k2"""
    out = []

    def rec(prefix, remaining, maxpart):
        if len(prefix) == k2:
            if remaining == 0:
                out.append(tuple(prefix))
            return
        for a in range(min(maxpart, remaining), -1, -1):
            rec(prefix + [a], remaining - a, a)
    rec([], d, 6)
    return out


def sanity():
    # Lambda^2(C^6) = so(6): no invariants.  Lambda^3 (x) Lambda^3: one invariant (pairing). S^2 has one.
    assert inv_dim((1, 1)) == 1          # Y (x) Y -> 1 invariant (the form)
    assert inv_dim((2,)) == 0            # Lambda^2 Y
    assert inv_dim((6,)) == 1            # Lambda^6 Y = det, trivial on SO(6)
    assert inv_dim((3, 3)) == 2          # Lambda^3 (x) Lambda^3: pairing a.b and wedge a^b in Lambda^6 = det
    assert inv_dim((1, 1, 1, 1)) == 3    # (Y^{x4})^{SO6}: 3 pairings (no eps for 4 slots)
    # (Y^{x6})^{SO(6)} = 15 pairings + 1 eps = 16  (Brauer: (6-1)!! = 15)
    assert inv_dim((1,) * 6) == 16, inv_dim((1,) * 6)
    # control: O(6)-invariants in Y^{x6} would be 15; SO(6) gives 16 (det appears)
    return True


if __name__ == "__main__":
    sanity()
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    res = {}
    for d in range(0, dmax + 1, 2):
        rows = []
        for mu in dominant_weights(2 * k, d):
            rows.append((mu, inv_dim(mu)))
        res[d] = rows
        tot = sum(r[1] for r in rows)
        print(f"k={k} degree {d}: {len(rows)} dominant weights, sum of inv dims over dominant weights = {tot}")
        for mu, v in rows:
            print("   ", mu, v)
    with open(f"inv_dims_k{k}_d{dmax}.json", "w") as f:
        json.dump({str(d): [[list(mu), v] for mu, v in rows] for d, rows in res.items()}, f)
