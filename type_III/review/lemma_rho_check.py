"""
lemma_rho_check.py  (review of Lemma lem:typeIII-rho; fresh code)

Exhibits a rational point of L(A) outside the identity component and checks that it acts by -1 on
the Weil classes W(A,K) for several K in D (so rho is onto {+-1}; ker rho is the identity component).
Setting: D = (a,c), T = <l, m l, n l> with l = i (pure), V = D^3 (left D-module, row vectors).
alpha = right multiplication by j on every coordinate: commutes with D, and
  T(x alpha, y alpha) = sum x_r j (t_r i) jbar ybar_r = -Nrd(j) T(x,y)   (j anticommutes with i),
so alpha alpha^dagger = nu = -Nrd(j) = c < 0.
With Milne's twisted action on H^6(A)(3) (alpha acts on H^1 = V^* by lambda -> lambda o alpha^{-1},
and on Q(1) by nu), alpha acts on the K-line W(A,K) = Lambda^6_K V^* (3) by nu^3 / det_K(alpha).
Homotheties t: nu = t^2, det_K = t^6 -> trivial (control).  An element of the identity component
(right multiplication by u in Q(i), u = x + y i, which preserves T up to Nrd(u)) acts trivially (control).
"""
from fractions import Fraction as Fr
from crit_check import Quat, KField, det_K, rank_Q, basis_V, vec, coords, from_coords, left, T_form


def right_mult(Qa, u, x):
    return [Qa.mul(xr, u) for xr in x]


def K_basis_any(Qa, q):
    chosen, rows = [], []
    for (r, beta) in basis_V():
        v = vec(r, beta)
        cand = rows + [coords(v), coords(left(Qa, q, v))]
        if rank_Q(cand) == len(cand):
            chosen.append(v); rows = cand
        if len(chosen) == 6:
            return chosen
    raise ValueError


def K_coords(Qa, q, Kb, w):
    """coordinates (in K = Q(q)) of w in the K-basis Kb: solve w = sum (u_r + v_r q) Kb[r] over Q"""
    cols = []
    for v in Kb:
        cols.append(coords(v)); cols.append(coords(left(Qa, q, v)))
    # solve 12x12 system  M sol = coords(w), M columns = cols
    n = 12
    M = [[Fr(cols[cidx][r]) for cidx in range(n)] + [Fr(coords(w)[r])] for r in range(n)]
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col] / M[col][col]
                M[r] = [x - f * y for x, y in zip(M[r], M[col])]
    sol = [M[r][n] / M[r][r] for r in range(n)]
    return [(sol[2 * t], sol[2 * t + 1]) for t in range(6)]


def detK_of_map(Qa, q, f):
    Kb = K_basis_any(Qa, q)
    Kf = KField(Qa.nrd(q))
    # matrix columns = K-coordinates of f(Kb[s])
    cols = [K_coords(Qa, q, Kb, f(v)) for v in Kb]
    G = [[cols[s][r] for s in range(6)] for r in range(6)]
    return det_K(G, Kf)


def nu_of(Qa, T, f):
    """nu with T(f x, f y) = nu T(x, y), checked on all basis pairs"""
    nu = None
    B = basis_V()
    for e1 in B:
        for e2 in B:
            x, y = vec(*e1), vec(*e2)
            t0 = T_form(Qa, T, x, y)
            t1 = T_form(Qa, T, f(x), f(y))
            if any(t0):
                r = next(Fr(t1[s]) / Fr(t0[s]) for s in range(4) if t0[s] != 0)
                if nu is None:
                    nu = r
                assert tuple(Fr(nu) * Fr(t) for t in t0) == tuple(Fr(t) for t in t1), "not a similitude"
            else:
                assert not any(t1)
    return nu


def run(a, c, mults, qs):
    Qa = Quat(a, c)
    i, j = (0, 1, 0, 0), (0, 0, 1, 0)
    z = (0, 0, 0, 0)
    T = [[Qa.scal(mults[0], i), z, z], [z, Qa.scal(mults[1], i), z], [z, z, Qa.scal(mults[2], i)]]
    out = []
    tests = [("alpha = right mult by j", lambda x: right_mult(Qa, j, x)),
             ("control: homothety 3", lambda x: [Qa.scal(3, xr) for xr in x]),
             ("control: right mult by 2+i (identity component)", lambda x: right_mult(Qa, (2, 1, 0, 0), x))]
    for name, f in tests:
        nu = nu_of(Qa, T, f)
        for q in qs:
            # f commutes with left multiplication by D
            for (r, beta) in basis_V()[:4]:
                v = vec(r, beta)
                for d in [(0, 1, 0, 0), (0, 0, 1, 0), (1, 2, 3, 4)]:
                    assert f(left(Qa, d, v)) == left(Qa, d, f(v))
            dK = detK_of_map(Qa, q, f)
            Kf = KField(Qa.nrd(q))
            nu3 = (Fr(nu) ** 3, Fr(0))
            rho = Kf.mul(nu3, Kf.inv(dK))
            out.append((name, q, nu, dK, rho))
            print(f"D=({a},{c}) T=<{mults[0]}i,{mults[1]}i,{mults[2]}i> {name}: K=Q(q), q={q}, Nrd(q)={Qa.nrd(q)}: "
                  f"nu={nu}, det_K={dK[0]}+{dK[1]}q, rho = nu^3/det_K = {rho[0]}+{rho[1]}q")
    return out


if __name__ == "__main__":
    res = []
    res += run(-1, -1, (1, 1, 1), [(0, 1, 0, 0), (0, 1, 1, 0), (0, 1, 1, 1), (0, 0, 1, 2)])
    res += run(-2, -5, (1, 2, 5), [(0, 1, 0, 0), (0, 0, 1, 0), (0, 1, 1, 1)])
    res += run(-3, -7, (1, 3, 1), [(0, 1, 0, 0), (0, 1, 2, 0), (0, 2, 1, 1)])
    bad = [r for r in res if r[0].startswith("alpha") and r[4] != (-1, 0)]
    badc = [r for r in res if r[0].startswith("control") and r[4] != (1, 0)]
    print(f"SUMMARY: non-identity-component element acts by -1 on W(A,K) in "
          f"{sum(1 for r in res if r[0].startswith('alpha')) - len(bad)}/{sum(1 for r in res if r[0].startswith('alpha'))} cases; "
          f"controls act trivially in {sum(1 for r in res if r[0].startswith('control')) - len(badc)}/"
          f"{sum(1 for r in res if r[0].startswith('control'))} cases")
