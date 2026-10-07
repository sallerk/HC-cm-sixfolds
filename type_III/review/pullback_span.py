"""
pullback_span.py  (review of Theorem C / Proposition prop:typeIII-gen, fresh code)

Checks, with an actual definite quaternion algebra D = (a,c)_Q acting on V = D^3 (H_1 of A), an
imaginary quadratic K = Q(q) in D and its Weil classes w in W(A,K) = Lambda^6_K H^1(A):
  the span of the pull-backs f^*w, f = (d_1..d_k) in D^k = Hom(A^k, A) (x) Q, has dimension
  dim Sym^6(C^{2k}) = C(2k+5, 6)
(= the full Lambda^6 Y (x) Sym^6 U' piece; it is contained in it because each f^*w is
 SL(Y)-invariant, Y the multiplicity space of D (x) C on V (x) C).
Method: reduce mod a prime P with D (x) F_P split and -b a square mod P; then w_tau = wedge of
the 6 covectors spanning the tau-eigenspace of K on V^* (x) F_P.  For decomposable test 6-vectors
Xi (12k x 6 matrices) the pairing <f^* w_tau, Xi> = det(B_f Xi), B_f = [B L_{d_1} | ... | B L_{d_k}].
rank_P of the matrix (det(B_f Xi))_{f, Xi} is a lower bound for the dimension of the span over the
number field; it is compared with C(2k+5, 6).
Controls: (i) only f in D x {0} x ... (endomorphisms composed with projections): span 7k;
          (ii) a random 6-dim covector space instead of an eigenspace of K: rank exceeds C(2k+5,6).
"""
import sys, random, json, time
import numpy as np
import flint


def quat_mult_table(a, c):
    # basis 1, i, j, k ; i^2 = a, j^2 = c, ij = -ji = k
    # mult[x][y] = list of (coef, index)
    T = {}
    T[(0, 0)] = (1, 0); T[(0, 1)] = (1, 1); T[(0, 2)] = (1, 2); T[(0, 3)] = (1, 3)
    T[(1, 0)] = (1, 1); T[(1, 1)] = (a, 0); T[(1, 2)] = (1, 3); T[(1, 3)] = (a, 2)
    T[(2, 0)] = (1, 2); T[(2, 1)] = (-1, 3); T[(2, 2)] = (c, 0); T[(2, 3)] = (-c, 1)
    T[(3, 0)] = (1, 3); T[(3, 1)] = (-a, 2); T[(3, 2)] = (c, 1); T[(3, 3)] = (-a * c, 0)
    return T


def qmul(x, y, T):
    z = [0, 0, 0, 0]
    for s in range(4):
        if x[s] == 0:
            continue
        for t in range(4):
            if y[t] == 0:
                continue
            cf, u = T[(s, t)]
            z[u] += cf * x[s] * y[t]
    return z


def left_mult_matrix(d, T):
    """4x4 integer matrix of x -> d x (columns = images of basis vectors)"""
    M = np.zeros((4, 4), dtype=object)
    for t in range(4):
        e = [0, 0, 0, 0]; e[t] = 1
        img = qmul(d, e, T)
        for s in range(4):
            M[s, t] = img[s]
    return M


def L_on_V(d, T):
    M4 = left_mult_matrix(d, T)
    M = np.zeros((12, 12), dtype=object)
    for r in range(3):
        M[4 * r:4 * r + 4, 4 * r:4 * r + 4] = M4
    return M


def sqrt_mod(n, P):
    n %= P
    if n == 0:
        return 0
    if pow(n, (P - 1) // 2, P) != 1:
        return None
    # Tonelli-Shanks
    q, s = P - 1, 0
    while q % 2 == 0:
        q //= 2; s += 1
    z = 2
    while pow(z, (P - 1) // 2, P) != P - 1:
        z += 1
    m, c, t, r = s, pow(z, q, P), pow(n, q, P), pow(n, (q + 1) // 2, P)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % P; i += 1
        b = pow(c, 1 << (m - i - 1), P)
        m, c, t, r = i, b * b % P, t * b * b % P, r * b % P
    return r


def is_prime(n):
    if n < 2:
        return False
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def left_nullspace_mod(M, P):
    """basis (rows) of {x : x M = 0} over F_P, M given as list of lists"""
    A = flint.nmod_mat([[int(v) % P for v in row] for row in np.array(M).T.tolist()], P)  # M^T
    # nullspace of M^T  (column vectors) = left nullspace of M
    N, nullity = A.nullspace()
    rows = []
    for j in range(nullity):
        rows.append([int(N[i, j]) for i in range(N.nrows())])
    return rows


def batch_det_mod(Ms, P):
    """determinants mod P of a batch of n x n int64 matrices (entries in [0,P)), P < 2^31"""
    A = Ms.copy() % P
    N, n, _ = A.shape
    det = np.ones(N, dtype=np.int64)
    alive = np.ones(N, dtype=bool)
    for col in range(n):
        # pivot: first row >= col with nonzero entry
        sub = A[:, col:, col] != 0
        has = sub.any(axis=1)
        alive &= has
        piv = np.argmax(sub, axis=1) + col
        idxN = np.arange(N)
        swap = piv != col
        if swap.any():
            rows_col = A[idxN, col, :].copy()
            rows_piv = A[idxN, piv, :].copy()
            A[idxN, col, :] = rows_piv
            A[idxN, piv, :] = rows_col
            det[swap] = (-det[swap]) % P
        pv = A[:, col, col] % P
        pv_safe = np.where(pv == 0, 1, pv)
        det = (det * pv_safe) % P
        # inverse via Fermat
        inv = np.ones(N, dtype=np.int64); base = pv_safe.copy(); e = P - 2
        while e:
            if e & 1:
                inv = (inv * base) % P
            base = (base * base) % P
            e >>= 1
        for r in range(col + 1, n):
            f = (A[:, r, col] * inv) % P
            A[:, r, :] = (A[:, r, :] - (f[:, None] * A[:, col, :]) % P) % P
    det[~alive] = 0
    return det % P


def run_case(a, c, q, k, seed, P=None, n_extra=12, control=None, nf=None, nxi=None):
    T = quat_mult_table(a, c)
    b = -qmul(q, q, T)[0]
    assert qmul(q, q, T)[1:] == [0, 0, 0] and b > 0
    rnd = random.Random(seed)
    if P is None:
        P = 50_000_017          # P < 2^26: all int64 products below stay < 2^63
        while not (is_prime(P) and sqrt_mod(-b, P) is not None and (2 * a * c * b) % P != 0):
            P += 2
    s = sqrt_mod(-b, P)
    Lq = L_on_V(q, T)
    # covectors lambda (rows) with lambda L_q = s lambda  <=>  lambda (L_q - s I) = 0
    Mtau = [[(int(Lq[i, j]) - (s if i == j else 0)) % P for j in range(12)] for i in range(12)]
    Mbar = [[(int(Lq[i, j]) + (s if i == j else 0)) % P for j in range(12)] for i in range(12)]
    Btau = left_nullspace_mod(Mtau, P)
    Bbar = left_nullspace_mod(Mbar, P)
    assert len(Btau) == 6 and len(Bbar) == 6, (len(Btau), len(Bbar))
    if control == "random_subspace":
        Btau = [[rnd.randrange(P) for _ in range(12)] for _ in range(6)]
    Btau = np.array(Btau, dtype=np.int64); Bbar = np.array(Bbar, dtype=np.int64)
    target = 1
    for t in range(6):
        target = target * (2 * k + 5 - t) // (t + 1)      # C(2k+5, 6)
    if control == "projections":
        target_ctrl = 7 * k
    nf = nf or (target + n_extra)
    nxi = nxi or (target + n_extra)
    # random f in D^k with small integer coordinates (Hom(A^k,A) (x) Q)
    fs = []
    for _ in range(nf):
        if control == "projections":
            j = rnd.randrange(k)
            ds = [[0, 0, 0, 0] for _ in range(k)]
            ds[j] = [rnd.randint(-9, 9) for _ in range(4)]
        else:
            ds = [[rnd.randint(-9, 9) for _ in range(4)] for _ in range(k)]
        fs.append(ds)
    Ls = {}

    def Lmod(d):
        key = tuple(d)
        if key not in Ls:
            Ls[key] = np.array([[int(x) % P for x in row] for row in L_on_V(d, T)], dtype=np.int64)
        return Ls[key]

    def Bf(B, ds):
        return np.concatenate([(B @ Lmod(d)) % P for d in ds], axis=1)  # 6 x 12k
    Xis = np.array([[[rnd.randrange(P) for _ in range(6)] for _ in range(12 * k)] for _ in range(nxi)], dtype=np.int64)
    results = {}
    for label, B in (("tau", Btau), ("taubar", Bbar)):
        rows = []
        for ds in fs:
            Bfm = Bf(B, ds)                       # 6 x 12k
            # products B_f Xi for all Xi, mod P, done in two halves to avoid overflow: entries < 2^30
            prods = np.einsum('ij,njk->nik', Bfm, Xis % P) % P if 12 * k * P * P < 2 ** 63 else None
            if prods is None:
                prods = np.zeros((nxi, 6, 6), dtype=np.int64)
                for t in range(12 * k):
                    prods = (prods + (Bfm[:, t][None, :, None] * Xis[:, t, :][:, None, :]) % P) % P
            rows.append(batch_det_mod(prods, P).tolist())
        results[label] = rows
        if control == "random_subspace":
            break
    out = {}
    for label, rows in results.items():
        M = flint.nmod_mat(rows, P)
        out[label] = M.rank()
    if "taubar" in results:
        M = flint.nmod_mat(results["tau"] + results["taubar"], P)
        out["tau+taubar"] = M.rank()
    return dict(a=a, c=c, q=q, b=b, k=k, P=P, target=(target_ctrl if control == "projections" else target),
                dimSym6=target, nf=nf, nxi=nxi, control=control, ranks=out)


if __name__ == "__main__":
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    cases = [(-1, -1, [0, 1, 0, 0]),      # D=(-1,-1), K=Q(i)
             (-1, -1, [0, 1, 1, 1]),      # D=(-1,-1), K=Q(sqrt-3)
             (-2, -5, [0, 0, 1, 0]),      # disc 5, K=Q(sqrt-5)
             (-3, -7, [0, 1, 2, 0]),      # K = Q(sqrt(-3-28)) = Q(sqrt-31)
             (-11, -13, [0, 1, 1, 0])]    # K = Q(sqrt(-24))
    allres = []
    for (a, c, q) in cases:
        for k in range(1, kmax + 1):
            t0 = time.time()
            r = run_case(a, c, q, k, seed=17 * k + a)
            r["secs"] = round(time.time() - t0, 1)
            print(json.dumps(r)); sys.stdout.flush()
            allres.append(r)
    # controls
    for (a, c, q) in cases[:2]:
        for k in (2, 3):
            r = run_case(a, c, q, k, seed=99 + k, control="projections", nf=60, nxi=60)
            print("CONTROL projections only:", json.dumps(r)); sys.stdout.flush()
            allres.append(r)
            r = run_case(a, c, q, k, seed=199 + k, control="random_subspace")
            print("CONTROL random 6-dim subspace:", json.dumps(r)); sys.stdout.flush()
            allres.append(r)
    with open("pullback_span_out.json", "w") as f:
        json.dump(allres, f, indent=0)
