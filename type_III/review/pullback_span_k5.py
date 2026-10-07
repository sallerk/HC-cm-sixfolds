"""pullback_span_k5.py: k = 5 (target C(15,6) = 5005) for two of the five test cases of pullback_span.py,
tau-eigenspace only.  Same method as pullback_span.run_case, but the pairing matrix is kept as a numpy
int64 array and its rank mod P is computed by numpy Gaussian elimination (the first attempt via python
lists and flint ran out of memory and ended without output).  The elimination routine is checked first
on random matrices of known rank (control)."""
import sys, json, time, random
import numpy as np
from pullback_span import quat_mult_table, qmul, L_on_V, sqrt_mod, is_prime, left_nullspace_mod, batch_det_mod


def rank_mod_p(M, P):
    A = np.array(M, dtype=np.int64) % P
    nr, nc = A.shape
    r = 0
    for c in range(nc):
        if r == nr:
            break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        piv = r + int(nz[0])
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), P - 2, P)
        A[r, c:] = (A[r, c:] * inv) % P
        f = A[r + 1:, c]
        nzr = np.nonzero(f)[0]
        if len(nzr):
            rows = r + 1 + nzr
            A[rows, c:] = (A[rows, c:] - (f[nzr, None] * A[r, c:][None, :]) % P) % P
        r += 1
    return r


def control_rank(P, seed=5):
    rnd = np.random.default_rng(seed)
    ok = 0
    for (n, m, rk) in [(40, 50, 17), (60, 45, 45), (80, 80, 79), (30, 30, 0)]:
        if rk == 0:
            M = np.zeros((n, m), dtype=np.int64)
        else:
            X = rnd.integers(0, P, size=(n, rk), dtype=np.int64)
            Y = rnd.integers(0, P, size=(rk, m), dtype=np.int64)
            M = np.zeros((n, m), dtype=np.int64)
            for t in range(rk):
                M = (M + (X[:, t:t + 1] * Y[t:t + 1, :]) % P) % P
        got = rank_mod_p(M, P)
        ok += (got == rk)
        print(f"[control] random {n}x{m} matrix of rank {rk}: rank_mod_p = {got}")
    return ok


def run_case_np(a, c, q, k, seed, n_extra=12):
    T = quat_mult_table(a, c)
    b = -qmul(q, q, T)[0]
    rnd = random.Random(seed)
    P = 50_000_017
    while not (is_prime(P) and sqrt_mod(-b, P) is not None and (2 * a * c * b) % P != 0):
        P += 2
    s = sqrt_mod(-b, P)
    Lq = L_on_V(q, T)
    Mtau = [[(int(Lq[i, j]) - (s if i == j else 0)) % P for j in range(12)] for i in range(12)]
    Btau = np.array(left_nullspace_mod(Mtau, P), dtype=np.int64)
    assert Btau.shape == (6, 12)
    target = 1
    for t in range(6):
        target = target * (2 * k + 5 - t) // (t + 1)
    nf = nxi = target + n_extra
    Ls = {}

    def Lmod(d):
        key = tuple(d)
        if key not in Ls:
            Ls[key] = np.array([[int(x) % P for x in row] for row in L_on_V(d, T)], dtype=np.int64)
        return Ls[key]

    Xis = np.array([[[rnd.randrange(P) for _ in range(6)] for _ in range(12 * k)] for _ in range(nxi)],
                   dtype=np.int64)
    assert 12 * k * P * P < 2 ** 63
    M = np.zeros((nf, nxi), dtype=np.int64)
    for i in range(nf):
        ds = [[rnd.randint(-9, 9) for _ in range(4)] for _ in range(k)]
        Bfm = np.concatenate([(Btau @ Lmod(d)) % P for d in ds], axis=1)
        prods = np.einsum('ij,njk->nik', Bfm, Xis) % P
        M[i, :] = batch_det_mod(prods, P)
    t0 = time.time()
    rk = rank_mod_p(M, P)
    return dict(a=a, c=c, q=q, b=b, k=k, P=P, target=target, nf=nf, nxi=nxi, rank_tau=rk,
                rank_secs=round(time.time() - t0, 1))


if __name__ == "__main__":
    P0 = 50_000_017
    nok = control_rank(P0)
    print(f"[control] rank routine correct on {nok}/4 test matrices"); sys.stdout.flush()
    res = []
    for (a, c, q) in [(-1, -1, [0, 1, 1, 1]), (-2, -5, [0, 0, 1, 0])]:
        t0 = time.time()
        r = run_case_np(a, c, q, 5, seed=5 * 17 + a)
        r["secs"] = round(time.time() - t0, 1)
        print(json.dumps(r)); sys.stdout.flush()
        res.append(r)
    json.dump(res, open("pullback_span_k5_out.json", "w"), indent=0)
