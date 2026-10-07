# pullback_span_modp.py -- span of the pull-backs f^*(det_W (x) u^6) along random f = (d_1,...,d_k) in D^k
# (u an eigenvector of q, K = Q(q) fixed), i.e. the vectors (d_1 u, ..., d_k u)^6 in Sym^6(U^k), U = C^2.
# Computed mod a prime P in which a and -b are squares (s = sqrt(a), t = sqrt(-b) mod P), instead of over
# the number field L = Q(sqrt a, sqrt -b) (an exact computation over L was too slow). Rank mod P <= rank over L,
# so FULL rank mod P certifies full rank over L (one-sided certificate, which is the direction needed).
# Writes pullback_span_modp_out.txt next to this script.
import itertools, random, os
HERE = os.path.dirname(os.path.abspath(__file__))
random.seed(7)

def sqrt_mod(x, P):
    x %= P
    for r in range(P):
        if r*r % P == x: return r
    return None

def rank_mod(M, P):
    M = [row[:] for row in M]; r = 0; ncols = len(M[0])
    for col in range(ncols):
        piv = next((i for i in range(r, len(M)) if M[i][col] % P), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][col], P-2, P)
        M[r] = [x*inv % P for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][col] % P:
                f = M[i][col]
                M[i] = [(x - f*y) % P for x, y in zip(M[i], M[r])]
        r += 1
    return r

out = []
for (a, c, qv) in [(-1, -1, (1, 1, 0)), (-2, -5, (1, 0, 1)), (-3, -10, (0, 1, 1)), (-7, -15, (1, 2, 1))]:
    b = -(a*qv[0]**2 + c*qv[1]**2 - a*c*qv[2]**2)
    P = next(p for p in range(10007, 200000) if all(p % d for d in range(2, int(p**0.5)+1))
             and sqrt_mod(a, p) is not None and sqrt_mod(-b, p) is not None)
    s = sqrt_mod(a, P); t = sqrt_mod(-b, P)
    def m2(x):
        x0, x1, x2, x3 = x
        I = [[s, 0], [0, -s]]; J = [[0, c], [1, 0]]
        K = [[sum(I[i][k]*J[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
        return [[(x0*(i == j) + x1*I[i][j] + x2*J[i][j] + x3*K[i][j]) % P for j in range(2)] for i in range(2)]
    q = m2((0,) + qv)
    q2 = [[sum(q[i][k]*q[k][j] for k in range(2)) % P for j in range(2)] for i in range(2)]
    assert q2 == [[(-b) % P, 0], [0, (-b) % P]]
    # eigenvector of q for eigenvalue t: (q - t) u = 0
    A = [[(q[0][0]-t) % P, q[0][1]], [q[1][0], (q[1][1]-t) % P]]
    u = [A[0][1] % P, (-A[0][0]) % P] if (A[0][0] or A[0][1]) else [(A[1][1]) % P, (-A[1][0]) % P]
    assert all((A[i][0]*u[0] + A[i][1]*u[1]) % P == 0 for i in range(2)) and any(u)
    for k in (1, 2, 3):
        n = 2*k
        mons = [e for e in itertools.product(range(7), repeat=n) if sum(e) == 6]
        rows = []
        for trial in range(len(mons) + 20):
            up = []
            for _ in range(k):
                d = m2([random.randint(-5, 5) for _ in range(4)])
                up += [(d[0][0]*u[0] + d[0][1]*u[1]) % P, (d[1][0]*u[0] + d[1][1]*u[1]) % P]
            row = []
            for e in mons:
                v = 1
                for i in range(n): v = v * pow(up[i], e[i], P) % P
                row.append(v)
            rows.append(row)
        r = rank_mod(rows, P)
        line = f'a={a} c={c} q={qv} b={b} P={P} k={k}: dim Sym^6(U^k) = {len(mons)}, rank mod P of pullback span = {r}'
        print(line); out.append(line)
open(os.path.join(HERE, 'pullback_span_modp_out.txt'), 'w').write('\n'.join(out) + '\n')
