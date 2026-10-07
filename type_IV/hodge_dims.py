# hodge_dims.py -- dimensions of Hodge classes on A^n for the VERY GENERAL member A of
# the family (Hg = ker(det_K) in Lef = R_{F/Q} U(H); over C: Hg = {(g1,g2) in GL3 x GL3 : det g1 det g2 = 1}),
# compared with Lefschetz classes (Lef-invariants).  Method: Weyl character formula (multiplicity of
# det^j (x) det^j in the graded character of  /\*( (St1 + St1^v + St2 + St2^v) (x) C^n )).
#   dim B^k(A^n) = sum_j  sum_{k1+k2=2k} c_n(j, k1) c_n(j, k2),  c_n(j,k) = mult of det^j in degree k of
#   /\*((St + St^v) (x) C^n)  for GL3.   j = 0 gives the Lefschetz classes D^k(A^n) (Milne 1999).
import sys, itertools
from collections import defaultdict

def poly_mul(P, Q):
    R = defaultdict(int)
    for k1, v1 in P.items():
        for k2, v2 in Q.items():
            R[tuple(a + b for a, b in zip(k1, k2))] += v1 * v2
    return {k: v for k, v in R.items() if v}

def gl3_char(n):
    # monomials (q, e1, e2, e3): prod_a ((1 + q x_a)(1 + q / x_a))^n
    P = {(0, 0, 0, 0): 1}
    for a in range(3):
        f = {(0, 0, 0, 0): 1}
        for _ in range(n):
            e = [0, 0, 0, 0]; e[1 + a] = 1
            g = {(0, 0, 0, 0): 1, (1,) + tuple(e[1:]): 1}
            e2 = [0, 0, 0]; e2[a] = -1
            h = {(0, 0, 0, 0): 1, (1,) + tuple(e2): 1}
            f = poly_mul(poly_mul(f, g), h)
        P = poly_mul(P, f)
    return P

def mults(n):
    P = gl3_char(n)
    # multiply by Vandermonde prod_{a<b}(x_a - x_b)
    V = {(0, 0, 0, 0): 1}
    for a, b in [(0, 1), (0, 2), (1, 2)]:
        ea = [0, 0, 0, 0]; ea[1 + a] = 1
        eb = [0, 0, 0, 0]; eb[1 + b] = 1
        V = poly_mul(V, {tuple(ea): 1, tuple(eb): -1})
    PV = poly_mul(P, V)
    c = defaultdict(int)      # c[(j, k)]
    for (k, e1, e2, e3), v in PV.items():
        # coefficient of x^{lambda + rho}, rho = (2,1,0), lambda = (j,j,j)
        if e1 - 2 == e2 - 1 == e3:
            c[(e3, k)] += v
    return c

def table(n):
    c = mults(n)
    js = sorted(set(j for j, _ in c))
    out = []
    for k in range(0, 6 * n + 1):            # H^{2k}(A^n), dim A^n = 6n
        B = 0; D = 0
        for j in js:
            s = sum(c.get((j, k1), 0) * c.get((j, 2 * k - k1), 0) for k1 in range(0, 2 * k + 1))
            B += s
            if j == 0:
                D += s
        out.append((k, B, D, B - D))
    return out, c

if __name__ == '__main__':
    for n in [1, 2, 3]:
        out, c = table(n)
        print('n = %d  (A^%d, dim %d):  k, dim B^k, dim D^k, extra' % (n, n, 6 * n))
        for row in out:
            print('   ', row)
        print('   nonzero det-weights j:', sorted(set(j for (j, k), v in c.items() if v)))
