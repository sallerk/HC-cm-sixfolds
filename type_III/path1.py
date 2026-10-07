# path1.py -- main computation (uses PARI). Quaternion structure constants over Q (Fractions).
# For each (D,T,q): K = Q(q), q^2=-b. K-Hermitian form H(u,v) = q * pi_K(T(u,v)),
# pi_K(w) = (w + q w q^{-1})/2, on the K-basis {e_m, z e_m} (z in D \ K) of V = D^3
# (left D-module, T(x,y) = sum x_r T_rs conj(y_s)).
# det_K(H) via PARI over K; split test: -det(H) in Nm(K^*) via Hilbert symbols (-b,-det)_p.
# Delta = -Nrd(T), Nrd via embedding D -> M_2(Q(sqrt a)).
# Criterion tested (note, Proposition 6.2): split  <=>  (-b, Delta)_p = (a,c)_p for all p.
# Reads data.json, writes path1_out.json (both next to this script).
import json, sys, os
from fractions import Fraction as F
from cypari import pari

HERE = os.path.dirname(os.path.abspath(__file__))

def mul(x, y, a, c):
    x0, x1, x2, x3 = x; y0, y1, y2, y3 = y
    return [x0*y0 + a*x1*y1 + c*x2*y2 - a*c*x3*y3,
            x0*y1 + x1*y0 - c*x2*y3 + c*x3*y2,
            x0*y2 + x2*y0 + a*x1*y3 - a*x3*y1,
            x0*y3 + x3*y0 + x1*y2 - x2*y1]
def conj(x): return [x[0], -x[1], -x[2], -x[3]]
def add(x, y): return [s+t for s, t in zip(x, y)]
def scal(s, x): return [s*t for t in x]
def nrd(x, a, c): return x[0]**2 - a*x[1]**2 - c*x[2]**2 + a*c*x[3]**2

def primes_of(*nums):
    S = {2}
    for n in nums:
        n = F(n)
        for m in (n.numerator, n.denominator):
            m = abs(m)
            if m > 1:
                for p in pari(f'factor({m})[,1]~'): S.add(int(p))
    return sorted(S)

def hilb(x, y, p):
    return int(pari(f'hilbert({x},{y},{p})'))

def nrd_matrix(T, a, c):
    # embed D -> M_2(Q(s)), s^2=a: i->[[s,0],[0,-s]], j->[[0,c],[1,0]], k->[[0,s*c],[-s,0]]
    rows = []
    for r in range(3):
        r0, r1 = [], []
        for s_ in range(3):
            x0, x1, x2, x3 = [str(F(t)) for t in T[r][s_]]
            m00 = f'({x0})+({x1})*s'; m01 = f'({x2})*{c}+({x3})*s*{c}'
            m10 = f'({x2})-({x3})*s'; m11 = f'({x0})-({x1})*s'
            r0 += [m00, m01]; r1 += [m10, m11]
        rows.append(r0); rows.append(r1)
    M = '[' + ';'.join(','.join(f'Mod({e},s^2-({a}))' for e in row) for row in rows) + ']'
    dt = pari(f'lift(matdet({M}))')
    assert int(pari(f'poldegree({dt},s)')) <= 0, dt
    return F(str(pari(f'polcoef({dt},0,s)')))

def run():
    data = json.load(open(os.path.join(HERE, 'data.json')))
    out = []
    for A in data['algebras']:
        a, c, disc = A['a'], A['c'], A['disc']
        for ti, Trec in enumerate(A['Ts']):
            T = [[[F(t) for t in e] for e in row] for row in Trec['T']]
            N = nrd_matrix(T, a, c)
            if N == 0:
                continue
            Delta = -N
            for q in A['qs']:
                q = [F(t) for t in q]
                b = nrd(q, a, c)          # q^2 = -b
                qinv = scal(F(1)/b, conj(q))
                # z not in K = Q + Qq
                z = None
                for cand in ([0,1,0,0], [0,0,1,0], [0,0,0,1]):
                    cand = [F(t) for t in cand]
                    # cand in K iff cand proportional to q (both pure)
                    cr = [q[2]*cand[3]-q[3]*cand[2], q[3]*cand[1]-q[1]*cand[3], q[1]*cand[2]-q[2]*cand[1]]
                    if any(cr): z = cand; break
                basis = [([F(1),0,0,0], m) for m in range(3)] + [(z, m) for m in range(3)]
                basis = [([F(t) for t in co], m) for co, m in basis]
                Hent = []
                for (cu, mu) in basis:
                    row = []
                    for (cv, mv) in basis:
                        w = mul(mul(cu, T[mu][mv], a, c), conj(cv), a, c)
                        piw = scal(F(1, 2), add(w, mul(mul(q, w, a, c), qinv, a, c)))
                        al = piw[0]
                        rest = add(piw, [-al, 0, 0, 0])      # = beta*q
                        be = mul(rest, conj(q), a, c)[0] / b
                        # check rest == be*q
                        assert all(rest[t] == be*q[t] for t in range(4))
                        # H = q*(al + be q) = -b*be + al*q   ; q -> t, t^2=-b
                        row.append((-b*be, al))
                    Hent.append(row)
                bnum = b  # rational
                M = '[' + ';'.join(','.join(f'Mod(({x})+({y})*t,t^2+({bnum}))' for (x, y) in row) for row in Hent) + ']'
                dt = pari(f'lift(matdet({M}))')
                assert int(pari(f'poldegree({dt},t)')) <= 0, ('det not rational', dt)
                detH = F(str(pari(f'polcoef({dt},0,t)')))
                assert detH != 0
                # split <=> -detH is a norm from Q(sqrt(-b)) <=> (-b,-detH)_v = 1 all v
                P = primes_of(b, detH, Delta, a, c)
                split = (-detH > 0) and all(hilb(-b, -detH, p) == 1 for p in P)
                pred = all(hilb(-b, Delta, p) == hilb(a, c, p) for p in P)
                # detH and NrdT are recorded; the determinant formula of Lemma 6.1 is tested in check_lemma2.py
                out.append({'disc': disc, 'a': a, 'c': c, 'ti': ti, 'kind': Trec['kind'], 'q': [str(t) for t in q],
                            'b': str(b), 'detH': str(detH), 'NrdT': str(N), 'Delta': str(Delta),
                            'split_hilbert': split, 'pred_formula': pred, 'sign_ok': detH < 0})
        print('done disc', disc, file=sys.stderr)
    json.dump(out, open(os.path.join(HERE, 'path1_out.json'), 'w'))
    n = len(out); agree = sum(r['split_hilbert'] == r['pred_formula'] for r in out)
    nsplit = sum(r['split_hilbert'] for r in out)
    print(f'path1: {n} (D,T,K) cases; split {nsplit}; formula agrees {agree}/{n}; detH<0 in {sum(r["sign_ok"] for r in out)}')

if __name__ == '__main__':
    run()
