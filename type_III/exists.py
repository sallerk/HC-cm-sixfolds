# exists.py -- test "exists K in D with split K-Hermitian form  <=>  Q(sqrt Delta) embeds in D"
# (note, Proposition 6.2, second statement), and run the construction used in its proof: when
# Q(sqrt Delta) embeds, find pure p0 with p0^2 = Delta*s^2, take pure x anticommuting with p0
# (trace-orthogonal), K := Q(x), and check by the direct Gram computation (path1 code) that this K is split.
# Two embedding tests for Q(sqrt Delta) -> D: (i) local: Delta not a square in Q_p for all p | disc D;
# (ii) global: qfsolve finds an isotropic vector of <-a,-c,ac,Delta> (pure-norm form represents -Delta).
# Reads data.json, writes exists_out.json (both next to this script).
import json, sys, itertools, random, os
from fractions import Fraction as F
from cypari import pari
from path1 import mul, conj, add, scal, nrd, nrd_matrix, primes_of, hilb

HERE = os.path.dirname(os.path.abspath(__file__))

def is_sq_Qp(x, p):
    x = F(x); n = x.numerator * x.denominator
    v = 0
    while n % p == 0: n //= p; v += 1
    if v % 2: return False
    if p == 2: return n % 8 == 1
    return pow(n % p, (p-1)//2, p) == 1

def gram_split(T, q, a, c):
    """direct computation as in path1: returns (b, detH, split)"""
    b = nrd(q, a, c); qinv = scal(F(1)/b, conj(q))
    z = None
    for cand in ([0,1,0,0], [0,0,1,0], [0,0,0,1]):
        cand = [F(t) for t in cand]
        cr = [q[2]*cand[3]-q[3]*cand[2], q[3]*cand[1]-q[1]*cand[3], q[1]*cand[2]-q[2]*cand[1]]
        if any(cr): z = cand; break
    basis = [([F(1),F(0),F(0),F(0)], m) for m in range(3)] + [(z, m) for m in range(3)]
    rows = []
    for (cu, mu) in basis:
        row = []
        for (cv, mv) in basis:
            w = mul(mul(cu, T[mu][mv], a, c), conj(cv), a, c)
            piw = scal(F(1, 2), add(w, mul(mul(q, w, a, c), qinv, a, c)))
            al = piw[0]; rest = add(piw, [-al, 0, 0, 0]); be = mul(rest, conj(q), a, c)[0] / b
            row.append(f'Mod(({-b*be})+({al})*t,t^2+({b}))')
        rows.append(','.join(row))
    dt = pari(f'lift(matdet([{";".join(rows)}]))')
    detH = F(str(pari(f'polcoef({dt},0,t)')))
    P = primes_of(b, detH, a, c)
    return b, detH, (-detH > 0) and all(hilb(-b, -detH, p) == 1 for p in P)

def disc_primes(a, c):
    return [p for p in primes_of(a, c) if hilb(a, c, p) == -1]

def construct_K(a, c, Delta):
    """construction from the proof: pure p0 with Nrd(p0) = -Delta*s^2, then pure x orthogonal to p0."""
    M = f'matdiagonal([{-a},{-c},{a*c},{Delta}])'
    sol = pari(f'qfsolve({M})')
    if str(pari(f'type({sol})')) != 't_COL':
        return None
    x1, x2, x3, y = [F(str(t)) for t in sol]
    assert y != 0
    p0 = [F(0), x1/y, x2/y, x3/y]
    assert nrd(p0, a, c) == -F(Delta)
    # pure x with Trd(p0 * x) = 0 : bilinear form of Nrd on pure part: -a x1 y1 - c x2 y2 + ac x3 y3 = 0
    g = [-a*p0[1], -c*p0[2], a*c*p0[3]]
    ker = [v for v in ([g[1], -g[0], F(0)], [g[2], F(0), -g[0]], [F(0), g[2], -g[1]]) if any(v)]
    k1 = ker[0]; k2 = next(v for v in ker[1:] if any(k1[i]*v[j] - k1[j]*v[i] for i in range(3) for j in range(3)))
    cands = []
    for s, t in itertools.product(range(-3, 4), repeat=2):
        v = [s*k1[i] + t*k2[i] for i in range(3)]
        if any(v):
            assert sum(gi*vi for gi, vi in zip(g, v)) == 0
            cands.append([F(0)] + v)
    x = min(cands, key=lambda w: abs(nrd(w, a, c)))
    # check anticommutation
    assert mul(p0, x, a, c) == scal(-1, mul(x, p0, a, c))
    return p0, x

def main():
    data = json.load(open(os.path.join(HERE, 'data.json')))
    res = []
    stats = {'emb_and_found_split': 0, 'emb_no_split_in_list_but_constructed': 0, 'emb_construct_FAIL': 0,
             'notemb_and_no_split': 0, 'notemb_but_split_FOUND': 0, 'local_vs_qfsolve_mismatch': 0}
    for A in data['algebras']:
        a, c, disc = A['a'], A['c'], A['disc']
        dp = disc_primes(a, c); assert eval('*'.join(map(str, dp))) == disc
        for ti, Trec in enumerate(A['Ts']):
            T = [[[F(t) for t in e] for e in row] for row in Trec['T']]
            N = nrd_matrix(T, a, c)
            if N == 0: continue
            Delta = -N
            emb_local = all(not is_sq_Qp(Delta, p) for p in dp)
            sol = pari(f'qfsolve(matdiagonal([{-a},{-c},{a*c},{Delta}]))')
            emb_global = str(pari(f'type({sol})')) == 't_COL'
            if emb_local != emb_global: stats['local_vs_qfsolve_mismatch'] += 1
            splits = []
            for qv in A['qs']:
                q = [F(t) for t in qv]
                b, dH, sp = gram_split(T, q, a, c)
                if sp: splits.append(str(b))
            rec = {'disc': disc, 'ti': ti, 'Delta': str(Delta), 'emb': emb_local, 'n_split_listed': len(splits),
                   'split_b': sorted(set(splits), key=lambda s: F(s))[:10]}
            if emb_local:
                if splits: stats['emb_and_found_split'] += 1
                ck = construct_K(a, c, Delta)
                p0, x = ck
                b, dH, sp = gram_split(T, x, a, c)
                rec['constructed_b'] = str(b); rec['constructed_split'] = sp
                if not sp: stats['emb_construct_FAIL'] += 1
                elif not splits: stats['emb_no_split_in_list_but_constructed'] += 1
            else:
                if splits: stats['notemb_but_split_FOUND'] += 1
                else: stats['notemb_and_no_split'] += 1
            res.append(rec)
        print('done', disc, file=sys.stderr)
    json.dump({'stats': stats, 'records': res}, open(os.path.join(HERE, 'exists_out.json'), 'w'), indent=0)
    print(stats)

if __name__ == '__main__':
    main()
