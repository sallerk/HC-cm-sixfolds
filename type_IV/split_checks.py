# split_checks.py -- tests of the SPLIT LEMMA (note, Sec. 7) on random general L-Hermitian forms.
#
#  LEMMA.  L = K F (K = Q(sqrt-d), F = Q(sqrt m)), V = L^3, H L-Hermitian with signature (2,1) at f1 and
#  (1,2) at f2.  Then some totally positive a in F makes h_a = Tr_{L/K}(a H) split (hyperbolic) over K.
#
#  Main computation A (proof 1, local-global): c = -N(Delta) d_F; at every "forced" prime p (-dm square in Q_p,
#     -d not) check (c,-d)_p = 1 and its two ingredients (N Delta,-d)_p = 1, (-d_F,-d)_p = 1; then solve
#     u^2 + d v^2 - c x^2 + c m y^2 = 0 with PARI qfsolve (Simon) and put a = +-(x + y sqrt m).
#  Main computation B (proof 2, explicit): find an H-isotropic e (diagonalise H, solve a relative norm equation
#     with rnfisnorm), complete to a hyperbolic L-plane P = <e, f'>, g spans P^perp, delta = H(g,g);
#     a = sqrt(m)/delta.  Verify W = L e + K g is a 3-dim totally h_a-isotropic K-subspace.
#  Checker C (independent of A and B): explicit 6x6 Gram matrix of h_a over K from H and a;
#     (i) det == N(a)^3 N(Delta) (4m)^3 exactly; (ii) -det is a norm from K (Hilbert symbols);
#     (iii) Witt index of the 12-dim rational form x -> h_a(x,x) by iterated qfsolve == 6 (hyperbolic).
#     Control: a = 1 (split or not), consistency of (ii) and (iii).
#  Product (degenerate F = Q x Q) case and an obstruction control: see degenerate_cases.py.
import os, random, itertools, json, sys, time
sys.set_int_max_str_digits(0)
from fractions import Fraction as Fr
from cypari import pari
from biq import Ctx, witt_index

HERE = os.path.dirname(os.path.abspath(__file__))   # JSON output is written next to this script
rng = random.Random(51005)

def rand_herm(C):
    H = [[None]*3 for _ in range(3)]
    for i in range(3):
        H[i][i] = C.Fel(rng.randint(-5, 5), rng.randint(-5, 5))
        for j in range(i + 1, 3):
            z = C.L(*[rng.randint(-3, 3) for _ in range(4)]) if rng.random() < 0.8 else C.L()
            H[i][j] = z; H[j][i] = C.conj(z)
    return H

def good_herm(C, tries=200):
    for _ in range(tries):
        H = rand_herm(C)
        try:
            D = C.detL(H)
        except AssertionError:
            continue
        if C.NF(D) == 0:
            continue
        s1, s2 = C.signature(H, 1), C.signature(H, 2)
        if (s1, s2) == ((2, 1), (1, 2)):
            return H
        if (s1, s2) == ((1, 2), (2, 1)):
            return [[C.rho(H[i][j]) for j in range(3)] for i in range(3)]   # relabel the two real places
    return None

# ---------------- Main computation A ----------------
def compute_A(C, H):
    d, m = C.d, C.m
    D = C.detL(H); ND = C.NF(D)
    assert ND < 0, 'N(Delta) must be negative for signatures (2,1),(1,2)'
    c = -ND * C.dF
    forced = C.forced_primes(c)
    fchk = []
    for p in forced:
        e1 = C.hilbert(ND, -d, p); e2 = C.hilbert(Fr(-C.dF), -d, p); e = C.hilbert(c, -d, p)
        fchk.append((p, e1, e2, e))
        if not (e1 == e2 == e == 1):
            return dict(ok=False, why='forced prime %d symbols %s' % (p, (e1, e2, e)))
    den = c.denominator
    Qm = pari.matdiagonal([den*den, den*den*d, -c.numerator*den, c.numerator*den*m])
    v = pari.qfsolve(Qm)
    if str(pari.type(v)) == 't_MAT':
        v = v[0]
    if str(pari.type(v)) != 't_COL':
        return dict(ok=False, why='qfsolve obstruction %s' % v, forced=fchk)
    u, vv, x, y = [Fr(int(pari.numerator(t)), int(pari.denominator(t))) for t in v]
    a = C.Fel(x, y)
    assert x != 0 or y != 0
    if not C.totpos(a):
        a = C.neg(a)
    assert C.totpos(a)
    # claim of computation A: c N(a) = u'^2 + d v'^2 (with the scaling by den)
    return dict(ok=True, a=a, forced=fchk, c=c)

# ---------------- Main computation B ----------------
def lin_comb(C, coeffs, vecs):
    out = [C.L()]*3
    for c, v in zip(coeffs, vecs):
        out = [C.add(o, C.mul(c, vi)) for o, vi in zip(out, v)]
    return out

def gram_schmidt(C, H):
    E = [[C.L(1) if i == j else C.L() for j in range(3)] for i in range(3)]
    basis, deltas = [], []
    for v in E:
        w = v
        for b, db in zip(basis, deltas):
            coef = C.div(C.Hform(H, w, b), db)
            w = [C.sub(wi, C.mul(coef, bi)) for wi, bi in zip(w, b)]
        dw = C.Hform(H, w, w)
        if dw == C.L():
            return None, w          # found an isotropic vector directly
        basis.append(w); deltas.append(dw)
    return (basis, deltas), None

def small_L(B):
    rng_ = range(-B, B + 1)
    els = [Ctx.L(*t) for t in itertools.product(rng_, repeat=4)]
    els.sort(key=lambda z: sum(abs(c) for c in z))
    return els

SMALL = small_L(2)

def isotropic_vector(C, H):
    gs, iso = gram_schmidt(C, H)
    if iso is not None:
        return iso, 'gram-schmidt pivot'
    basis, deltas = gs
    for (i, j, k) in [(0, 1, 2), (1, 0, 2), (0, 2, 1), (2, 0, 1), (1, 2, 0), (2, 1, 0)]:
        for x2 in SMALL:
            n2 = C.mul(x2, C.conj(x2))
            t = C.neg(C.div(C.add(deltas[k], C.mul(deltas[j], n2)), deltas[i]))
            if t == C.L():
                continue
            if not C.totpos(t):
                continue
            # cheap necessary condition before the (sometimes very slow) rnfisnorm call:
            # t = N_{L/F}(x) implies N_{F/Q}(t) = N_{K/Q}(N_{L/K} x) is a norm from K; also skip huge heights
            Nt = C.NF(t)
            if abs(Nt.numerator) * Nt.denominator > 10**40 or not C.is_norm_K(Nt):
                continue
            x1 = C.norm_LF_solve(t)
            if x1 is None:
                continue
            coeffs = [None]*3
            coeffs[i], coeffs[j], coeffs[k] = x1, x2, C.L(1)
            e = lin_comb(C, coeffs, basis)
            assert C.Hform(H, e, e) == C.L()
            return e, 'norm equation'
    return None, 'not found'

def kcoords(v):
    out = []
    for z in v:
        out += [(z[0], z[1]), (z[2], z[3])]
    return out

def krank(C, vecs):
    rows = [','.join('Mod(%s+(%s)*t,t^2+%d)' % (str(x), str(y), C.d) for x, y in kcoords(v)) for v in vecs]
    return int(pari.matrank(pari('[' + ';'.join(rows) + ']')))

def cross(C, a, b):
    return [C.sub(C.mul(a[1], b[2]), C.mul(a[2], b[1])),
            C.sub(C.mul(a[2], b[0]), C.mul(a[0], b[2])),
            C.sub(C.mul(a[0], b[1]), C.mul(a[1], b[0]))]

def compute_B(C, H):
    e, how = isotropic_vector(C, H)
    if e is None:
        return dict(ok=False, why='no isotropic vector found in search box')
    j = next(j for j in range(3) if C.Hform(H, e, [C.L(1) if k == j else C.L() for k in range(3)]) != C.L())
    ej = [C.L(1) if k == j else C.L() for k in range(3)]
    lam = C.conj(C.inv(C.Hform(H, e, ej)))
    f = [C.mul(lam, z) for z in ej]
    assert C.Hform(H, e, f) == C.L(1)
    hff = C.Hform(H, f, f)
    half = tuple(x / 2 for x in hff)
    fp = [C.sub(fi, C.mul(half, ei)) for fi, ei in zip(f, e)]
    assert C.Hform(H, fp, fp) == C.L() and C.Hform(H, e, fp) == C.L(1)
    # g with H(g,e) = H(g,f') = 0: H(g,x) = sum_i g_i (H conj(x))_i  -> g = cross(H conj(e), H conj(f'))
    ce = [C.L()]*3; cf = [C.L()]*3
    for i in range(3):
        for jj in range(3):
            ce[i] = C.add(ce[i], C.mul(H[i][jj], C.conj(e[jj])))
            cf[i] = C.add(cf[i], C.mul(H[i][jj], C.conj(fp[jj])))
    g = cross(C, ce, cf)
    assert C.Hform(H, g, e) == C.L() and C.Hform(H, g, fp) == C.L()
    delta = C.Hform(H, g, g)
    assert delta[1] == 0 and delta[3] == 0 and delta != C.L()
    s1, s2 = C.realF(delta, 1), C.realF(delta, 2)
    mixed = s1 > 0 > s2
    a = C.div(C.Fel(0, 1), delta)               # a = sqrt(m)/delta
    tp = C.totpos(a)
    we = [C.mul(C.Fel(0, 1), z) for z in e]
    W = [e, we, g]
    iso = all(C.trLK(C.mul(a, C.Hform(H, x, y))) == (0, 0) for x in W for y in W)
    rk = krank(C, W)
    return dict(ok=bool(mixed and tp and iso and rk == 3), a=a, mixed=mixed, totpos=tp, W_isotropic=iso,
                W_rank=rk, how=how)

# ---------------- Checker C ----------------
def checker(C, H, a):
    G = C.gram_K(H, a)
    det = C.detK(G)
    formula = C.NF(a)**3 * C.NF(C.detL(H)) * (4*C.m)**3
    hil = C.is_norm_K(-det)
    wi, why = witt_index(C.gram_Q(G))
    return dict(det_ok=(det == formula), hilbert_split=hil, witt=wi, witt_reason=why)

def F2(x): return str(x)

def run():
    pairs = [(d, m) for d in [1, 2, 3, 5, 6, 7, 11, 15, 19, 23, 31, 35]
             for m in [2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 21, 30]]
    per = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    # optional chunking: argv[2] = chunk index i, argv[3] = number of chunks k; chunk i takes pairs[i::k]
    # and reseeds rng with 51005 + i (so H's differ from the unchunked run; all pairs still covered)
    tag = ''
    if len(sys.argv) > 3:
        i, k = int(sys.argv[2]), int(sys.argv[3])
        pairs = pairs[i::k]
        rng.seed(51005 + i)
        tag = '_%d' % i
    stats = dict(cases=0, A_ok=0, B_ok=0, C_A_ok=0, C_B_ok=0, h1_split=0, h1_consistent=0,
                 forced_cases=0, forced2_cases=0, B_notfound=0, det_formula_ok=0)
    fails = []
    examples = []
    t0 = time.time()
    for (d, m) in pairs:
        C = Ctx(d, m)
        for _ in range(per):
            H = good_herm(C)
            if H is None:
                continue
            stats['cases'] += 1
            A = compute_A(C, H)
            B = compute_B(C, H)
            rec = dict(d=d, m=m)
            if A['ok']:
                stats['A_ok'] += 1
                if A['forced']:
                    stats['forced_cases'] += 1
                    if any(p == 2 for p, *_ in A['forced']):
                        stats['forced2_cases'] += 1
                cA = checker(C, H, A['a'])
                stats['det_formula_ok'] += cA['det_ok']
                if cA['det_ok'] and cA['hilbert_split'] and cA['witt'] == 6:
                    stats['C_A_ok'] += 1
                else:
                    fails.append(('C_A', d, m, cA))
            else:
                fails.append(('A', d, m, A.get('why')))
            if B['ok']:
                stats['B_ok'] += 1
                cB = checker(C, H, B['a'])
                if cB['det_ok'] and cB['hilbert_split'] and cB['witt'] == 6:
                    stats['C_B_ok'] += 1
                else:
                    fails.append(('C_B', d, m, cB))
            else:
                if B.get('why', '').startswith('no isotropic'):
                    stats['B_notfound'] += 1
                else:
                    fails.append(('B', d, m, B))
            c1 = checker(C, H, C.Fel(1))
            stats['h1_split'] += int(c1['hilbert_split'])
            stats['h1_consistent'] += int(c1['hilbert_split'] == (c1['witt'] == 6))
            print('case', stats['cases'], d, m, 'A', A['ok'], 'B', B['ok'], 'h1', c1['witt'], round(time.time() - t0), flush=True)
            if len(examples) < 12 and A['ok'] and B['ok']:
                examples.append(dict(d=d, m=m, H=[[[F2(x) for x in H[i][j]] for j in range(3)] for i in range(3)],
                                     a_A=[F2(A['a'][0]), F2(A['a'][2])], a_B=[F2(B['a'][0]), F2(B['a'][2])],
                                     forced=A['forced'], h1_witt=c1['witt']))
    stats['seconds'] = round(time.time() - t0, 1)
    print('STATS', json.dumps(stats))
    print('FAILS', len(fails))
    for f in fails[:40]:
        print('  ', f)
    json.dump(dict(stats=stats, fails=[str(f) for f in fails], examples=examples),
              open(os.path.join(HERE, 'split_checks%s.json' % tag), 'w'), indent=1)

if __name__ == '__main__':
    run()
