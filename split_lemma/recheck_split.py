# recheck_split.py -- independent check (no PARI): re-checks every witness in split_witnesses.json
# (output of verify_split.py, read from the directory of this script).
# For each case and witness a = (a_i in F_i = Q[x]/f_i):
#   1. signs of a at the real roots of f_i (mpmath, 60 digits) -> sign vector has exactly n minus signs,
#      and the set of witnessed sign vectors covers ALL C(2n, n) patterns (place ordering may differ
#      from PARI's, so we compare as a set after a per-factor root sort, then check coverage);
#   2. N(a_i) = Res(f_i, a_i) (sympy, exact; f_i monic) and agreement with the stored N(a);
#   3. delta' = (-1)^n N(a) prod disc(f_i)  (polynomial discriminants: same class mod squares as d_F);
#      delta'/delta is a rational square;
#   4. exact identity X^2 + d Y^2 = delta Z^2 with Z != 0  (so delta in Nm(Q(sqrt -d)^x));
#   5. determinant formula: det( Tr_{F/Q}(a b_j b_k) ) = N(a) * disc(power basis) (exact, sympy, power
#      basis of each factor; block-diagonal for products) for the first witness of each case;
#   6. Hilbert symbols t_p = ((-1)^n dF, -d)_p recomputed with an independent implementation.
import json, os, itertools, sys
from fractions import Fraction
import sympy as sp
import mpmath as mp
mp.mp.dps = 60
x = sp.symbols('x')

def is_rat_square(q):
    q = Fraction(q)
    if q <= 0: return False
    a, b = q.numerator, q.denominator
    return sp.integer_nthroot(a, 2)[1] and sp.integer_nthroot(b, 2)[1]

def vp(n, p):
    k = 0
    while n % p == 0:
        n //= p; k += 1
    return k, n

def hilbert(a, b, p):
    # Serre, A Course in Arithmetic, Ch. III Thm 1; a, b nonzero integers
    if p == 0:  # infinity
        return -1 if (a < 0 and b < 0) else 1
    al, u = vp(a, p); be, v = vp(b, p)
    if p != 2:
        leg = lambda t: sp.legendre_symbol(t % p, p)
        eps = (p - 1) // 2
        s = (-1) ** (al * be * eps)
        return s * leg(u) ** be * leg(v) ** al
    e = lambda t: ((t - 1) // 2) % 2
    w = lambda t: ((t * t - 1) // 8) % 2
    return (-1) ** ((e(u) * e(v) + al * w(v) + be * w(u)) % 2)

data = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'split_witnesses.json')))
tot_w = 0; bad = []
for case in data:
    n, d = case['n'], case['d']
    fs = [sp.Poly(sp.sympify(p.replace('^', '**')), x) for p in case['polys']]
    roots = [sorted([mp.mpf(r) for r in mp.polyroots([mp.mpf(int(c)) for c in f.all_coeffs()],
                                                     maxsteps=400, extraprec=400)]) if f.degree() > 1
             else [mp.mpf(-int(f.all_coeffs()[1]))] for f in fs]
    discs = [int(sp.discriminant(f.as_expr(), x)) if f.degree() > 1 else 1 for f in fs]
    pats = set()
    first = True
    for pat in case['patterns']:
        wit = pat['witness']
        if wit is None:
            bad.append((case['label'], d, 'no witness', pat['signs'])); continue
        tot_w += 1
        als = [sp.Poly(sp.sympify(a.replace('^', '**')), x) for a in wit['a_alg']]
        sv = []
        N = sp.Integer(1)
        for f, a, rts in zip(fs, als, roots):
            vals = [a.eval(sp.Float(str(r), 60)) if a.degree() >= 0 else 0 for r in rts]
            sv += [1 if v > 0 else -1 for v in vals]
            N *= (sp.resultant(f.as_expr(), a.as_expr(), x) if f.degree() > 1
                  else sp.sympify(a.as_expr()).subs(x, sp.Integer(int(mp.nint(rts[0])))))
        if sv.count(-1) != n:
            bad.append((case['label'], d, 'sign count', sv))
        pats.add(tuple(sv))
        if sp.Integer(int(wit['N_a'])) != N:
            bad.append((case['label'], d, 'norm mismatch', str(N), wit['N_a']))
        delta = Fraction(int(wit['delta']))
        dprime = Fraction((-1) ** n * int(N))
        for D in discs: dprime *= D
        if not is_rat_square(dprime / delta):
            bad.append((case['label'], d, 'delta class mismatch'))
        X, Y, Z = [Fraction(sp.Rational(c)) for c in wit['XYZ']]
        if Z == 0 or X * X + d * Y * Y != delta * Z * Z:
            bad.append((case['label'], d, 'norm equation fails'))
        if first:
            first = False
            # determinant of the twisted trace form on the power bases (block diagonal over factors)
            det = sp.Integer(1); prodN = sp.Integer(1); prodD = sp.Integer(1)
            for f, a in zip(fs, als):
                k = f.degree()
                if k == 1:
                    val = a.as_expr().subs(x, -f.all_coeffs()[1]); det *= val; prodN *= val; continue
                # Tr(g) for g in Q[x]/f via companion matrix trace
                Cm = sp.Matrix(k, k, lambda i, j: 0)
                co = f.all_coeffs()[::-1]  # c0..ck
                for i in range(1, k): Cm[i, i - 1] = 1
                for i in range(k): Cm[i, k - 1] = -co[i]
                def tr(g):
                    M = sp.zeros(k, k); P = sp.eye(k)
                    for c in sp.Poly(g, x).all_coeffs()[::-1]:
                        M += c * P; P = P * Cm
                    return M.trace()
                G = sp.Matrix(k, k, lambda i, j: tr(sp.rem(a.as_expr() * x ** (i + j), f.as_expr(), x)))
                det *= G.det(); prodN *= sp.resultant(f.as_expr(), a.as_expr(), x)
                prodD *= sp.discriminant(f.as_expr(), x)
            if det != prodN * prodD:
                bad.append((case['label'], d, 'det formula fails', str(det), str(prodN * prodD)))
    allp = set(s for s in itertools.product((1, -1), repeat=2 * n) if s.count(-1) == n)
    if pats != allp:
        bad.append((case['label'], d, f'pattern coverage {len(pats)}/{len(allp)}'))
    dF = int(case['dF'])
    for r in case['forced']:
        h = hilbert((-1) ** n * dF, -d, r['p'])
        if h != r['t_p']:
            bad.append((case['label'], d, 'hilbert mismatch', r))
print('cases', len(data), 'witnesses rechecked', tot_w, 'problems', len(bad))
for b in bad[:40]: print('  ', b)
