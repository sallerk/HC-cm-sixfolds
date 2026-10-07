"""Search totally real sextic fields F with each Galois group 6Tj (j = 1..16), PARI polgalois (exact for deg 6).
Families: (A) even g(x^2); (B) norm forms A^2 - d B^2 (quadratic subfield); (C) random dense;
(D) resolvent sextics: sum-of-pairs of totally real quartics (6T4/6T7/6T8), Cayley-type sextic of totally
real quintics (6T12/6T14). Keeps the smallest |disc| (polredabs) per j. Output: search_F.json."""
import json, random, itertools, time, os
from cypari import pari
pari.allocatemem(2*10**9)
pari('default(new_galois_format,1)')
best = {}
def consider(f, tag):
    f = pari(f)
    if pari.poldegree(f) != 6 or not pari.polisirreducible(f) or int(pari.polsturm(f)) != 6:
        return
    try:
        f = pari.polredabs(f)
    except Exception:
        return
    D = abs(int(pari.poldisc(f)))  # = field disc for polredabs only if monogenic; use nfdisc below
    D = abs(int(pari.nfdisc(f)))
    j = int(pari.polgalois(f)[2])
    if j not in best or D < best[j][0]:
        best[j] = (D, str(f), tag)
random.seed(1)
t0 = time.time()
x = pari('x')
# (A) even sextics x^6 - a x^4 + b x^2 - c
for a in range(1, 25):
    for b in range(1, 60):
        for c in range(1, 40):
            consider(f'x^6-{a}*x^4+{b}*x^2-{c}', 'even')
print('A done', time.time() - t0, sorted((j, v[0]) for j, v in best.items()), flush=True)
# (B) norm forms
for d in (2, 3, 5, 6, 7, 13, 17):
    for _ in range(40000):
        A = [random.randint(-6, 6) for _ in range(3)]
        B = [random.randint(-3, 3) for _ in range(3)]
        if not any(B): continue
        f = f'(x^3+({A[0]})*x^2+({A[1]})*x+({A[2]}))^2-{d}*(({B[0]})*x^2+({B[1]})*x+({B[2]}))^2'
        consider(f, 'norm%d' % d)
print('B done', time.time() - t0, sorted((j, v[0]) for j, v in best.items()), flush=True)
# (C) random dense
for _ in range(150000):
    c = [random.randint(-8, 8) for _ in range(6)]
    consider('x^6+' + '+'.join(f'({c[i]})*x^{i}' for i in range(6)), 'dense')
print('C done', time.time() - t0, sorted((j, v[0]) for j, v in best.items()), flush=True)
# (D) resolvents (numerical roots, exact rounding checked by integrality tolerance)
pari.set_real_precision(120)
def resolvent(roots, exprs):
    vals = [e(roots) for e in exprs]
    P = pari(1)
    for v in vals:
        P = P * (x - v)
    co = [P.polcoef(k) for k in range(7)]
    rc = [pari.real(c).round() for c in co]
    if max(abs(float(pari.real(c) - r)) + abs(float(pari.imag(c))) for c, r in zip(co, rc)) > 1e-40 * (1 + max(abs(float(r)) for r in rc)):
        return None
    return sum(rc[k] * x**k for k in range(7))
pairs = list(itertools.combinations(range(4), 2))
for _ in range(20000):
    c = [random.randint(-12, 12) for _ in range(4)]
    g = pari('x^4+' + '+'.join(f'({c[i]})*x^{i}' for i in range(4)))
    if not pari.polisirreducible(g) or int(pari.polsturm(g)) != 4: continue
    r = pari.polroots(g, precision=600)
    R = resolvent(r, [lambda r, i=i, j=j: r[i] + r[j] for i, j in pairs])
    if R is not None: consider(R, 'quartic-pairs')
    # C4-invariant sextic: x1 x2^2 + x2 x3^2 + x3 x4^2 + x4 x1^2 over S4 cosets of <(1234)>
    exprs = []
    for perm in itertools.permutations(range(4)):
        if perm[0] != 0: continue
        exprs.append(lambda r, p=perm: sum(r[p[k]] * r[p[(k + 1) % 4]]**2 for k in range(4)))
    R = resolvent(r, exprs)
    if R is not None: consider(R, 'quartic-C4')
print('D4 done', time.time() - t0, sorted((j, v[0]) for j, v in best.items()), flush=True)
# quintic -> sextic: F20-invariant theta = (x1x2+x2x3+x3x4+x4x5+x5x1 - x1x3-x3x5-x5x2-x2x4-x4x1)^2
# cosets of F20 in S5: 6 cyclic orderings up to rotation/reflection
cyc = []
for perm in itertools.permutations(range(1, 5)):
    c5 = (0,) + perm
    key = min(min(tuple(c5[(s + k) % 5] for k in range(5)), tuple(c5[(s - k) % 5] for k in range(5))) for s in range(5))
    if key not in cyc: cyc.append(key)
assert len(cyc) == 12
def th(r, c):
    adj = sum(r[c[k]] * r[c[(k + 1) % 5]] for k in range(5))
    nonadj = sum(r[c[k]] * r[c[(k + 2) % 5]] for k in range(5))
    return (adj - nonadj)**2
# cyc has 12 = 5-cycles up to rotation and reflection... S5/D5 = 12; theta is F20-invariant so 6 distinct values
for _ in range(20000):
    c = [random.randint(-15, 15) for _ in range(5)]
    g = pari('x^5+' + '+'.join(f'({c[i]})*x^{i}' for i in range(5)))
    if not pari.polisirreducible(g) or int(pari.polsturm(g)) != 5: continue
    r = pari.polroots(g, precision=600)
    vals = [th(r, cc) for cc in cyc]
    uniq = []
    for v in vals:
        if not any(abs(complex(v - u)) < 1e-40 * (1 + abs(complex(u))) for u in uniq): uniq.append(v)
    if len(uniq) != 6: continue
    R = resolvent(r, [lambda r, v=v: v for v in uniq])
    if R is not None: consider(R, 'quintic-sextic')
print('D5 done', time.time() - t0, sorted((j, v[0]) for j, v in best.items()), flush=True)
for _ in range(3000000):
    c = [random.randint(-20, 20) for _ in range(6)]
    f = pari('x^6+' + '+'.join(f'({c[i]})*x^{i}' for i in range(6)))
    dsc = pari.poldisc(f)
    if dsc <= 0 or not dsc.issquare(): continue
    consider(f, 'dense-squaredisc')
    if 15 in best: break
print('E done', time.time() - t0, sorted((j, v[0]) for j, v in best.items()), flush=True)
json.dump({str(j): best[j] for j in sorted(best)}, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'search_F.json'), 'w'), indent=1)
for j in sorted(best): print(j, best[j])
