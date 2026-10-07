"""INDEPENDENT exact certificate that E = Q(pi) = F.K with K = Q(sqrt(-d)) imaginary quadratic (no PARI).
F = Q(y), y = pi + q/pi, h(y) = real Weil polynomial; E = F(sqrt(delta)), delta = y^2 - 4q.
Q(sqrt(-d)) is contained in E  <=>  -d*delta is a square in F. Search d in a small list, numerically interpolate
a square root c(y) (degree <= 5) over the 6 real embeddings and all sign patterns, round disc(h)*c to integers,
then check EXACTLY (FLINT fmpq_poly) that c(y)^2 + d*(y^2 - 4q) == 0 mod h(y).
Consequence (with GaloisType(h) = 6Tj from verify_gap.py): Gal(E) = 6Tj x C2 in the product action (F~ is totally
real so F~ and K are disjoint); its 12T label is read off combi_rows.out (next to this script).
Usage: python verify_K.py <examples.json> <verify_gap output.json> <output.json>"""
import json, sys, itertools, os
import mpmath, flint
from verify import real_weil_poly
HERE = os.path.dirname(os.path.abspath(__file__))
ROWS = {}
for line in open(os.path.join(HERE, 'combi_rows.out')):
    parts = dict(kv.split('=', 1) for kv in line.split() if '=' in kv and kv.split('=', 1)[0] in ('j', 'TI'))
    if 'j' in parts: ROWS[int(parts['j'])] = parts['TI']
GAP = {r['name']: r for r in json.load(open(sys.argv[2]))}
out = []
for ex in json.load(open(sys.argv[1])):
    a = [int(c) for c in ex['weil_poly']]; q = int(ex['q'])
    h = real_weil_poly(a, q)
    H = flint.fmpz_poly(h); Dh = abs(int(H.discriminant()))
    mpmath.mp.dps = 3000
    ys = sorted(mpmath.polyroots([mpmath.mpf(c) for c in reversed(h)], maxsteps=2000, extraprec=20000))
    assert all(abs(mpmath.im(y)) < mpmath.mpf(10) ** -2000 for y in ys)
    ys = [mpmath.re(y) for y in ys]
    V = mpmath.matrix([[y ** k for k in range(6)] for y in ys])
    found = None
    for d in (1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 19, 21, 22, 23, 31, 35, 39, 43):
        s = [mpmath.sqrt(d * (4 * q - y ** 2)) for y in ys]
        for sg in itertools.product((1, -1), repeat=5):
            rhs = mpmath.matrix([s[0]] + [e * x for e, x in zip(sg, s[1:])])
            c = mpmath.lu_solve(V, rhs)
            cz = [Dh * c[k] for k in range(6)]
            if max(abs(z - mpmath.nint(z)) for z in cz) > mpmath.mpf(10) ** -100: continue
            C = flint.fmpq_poly([int(mpmath.nint(z)) for z in cz], Dh)
            Hq = flint.fmpq_poly(h)
            rem = (C * C + flint.fmpq_poly([-4 * q * d, 0, d])) % Hq
            if rem == 0:
                found = (d, [str(x) for x in C.coeffs()]); break
        if found: break
    gt = GAP.get(ex['name'], {})
    jj = int(gt['GaloisType_h'][2:]) if gt else None
    lab = ROWS.get(jj) if (found and jj) else None
    rec = dict(name=ex['name'], d=found[0] if found else None, sqrt_cert=found[1] if found else None,
               GaloisType_h=gt.get('GaloisType_h'), label_12T=lab, claimed_row=ex.get('row'),
               label_matches=(lab == ex.get('row')))
    out.append(rec)
    print(json.dumps({k: v for k, v in rec.items() if k != 'sqrt_cert'}), flush=True)
json.dump(out, open(sys.argv[3], 'w'), indent=1)
