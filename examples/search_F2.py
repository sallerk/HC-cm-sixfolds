"""Extra search: totally real sextics with Galois group 6T12 (A5) and 6T14 (S5) as F20-resolvent sextics of
totally real quintics (A5: square discriminant). Keeps up to 8 smallest-disc fields per group. -> extra_F2.json"""
import json, random, itertools, os
from cypari import pari
pari('default(new_galois_format,1)')
random.seed(7)
x = pari('x')
cyc = []
for perm in itertools.permutations(range(1, 5)):
    c5 = (0,) + perm
    key = min(min(tuple(c5[(s + k) % 5] for k in range(5)), tuple(c5[(s - k) % 5] for k in range(5))) for s in range(5))
    if key not in cyc: cyc.append(key)
def th(r, c):
    adj = sum(r[c[k]] * r[c[(k + 1) % 5]] for k in range(5))
    nonadj = sum(r[c[k]] * r[c[(k + 2) % 5]] for k in range(5))
    return (adj - nonadj) ** 2
found = {12: {}, 14: {}}
tries = 0
while (len(found[12]) < 8 or len(found[14]) < 8) and tries < 400000:
    tries += 1
    c = [random.randint(-25, 25) for _ in range(4)]
    g = pari(f'x^5+({c[3]})*x^3+({c[2]})*x^2+({c[1]})*x+({c[0]})')
    d = pari.poldisc(g)
    sq = d > 0 and d.issquare()
    if not sq and len(found[14]) >= 8: continue
    if not pari.polisirreducible(g) or int(pari.polsturm(g)) != 5: continue
    r = pari.polroots(g, precision=800)
    vals = []
    for cc in cyc:
        v = th(r, cc)
        if not any(abs(complex(v - u)) < 1e-60 * (1 + abs(complex(u))) for u in vals): vals.append(v)
    if len(vals) != 6: continue
    Pp = pari(1)
    for v in vals: Pp = Pp * (x - v)
    co = [Pp.polcoef(k) for k in range(7)]
    rc = [cc.real().round() for cc in co]
    if max(abs(float(cc.real() - rr)) + abs(float(cc.imag())) for cc, rr in zip(co, rc)) > 1e-30 * (1 + max(abs(float(rr)) for rr in rc)):
        continue
    R = sum(rc[k] * x ** k for k in range(7))
    if not pari.polisirreducible(R) or int(pari.polsturm(R)) != 6: continue
    R = pari.polredabs(R)
    j = int(pari.polgalois(R)[2])
    if j not in found: continue
    D = abs(int(pari.nfdisc(R)))
    found[j][str(R)] = D
json.dump({j: sorted(v.items(), key=lambda t: t[1])[:8] for j, v in found.items()}, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'extra_F2.json'), 'w'), indent=1)
for j in found: print(j, sorted(found[j].items(), key=lambda t: t[1])[:8])
print('tries', tries)
