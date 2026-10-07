# fields_pari.py -- explicit CM fields (local Python + cypari/PARI).
# Builds explicit CM fields and the Galois action on their complex embeddings, in the cmenum labelling
# (points 0..11, rho(i) = i+6 = complex conjugation), for candidate OPEN non-simple CM sixfolds:
#   X1  (5,1):  E = E0(i), E0 = Q[x]/(x^5-5x^3-x^2+3x+1) totally real with Galois group S5; K = Q(i).
#               Gal = S5 x C2 (combinatorial model, justified by polgalois + disc > 0).
#   D4 family:  L = Q(sqrt(3+sqrt2), sqrt(3-sqrt2), i), Gal(L/Q) of order 16 (computed with nfgaloisconj).
#               All CM subfields of degree 8, 4, 2 of L; the action of Gal(L/Q) on their embeddings is computed
#               exactly (images of a generator under the 16 automorphisms) and complex conjugation is identified
#               numerically from a fixed complex embedding of L.
# Output: fields_out.json  (configs: name, description, gens (list of 12-perms), point labels).
import json, itertools, os
from cypari import pari

HERE = os.path.dirname(os.path.abspath(__file__))

pari.allocatemem(2 * 10**9)
out = dict(configs=[], facts={})

# ---------------- X1: S5 x C2 ----------------
f = pari('x^5-5*x^3-x^2+3*x+1')
gal = pari.polgalois(f)
out['facts']['X1_quintic'] = dict(poly=str(f), polgalois=str(gal), polsturm=int(pari.polsturm(f)),
                                  disc=int(pari.poldisc(f)), disc_factor=str(pari.factor(pari.poldisc(f))),
                                  nfdisc=int(pari.nfdisc(f)))
assert int(gal[0]) == 120 and int(pari.polsturm(f)) == 5 and int(pari.poldisc(f)) > 0
# labels: 0..4 = (sigma_k, +), 5 = (K,+), 6..10 = (sigma_k, -), 11 = (K,-); '+' = restriction i -> +i
def p_shift(sig):   # sig: permutation of 0..4 acting on sigma_k
    return [sig[x] if x < 5 else (5 if x == 5 else (6 + sig[x - 6] if x < 11 else 11)) for x in range(12)]
c5 = [1, 2, 3, 4, 0]
t2 = [1, 0, 2, 3, 4]
cc = [(x + 6) % 12 for x in range(12)]
lab = ['E(s%d,+i)' % k for k in range(5)] + ['K(+i)'] + ['E(s%d,-i)' % k for k in range(5)] + ['K(-i)']
out['configs'].append(dict(name='X1_S5xC2_(5,1)', desc='A_E (E=E0(i), E0 S5 quintic) x E_K, K=Q(i)',
                           gens=[p_shift(c5), p_shift(t2), cc], labels=lab,
                           note='types: Phi = {0,1,2,3,10} u {5 or 11}: K-signature (4,1) on A_E'))

# ---------------- D4 family ----------------
P = pari('x^4-6*x^2+7')
out['facts']['D4_quartic'] = dict(poly=str(P), polgalois=str(pari.polgalois(P)), polsturm=int(pari.polsturm(P)),
                                  disc_factor=str(pari.factor(pari.poldisc(P))))
L0 = pari.nfsplitting(P)
Lpol = pari.polredbest(pari.polcompositum(L0, pari('x^2+1'))[0])
nf = pari.nfinit(Lpol)
auts = pari.nfgaloisconj(nf)
nA = len(auts)
out['facts']['L'] = dict(poly=str(Lpol), degree=int(pari.poldegree(Lpol)), n_auts=nA,
                         polsturm=int(pari.polsturm(Lpol)))
assert nA == 16 and int(pari.poldegree(Lpol)) == 16
z = pari.polroots(Lpol)[0]
def ev(a):           # a: polynomial in x (element of L), evaluate at the chosen complex root z
    return complex(pari.subst(pari.lift(pari.Mod(a, Lpol)), 'x', z))
def app(s, a):       # apply automorphism s (polynomial) to element a
    return pari.lift(pari.Mod(pari.subst(pari.lift(pari.Mod(a, Lpol)), 'x', s), Lpol))
# complex conjugation
zc = complex(z).conjugate()
cidx = [k for k in range(nA) if abs(ev(auts[k]) - zc) < 1e-8 * (1 + abs(zc))]
assert len(cidx) == 1
sc = auts[cidx[0]]
# check central
for s in auts:
    assert app(s, app(sc, pari('x'))) == app(sc, app(s, pari('x')))

def cm_subfields(deg):
    res = []
    for pol, emb in pari.nfsubfields(nf, deg):
        if int(pari.polsturm(pol)) == 0:
            res.append((pari.polredabs(pol), emb))
    # isomorphism classes (polredabs is canonical)
    seen = {}
    for pol, emb in res:
        seen.setdefault(str(pol), (pol, emb))
    return list(seen.values())

def embeddings(emb):
    imgs = []
    for s in auts:
        e = app(s, emb)
        if all(e != i for i in imgs):
            imgs.append(e)
    return imgs

fields = {}
for deg in (2, 4, 8):
    for pol, emb in cm_subfields(deg):
        imgs = embeddings(emb)
        assert len(imgs) == deg
        fields[str(pol)] = dict(deg=deg, pol=pol, emb=emb, imgs=imgs,
                                galois=str(pari.polgalois(pol)) if deg <= 7 else None,
                                real_sub=str(pari.polredabs(pari.nfsubfields(pol, deg // 2)[0][0])) if deg > 2 else None)
out['facts']['cm_subfields'] = {k: dict(deg=v['deg'], polgalois=v['galois']) for k, v in fields.items()}

def build(field_keys, name, desc):
    # points: embeddings of each field, paired under complex conjugation
    pts = []
    for fk in field_keys:
        F = fields[fk]
        for t in F['imgs']:
            pts.append((fk, t))
    n = len(pts)
    assert n == 12
    def find(fk, t):
        for i, (gk, u) in enumerate(pts):
            if gk == fk and u == t:
                return i
        raise ValueError
    conj = [find(fk, app(sc, t)) for fk, t in pts]
    reps, seen = [], set()
    for i in range(n):
        if i not in seen:
            reps.append(i); seen.add(i); seen.add(conj[i])
    assert len(reps) == 6
    order = reps + [conj[i] for i in reps]          # new label j -> old index order[j]
    newlab = {old: j for j, old in enumerate(order)}
    gens = []
    for s in auts:
        perm = [newlab[find(pts[order[j]][0], app(s, pts[order[j]][1]))] for j in range(n)]
        gens.append(perm)
    assert gens[cidx[0]] == [(j + 6) % 12 for j in range(12)]
    labels = []
    for j in range(n):
        fk, t = pts[order[j]]
        v = ev(t)
        labels.append('%s|%s|%.4f%+.4fi' % (fields[fk]['deg'], fk, v.real, v.imag))
    out['configs'].append(dict(name=name, desc=desc, gens=gens, labels=labels, fields=field_keys))

octics = [k for k, v in fields.items() if v['deg'] == 8]
quartics = [k for k, v in fields.items() if v['deg'] == 4]
quads = [k for k, v in fields.items() if v['deg'] == 2]
for o in octics:
    for q in quartics:
        build([o, q], 'D4_(4,2)_%s_%s' % (octics.index(o), quartics.index(q)), 'A(octic %s) x S(quartic %s)' % (o, q))
    # (4,1,1): octic x two imaginary quadratic subfields of the octic
    subs_o = [str(pari.polredabs(p)) for p, e in pari.nfsubfields(fields[o]['pol'], 2)]
    qin = [k for k in quads if k in subs_o]
    for a, b in itertools.combinations(qin, 2):
        build([o, a, b], 'D4_(4,1,1)_%s_%s_%s' % (octics.index(o), quads.index(a), quads.index(b)),
              'A(octic %s) x E(%s) x E(%s)' % (o, a, b))
out['facts']['octics'] = octics; out['facts']['quartics'] = quartics; out['facts']['quads'] = quads
json.dump(out, open(os.path.join(HERE, 'fields_out.json'), 'w'), indent=1, default=str)
print(json.dumps(out['facts'], indent=1, default=str))
print('configs:', [c['name'] for c in out['configs']])
