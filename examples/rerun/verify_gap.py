"""INDEPENDENT Galois-group check. GaloisType(P) (degree 12) is run only for entries with try12=true
(it stalled >8 min on a 12T7 polynomial with large coefficients); GaloisType(h) (degree 6) always.
 (GAP library GaloisType, Hulpke's resolvent method; no PARI).
Run in sagemath docker: sage -python verify_gap.py <in.json> <out.json>. For each Weil polynomial P (deg 12)
and its real Weil polynomial h (deg 6): GaloisType -> transitive group number, and the order of TransitiveGroup(6,k)
(GaloisType(P) only for entries with try12=true). The 12T label 6Tj x C2 is assigned in verify_K.py
(via combi_rows.out), not here."""
import sys, json, time
from sage.all import libgap
def real_weil(a, q):
    s = [[2], [0, 1]]
    for k in range(1, 6):
        nxt = [0] + s[k]; prev = s[k - 1] + [0] * (len(nxt) - len(s[k - 1]))
        s.append([nxt[i] - q * prev[i] for i in range(len(nxt))])
    h = [0] * 7; h[0] += a[6]
    for k in range(1, 7):
        for i, c in enumerate(s[k]): h[i] += a[6 - k] * c
    return h
libgap.eval('x := Indeterminate(Rationals, "x");')
exs = json.load(open(sys.argv[1])); out = []
for ex in exs:
    a = [int(c) for c in ex['weil_poly']]; q = int(ex['q'])
    Ps = '+'.join(f'({c})*x^{12-i}' for i, c in enumerate(a))
    h = real_weil(a, q)
    hs = '+'.join(f'({c})*x^{i}' for i, c in enumerate(h))
    t0 = time.time()
    k12 = None
    if ex.get('try12', False):
        k12 = int(libgap.eval(f'GaloisType({Ps})'))
    k6 = int(libgap.eval(f'GaloisType({hs})'))
    t1 = time.time()
    order6 = int(libgap.eval(f'Size(TransitiveGroup(6,{k6}))'))
    rec = dict(name=ex['name'], GaloisType_P=(f'12T{k12}' if k12 else 'not run'), GaloisType_h=f'6T{k6}',
               order_h=order6, sec=round(t1 - t0, 1))
    print(json.dumps(rec), flush=True); out.append(rec)
json.dump(out, open(sys.argv[2], 'w'), indent=1)
