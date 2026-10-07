"""Combinatorics (Sage/libgap, run in docker): which degree-12 CM Galois data arise as E = F.K
(F totally real sextic with Gal 6Tj, K imaginary quadratic), and which CM types on E are degenerate.
Points: (a,b) -> a + 6b, a in 0..5 (embedding of F), b in {0,1} (embedding of K); rho(i) = i+6 mod 12
(same convention as cmenum.py in the enumeration folder of this repository).
Balanced type for S (|S|=3): Phi_S = S u {a+6 : a not in S}.
For each j: TransitiveIdentification of G12 = 6Tj x C2 (product action); for every one of the 64 CM types:
primitivity (right-stabiliser test), rank of M[g][i] = 1_Phi(g(i)) (= angle rank + 1), and for degenerate
primitive ones: GAP conjugacy test of (G12, Phi) against the APFV row with the same 12T label.
Input: apfv_data/6-0.000_..._1.00.m (APFV data file, g = 6, ordinary slopes).
Output: combi_rows.json (written next to this script) and the summary lines on stdout (stored as combi_rows.out)."""
import json, re, itertools, os
HERE = os.path.dirname(os.path.abspath(__file__))
from sage.all import libgap, matrix, QQ

def parse_apfv(fn):
    t = open(fn).read(); rows = []
    for m in re.finditer(r'<"([^"]+)",\s*\[(.*?)\]\s*,\s*(\d+),\s*(true|false)>', t, re.S):
        L = [list(map(int, x.split(','))) for x in re.findall(r'\[([^\[\]]+)\]', m.group(2))]
        rows.append((m.group(1), L, int(m.group(3)), m.group(4)))
    return rows

def signed_to_perm(L, g=6):
    idx = lambda s: s - 1 if s > 0 else g + (-s) - 1
    p = [None] * (2 * g)
    for i in range(1, g + 1):
        p[idx(i)] = idx(L[i - 1]); p[idx(-i)] = idx(L[g + i - 1])
    return p

APFV = {}
for lab, L, ar, ex in parse_apfv(os.path.join(HERE, 'apfv_data', '6-0.000_0.000_0.000_0.000_0.000_0.000_1.00_1.00_1.00_1.00_1.00_1.00.m')):
    G = libgap.Group([libgap.PermList([x + 1 for x in signed_to_perm(l)]) for l in L])
    APFV[lab.split('.')[0]] = (lab, G, ar, ex)

S12 = libgap.SymmetricGroup(12)
def apfv_match(G12, Phi, tlabel):
    if tlabel not in APFV: return None
    lab, GA, ar, ex = APFV[tlabel]
    x = libgap.RepresentativeAction(S12, G12, GA)
    if x == libgap.eval('fail'): return (lab, 'group not conjugate')
    Phi1 = libgap.OnSets(libgap.Set([i + 1 for i in Phi]), x)
    N = libgap.Normalizer(S12, GA)
    for tgt, nm in [(list(range(1, 7)), 'slope0'), (list(range(7, 13)), 'slope1')]:
        y = libgap.RepresentativeAction(N, Phi1, libgap.Set(tgt), libgap.OnSets)
        if y != libgap.eval('fail'):
            return (lab, 'MATCH(' + nm + ')')
    return (lab, 'no match')

out = {}
lines = []
for j in range(1, 17):
    G6 = libgap.TransitiveGroup(6, j)
    gens12 = []
    for g in libgap.GeneratorsOfGroup(G6):
        im = [int(x) - 1 for x in libgap.ListPerm(g, 6)]
        gens12.append([im[a % 6] + 6 * (a // 6) for a in range(12)])
    gens12.append([(a + 6) % 12 for a in range(12)])
    G12 = libgap.Group([libgap.PermList([x + 1 for x in p]) for p in gens12])
    TI = 'NA'
    if bool(libgap.IsTransitive(G12, list(range(1, 13)))):
        TI = '12T%d' % int(libgap.TransitiveIdentification(G12))
    elts = [[int(x) - 1 for x in libgap.ListPerm(e, 12)] for e in libgap.AsList(G12)]
    elts = [e + list(range(len(e), 12)) for e in elts]
    res = []
    for m in range(64):
        Phi = sorted(a + 6 * ((m >> a) & 1) for a in range(6))
        Ph = set(Phi)
        S = [a for a in range(6) if a in Ph]
        Phit = [e for e in elts if e[0] in Ph]
        Y = [y for y in range(12) if all(e[y] in Ph for e in Phit)]
        prim = (Y == [0])
        rk = matrix(QQ, [[1 if e[i] in Ph else 0 for i in range(12)] for e in elts]).rank()
        rec = dict(Phi=Phi, S=S, balanced=(len(S) == 3), primitive=prim, rank=int(rk), angle_rank=int(rk) - 1)
        if prim and rk < 7:
            rec['apfv'] = apfv_match(G12, Phi, TI)
        res.append(rec)
    prim_bal = [r for r in res if r['primitive'] and r['balanced']]
    prim_deg = [r for r in res if r['primitive'] and r['rank'] < 7]
    summ = dict(j=j, G6=str(libgap.StructureDescription(G6)), order12=len(elts), TI=TI,
                n_prim=sum(r['primitive'] for r in res), n_prim_balanced=len(prim_bal),
                n_prim_degenerate=len(prim_deg),
                prim_degenerate_all_balanced=all(r['balanced'] for r in prim_deg),
                prim_balanced_ranks=sorted(set(r['rank'] for r in prim_bal)),
                apfv=sorted(set(str(r.get('apfv')) for r in prim_deg)))
    out[j] = dict(summary=summ, types=res, gens12=gens12)
    line = ' '.join('%s=%s' % (k, v) for k, v in summ.items())
    print(line, flush=True); lines.append(line)
json.dump(out, open(os.path.join(HERE, 'combi_rows.json'), 'w'), indent=0)
