# controls.py -- positive controls at g = 4 and g = 6 (main computation, Sage+GAP).
#  (1) Milne 2112.12815 Ex. 1.12 with n=3 at g=4 (E = K.F sextic, signature (1,2), times E_K)
#  (2) APFV exceptional rows g=4 and g=6, located in our enumeration (needs apfv_data/, see fetch_apfv_data.py)
#  (3) E = Q(zeta_21), G=(Z/21)^x, Phi={1,2,4,5,8,10}
import sys, os, re, json
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sage.all import libgap
from cmenum import *
from locate import Stored

out = open(os.path.join(HERE, 'controls_out.txt'), 'w')
def P(*a):
    print(*a); print(*a, file=out); out.flush()

def show(tag, gd, ta):
    P("  %s: |G|=%d orbits=%s TI=%s Phi=%s dimMT=%d d=%d LU=%s primitive=%s reduced=%s" %
      (tag, gd['order'], [len(o) for o in gd['orbits']], [o['TI'] for o in gd['oinfo']], ta['Phi'], ta['dimMT'],
       ta['d'], ta['LU'], ta['primitive'], ta['reduced']))
    P("     admissible (K,J,B,|B|): %s ; rankW=%d index=%s KEYTEST=%s" %
      ([(a['K'], a['J'], a['B'], a['size']) for a in ta['adm']], ta['rankW'], ta['index'], ta['passed']))

def analyse_one(G, g, Phi):
    gd = group_data(G, g, None, False)
    for ot in gd['torbs']:
        if tuple(sorted(Phi)) in ot:
            return gd, type_analysis(gd, ot)

# ---------------- (1) Milne Example 1.12, n = 3 ----------------
P("=== (1) Milne 2112.12815 Ex. 1.12 (n=3): A (CM by sextic E=K.F, type {phi0, i.phi1, i.phi2}) x B (CM by K)")
# points: (+,1)(+,2)(+,3)(+K) = 0,1,2,3 ; (-,1)(-,2)(-,3)(-K) = 4,5,6,7 ; rho = sign flip
g = 4
rho = rho_perm(g)
S3gens = [libgap.eval("(1,2)(5,6)"), libgap.eval("(1,2,3)(5,6,7)")]   # S3 permuting j (1-based points)
C3gens = [libgap.eval("(1,2,3)(5,6,7)")]
PhiM = [0, 5, 6, 3]   # phi0=(+,1), iota.phi1=(-,2), iota.phi2=(-,3), and B of type (K, sigma0)=(+K)
st4 = Stored(4)
for name, gens in [("Gal = C2 x S3 (F non-Galois cubic)", S3gens), ("Gal = C2 x C3 (E cyclic sextic)", C3gens)]:
    G = libgap.Group(gens + [rho])
    gd, ta = analyse_one(G, g, PhiM)
    show(name, gd, ta)
    gid, case = st4.locate(G, PhiM)
    P("     located in enum_g4: gid=%s stored Phi=%s d=%s passed=%s reduced=%s" %
      (gid, case['Phi'], case['d'], case['passed'], case['reduced']))

# ---------------- (2) APFV rows ----------------
def parse_apfv(fn):
    t = open(fn).read()
    rows = []
    for m in re.finditer(r'<"([^"]+)",\s*\[(.*?)\]\s*,\s*(\d+),\s*(true|false)>', t, re.S):
        lab, gens, ar, ex = m.group(1), m.group(2), int(m.group(3)), m.group(4)
        L = [list(map(int, x.split(','))) for x in re.findall(r'\[([^\[\]]+)\]', gens)]
        rows.append((lab, L, ar, ex))
    return rows

def signed_to_perm(L, g):
    idx = lambda s: s - 1 if s > 0 else g + (-s) - 1
    p = [None] * (2 * g)
    for i in range(1, g + 1):
        p[idx(i)] = idx(L[i - 1])
        p[idx(-i)] = idx(L[g + i - 1])
        assert L[g + i - 1] == -L[i - 1]
    return p

apfv_dir = os.path.join(HERE, 'apfv_data', '')
apfv_res = {}
for g, fn in [(4, '4-0.000_0.000_0.000_0.000_1.00_1.00_1.00_1.00.m'),
              (6, '6-0.000_0.000_0.000_0.000_0.000_0.000_1.00_1.00_1.00_1.00_1.00_1.00.m')]:
    P("=== (2) APFV exceptional rows, g=%d (Phi = slope-0 set = {1..g} -> points 0..g-1)" % g)
    st = Stored(g)
    rho = rho_perm(g)
    for lab, L, ar, ex in parse_apfv(apfv_dir + fn):
        G = libgap.Group([libgap.PermList([x + 1 for x in signed_to_perm(l, g)]) for l in L])
        has_rho = bool(libgap.IsSubset(G, [rho]))
        Phi = list(range(g))
        gd, ta = analyse_one(G, g, Phi)
        gid, case = st.locate(G, Phi)
        tl = lab.split('.')[0]
        ok = (gd['oinfo'][0]['TI'] == tl and ta['dimMT'] - 1 == ar and case is not None and case['Phi'] is not None)
        P("  %-18s APFV angle rank %d | ours: rho in G=%s TI=%s |G|=%d dimMT=%d (=angle rank+1? %s) d=%d primitive=%s "
          "adm=%s KEYTEST=%s | located gid=%s storedPhi=%s Nclass=%s" %
          (lab, ar, has_rho, gd['oinfo'][0]['TI'], gd['order'], ta['dimMT'], ta['dimMT'] - 1 == ar, ta['d'],
           ta['primitive'], [(a['K'], a['B']) for a in ta['adm']], ta['passed'], gid, case['Phi'] if case else None,
           case['Nclass'] if case else None))
        apfv_res.setdefault(g, []).append(dict(label=lab, gid=gid, Phi=case['Phi'] if case else None, Nclass=case['Nclass'] if case else None,
                                               TI=gd['oinfo'][0]['TI'], d=ta['d'], passed=ta['passed'], ok=ok))
    # completeness: every reduced transitive d>=1 class in our list is an APFV row?
    ours = set((c['gid'], c['Nclass']) for c in st.D['cases'] if c['reduced'] and c['d'] >= 1 and st.groups[c['gid']]['transitive'])
    theirs = set((r['gid'], r['Nclass']) for r in apfv_res[g])
    P("  our transitive reduced degenerate (G,Phi) up to N_W(G): %d ; APFV rows: %d ; APFV rows found among ours: %d ; ours not in APFV: %s"
      % (len(ours), len(apfv_res[g]), len(theirs & ours), sorted(ours - theirs)))

# ---------------- (3) Q(zeta_21) ----------------
P("=== (3) E = Q(zeta_21), G = (Z/21)^x regular, rho = -1, Phi = {1,2,4,5,8,10}")
g = 6
units = [a for a in range(1, 21) if gcd(a, 21) == 1]
reps = [1, 2, 4, 5, 8, 10]
lab = {}
for i, a in enumerate(reps):
    lab[a] = i
    lab[21 - a] = i + 6
gens = []
for m in [2, 20, 4]:
    gens.append(libgap.PermList([lab[(m * a) % 21] + 1 for a in sorted(lab, key=lambda a: lab[a])]))
G = libgap.Group(gens)
assert int(libgap.Size(G)) == 12
Phi = [lab[a] for a in [1, 2, 4, 5, 8, 10]]
gd, ta = analyse_one(G, g, Phi)
show("Q(zeta21)", gd, ta)
inv = {v: k for k, v in lab.items()}
for kidx, K in enumerate(gd['Ks']):
    blk = sorted(inv[x] for x in K['blocks'][0])
    Bsub = set(blk)
    name = '?'
    if Bsub == set(a for a in units if a % 3 == 1):
        name = 'Q(sqrt-3)  [block = {a = 1 mod 3}]'
    elif Bsub == set(a for a in units if a % 3 == 2):
        name = 'Q(sqrt-3)  [block = {a = 2 mod 3}]'
    elif Bsub == set(a for a in units if a % 7 in (1, 2, 4)):
        name = 'Q(sqrt-7)  [block = {a mod 7 in squares}]'
    elif Bsub == set(a for a in units if a % 7 in (3, 5, 6)):
        name = 'Q(sqrt-7)  [block = {a mod 7 non-squares}]'
    t = 2 * len(set(Phi) & set(K['blocks'][0])) - len(K['blocks'][0])
    P("   K%d: block (as residues mod 21) %s -> %s ; |Phi cap B| = %d (balanced: %s)" %
      (kidx, blk, name, len(set(Phi) & set(K['blocks'][0])), t == 0))
st6 = Stored(6)
gid, case = st6.locate(G, Phi)
P("   located in enum_g6: gid=%s stored Phi=%s d=%s passed=%s reduced=%s TI=%s" %
  (gid, case['Phi'], case['d'], case['passed'], case['reduced'], st6.groups[gid]['orbit_TI']))
json.dump(apfv_res, open(os.path.join(HERE, 'controls_apfv.json'), 'w'), default=int)
