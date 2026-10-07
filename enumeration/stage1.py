# stage1.py -- transitive case g = 6 (CM fields of degree 12), route (i): TransitiveGroups(12) with a central
# fixed-point-free involution rho'.
# Conjugate each (G, rho') into W(B6) (rho' -> standard rho), analyse all CM types, and match with
# route (ii) = transitive classes in enum_g6.json (ConjugacyClassesSubgroups of W(B6)).
import sys, os, json, collections, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sage.all import libgap
from cmenum import *
from locate import Stored

t0 = time.time()
out = open(os.path.join(HERE, 'stage1_out.txt'), 'w')
def P(*a):
    print(*a); print(*a, file=out); out.flush()

g = 6
n = 12
S12 = libgap.SymmetricGroup(n)
rho0 = rho_perm(g)
st = Stored(g)
W = st.W
N = int(libgap.NrTransitiveGroups(n))
P("NrTransitiveGroups(12) = %d" % N)
pairs = []          # (k, rho' index, normalizer-class rep?)
groups_with = 0
matched = {}        # stored gid -> list of (k, j)
rows = []
for k in range(1, N + 1):
    T = libgap.TransitiveGroup(n, k)
    Z = libgap.Centre(T)
    invs = [z for z in libgap.Elements(Z) if int(libgap.Order(z)) == 2 and int(libgap.NrMovedPoints(z)) == n]
    if not invs:
        continue
    groups_with += 1
    NT = libgap.Normalizer(S12, T)
    # classes of central fpf involutions under N_{S12}(T)
    cls = []
    for z in invs:
        if not any(bool(libgap.IsConjugate(NT, z, y)) for y in cls):
            cls.append(z)
    for j, z in enumerate(cls):
        pi = libgap.RepresentativeAction(S12, z, rho0)
        G = libgap.ConjugateGroup(T, pi)
        assert bool(libgap.IsSubset(G, [rho0]))
        gd, res = analyse_group(G, g, W, True)
        ti = gd['oinfo'][0]['TI']
        assert ti == "12T%d" % k
        # locate this group among stored transitive reps
        hit = None
        for gid, Gd in st.groups.items():
            if not Gd['transitive'] or Gd['order'] != gd['order'] or Gd['orbit_TI'][0] != ti:
                continue
            if libgap.RepresentativeAction(W, G, st.gapgroup(gid)) != libgap.eval('fail'):
                assert hit is None
                hit = gid
        matched.setdefault(hit, []).append((k, j))
        nprim = sum(1 for ta in res if ta['primitive'][0])
        ndeg = sum(1 for ta in res if ta['primitive'][0] and ta['d'] >= 1)
        rows.append(dict(k=k, j=j, n_invs=len(invs), n_inv_classes=len(cls), order=gd['order'], gid=hit,
                         ntypes=len(res), nNtypes=len(set(ta['Nclass'] for ta in res)), nprim=nprim, ndeg=ndeg,
                         prim_N=len(set(ta['Nclass'] for ta in res if ta['primitive'][0])),
                         deg=[dict(Phi=ta['Phi'], Nclass=ta['Nclass'], d=ta['d'], dimMT=ta['dimMT'], LU=ta['LU'],
                                   adm=[(a['K'], a['B']) for a in ta['adm']], passed=ta['passed'],
                                   Ks=[(kk, K['blocks'][0], 2 * len(set(ta['Phi']) & set(K['blocks'][0])) - len(K['blocks'][0]))
                                       for kk, K in enumerate(gd['Ks'])])
                              for ta in res if ta['primitive'][0] and ta['d'] >= 1]))
P("12T groups with a central fpf involution: %d ; (G, rho') pairs up to N_S12(G): %d ; time %.1fs" %
  (groups_with, len(rows), time.time() - t0))
trans_gids = sorted(gid for gid, Gd in st.groups.items() if Gd['transitive'])
P("route (ii): transitive W(B6)-classes containing rho in enum_g6: %d" % len(trans_gids))
P("matching: stored gids hit = %d ; unmatched (i) = %d ; stored gids hit more than once = %d ; stored gids never hit = %s" %
  (len([x for x in matched if x is not None]), len(matched.get(None, [])),
   sum(1 for x, v in matched.items() if x is not None and len(v) > 1), sorted(set(trans_gids) - set(matched))))
lab_i = sorted(set(r['k'] for r in rows))
lab_ii = sorted(set(int(st.groups[gid]['orbit_TI'][0].split('T')[1]) for gid in trans_gids))
P("TransitiveIdentification label sets equal: %s (|(i)|=%d, |(ii)|=%d)" % (lab_i == lab_ii, len(lab_i), len(lab_ii)))
P("labels: %s" % lab_i)
tot = collections.Counter()
for r in rows:
    tot['types(G-orbits)'] += r['ntypes']; tot['types(N-orbits)'] += r['nNtypes']
    tot['primitive(G-orbits)'] += r['nprim']; tot['primitive(N-orbits)'] += r['prim_N']
    tot['prim.degenerate(G-orbits)'] += r['ndeg']; tot['prim.degenerate(N-orbits)'] += len(set(x['Nclass'] for x in r['deg']))
P("totals over (G,rho') pairs: %s" % dict(tot))
P("table of degenerate primitive types (one line per N_W(G)-class):")
P("  %-8s %-5s %-5s %-24s %-3s %-6s %-18s %s" % ("label", "|G|", "gid", "Phi", "d", "test", "LU", "balanced K: block (|Phi cap B|=3) ; other K: t=2|Phi cap B|-6"))
dcount = collections.Counter(); main_ok = True
for r in rows:
    seen = set()
    for x in r['deg']:
        if x['Nclass'] in seen:
            continue
        seen.add(x['Nclass'])
        bal = [(kk, B) for kk, B, t in x['Ks'] if t == 0]
        oth = [(kk, t) for kk, B, t in x['Ks'] if t != 0]
        uniqueK = (len(bal) == 1)
        spanned = (x['d'] == 1 and uniqueK and x['passed'])
        main_ok &= spanned
        dcount[x['d']] += 1
        P("  12T%-5d %-5d %-5s %-24s %-3d %-6s %-18s balanced=%s others=%s" %
          (r['k'], r['order'], r['gid'], x['Phi'], x['d'], x['passed'], x['LU'], bal, oth))
P("degenerate primitive N-classes by d: %s" % dict(dcount))
P("MAIN QUESTION (every primitive degenerate degree-12 type has d=1, Lambda_U = Z w_K, unique balanced K): %s" % main_ok)
json.dump(rows, open(os.path.join(HERE, 'stage1_rows.json'), 'w'), default=int)
P("time %.1fs" % (time.time() - t0))
