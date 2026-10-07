# final_tables.py -- local python: compile final tables from the enum_g*.json files (run_enum.py) and
# general_weil2_out_D12.json (general_weil2.py).
import json, collections, sys, os
from verify_enum import group_struct

HERE = os.path.dirname(os.path.abspath(__file__))
gw = json.load(open(os.path.join(HERE, 'general_weil2_out_D12.json')))
out = open(os.path.join(HERE, 'final_tables.txt'), 'w')
def P(*a):
    print(*a, file=out)


def gmaps_exist(gens, Oa, Ob):
    """is there a G-map Oa -> Ob (i.e. is the field of Ob a subfield of the field of Oa)?"""
    a0 = Oa[0]
    for b0 in Ob:
        f = {a0: b0}; st = [a0]; ok = True
        while st and ok:
            x = st.pop()
            for p in gens:
                y, fy = p[x], p[f[x]]
                if y in f:
                    if f[y] != fy:
                        ok = False; break
                else:
                    f[y] = fy; st.append(y)
        if ok:
            return True
    return False


for g in range(1, 7):
    D = json.load(open(os.path.join(HERE, 'enum_g%d.json' % g)))
    groups = {G['gid']: G for G in D['groups']}
    cases = D['cases']
    gwr = {(r['gid'], tuple(r['Phi'])): r for r in gw.get(str(g), [])}
    P("=" * 110)
    red = [c for c in cases if c['reduced']]
    P("g=%d  W-classes of G>=<rho>: %d (transitive %d) | (G,Phi) G-orbit cases: %d (N_W(G)-classes %d) | reduced: %d (N-classes %d)" %
      (g, len(groups), sum(G['transitive'] for G in groups.values()), len(cases),
       len(set((c['gid'], c['Nclass']) for c in cases)), len(red), len(set((c['gid'], c['Nclass']) for c in red))))
    for lab, sel in [('ALL', cases), ('REDUCED', red)]:
        P("  %-8s d-dist %s ; strict KEY TEST pass %d / fail %d" % (lab, dict(sorted(collections.Counter(c['d'] for c in sel).items())),
                                                                  sum(c['passed'] for c in sel), sum(not c['passed'] for c in sel)))
    nr = [c for c in cases if not c['reduced']]
    nrf = [c for c in nr if not c['passed']]
    P("  non-reduced: %d cases, dimMT(reduced model)=dimMT ok for %d ; strict failures %d, all of whose reduced models pass: %s" %
      (len(nr), sum(c['red_dimMT_check'] for c in nr), len(nrf), all(c['red_model']['passed'] and c['red_model']['reduced'] for c in nrf)))
    # reduced d>=1 by partition / tier
    tab = collections.defaultdict(collections.Counter)
    for c in red:
        if c['d'] < 1:
            continue
        G = groups[c['gid']]
        part = tuple(sorted((len(o) // 2 for o in G['orbits']), reverse=True))
        r = gwr[(c['gid'], tuple(c['Phi']))]
        tier = 'strict' if c['passed'] else ('IQmult<=6' if c['mult']['Dmin'] is not None and c['mult']['Dmin'] <= 6 else
                                             ('needs D=%s' % r['Dmin']))
        tab[(part, c['d'], tier)]['G-orbit cases'] += 1
        tab[(part, c['d'], tier)].setdefault('N', set())
    if tab:
        P("  reduced d>=1 by (factor dims, d, explanation tier):")
        for k in sorted(tab, key=lambda k: (-k[0][0], k)):
            ncl = len(set((c['gid'], c['Nclass']) for c in red if c['d'] >= 1 and
                          (tuple(sorted((len(o) // 2 for o in groups[c['gid']]['orbits']), reverse=True)), c['d'],
                           'strict' if c['passed'] else ('IQmult<=6' if c['mult']['Dmin'] is not None and c['mult']['Dmin'] <= 6 else
                                                         'needs D=%s' % gwr[(c['gid'], tuple(c['Phi']))]['Dmin'])) == k))
            P("    %-14s d=%d  %-12s G-orbit cases %3d   N-classes %3d" % (k[0], k[1], k[2], tab[k]['G-orbit cases'], ncl))
    # full list of reduced strict failures, one line per N-class + data
    fails = [c for c in red if not c['passed']]
    if fails:
        P("  REDUCED STRICT FAILURES (one entry per N_W(G)-class; %d G-orbit cases in %d N-classes):" %
          (len(fails), len(set((c['gid'], c['Nclass']) for c in fails))))
        seen = set()
        for c in fails:
            key = (c['gid'], c['Nclass'])
            if key in seen:
                continue
            seen.add(key)
            G = groups[c['gid']]
            gs = group_struct(G, g)
            r = gwr[(c['gid'], tuple(c['Phi']))]
            Ps = set(c['Phi'])
            kinfo = []
            for K in gs['Ks']:
                kinfo.append({str(o): 2 * len(Ps & set(B)) - len(B) for o, B in sorted(K.items())})
            quot = [(i, j) for i in range(len(G['orbits'])) for j in range(len(G['orbits']))
                    if i != j and len(G['orbits'][i]) > len(G['orbits'][j]) and gmaps_exist(G['gens'], G['orbits'][i], G['orbits'][j])]
            P("    gid=%d |G|=%d orbits(sizes)=%s TI=%s Phi=%s dimMT=%d d=%d" %
              (c['gid'], G['order'], [len(o) for o in G['orbits']], G['orbit_TI'], c['Phi'], c['dimMT'], c['d']))
            P("       Lambda_U (HNF, coords = points 0..g-1, [i+g] = -[i]) = %s ; strict admissible = %s ; rank Lambda_W(strict) = %d" %
              (c['LU'], [(a['J'], a['B']) for a in c['adm']], c['rankW']))
            P("       imag.quadratic K's: signature excess t_i=2|Phi cap B_i|-|B_i| per orbit: %s ; G-maps O_i->O_j (E_j subfield of E_i): %s" %
              (kinfo, quot))
            P("       IQ-multiplicity Dmin=%s (qrank %d) ; general CM-field Dmin=%s, fields used (index [G:M]=[K':Q], D)=%s" %
              (c['mult']['Dmin'], c['mult']['qrank'], r['Dmin'], r['fields_used']))
            P("       orbits=%s" % G['orbits'])
            P("       gens=%s" % G['gens'])
print(open(os.path.join(HERE, 'final_tables.txt')).read() if False else "written final_tables.txt")
