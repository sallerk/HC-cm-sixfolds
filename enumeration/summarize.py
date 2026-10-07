# summarize.py -- reads enum_g{g}.json (output of run_enum.py) and prints statistics / tables.  Local python.
#   usage:  python summarize.py 1,2,3,4,5,6 summary_strict.txt
import json, sys, os, collections

HERE = os.path.dirname(os.path.abspath(__file__))

def summarize(g, out):
    D = json.load(open(os.path.join(HERE, 'enum_g%d.json' % g)))
    groups = {G['gid']: G for G in D['groups']}
    cases = D['cases']
    P = lambda *a: print(*a, file=out)
    P("=" * 100)
    P("g = %d : groups (W-classes, rho in G) = %d ; transitive groups = %d ; (G,Phi) cases (G-orbits of types) = %d"
      % (g, len(groups), sum(G['transitive'] for G in groups.values()), len(cases)))
    ncls = len(set((c['gid'], c['Nclass']) for c in cases))
    P("   (G,Phi) up to N_W(G) (i.e. up to W-conjugacy of pairs): %d" % ncls)
    for label, sel in [("ALL", lambda c: True), ("REDUCED", lambda c: c['reduced']),
                       ("REDUCED & transitive", lambda c: c['reduced'] and groups[c['gid']]['transitive']),
                       ("REDUCED & non-simple", lambda c: c['reduced'] and not groups[c['gid']]['transitive'])]:
        cs = [c for c in cases if sel(c)]
        dd = collections.Counter(c['d'] for c in cs)
        fails = [c for c in cs if not c['passed']]
        P("  %-22s cases=%6d  d-distribution=%s  KEY TEST pass=%d fail=%d" %
          (label, len(cs), dict(sorted(dd.items())), len(cs) - len(fails), len(fails)))
    # non-reduced failures
    nrf = [c for c in cases if not c['reduced'] and not c['passed']]
    P("  non-reduced failures: %d ; of which reduced model passes: %d ; reduced model is reduced: %d"
      % (len(nrf), sum(c['red_model']['passed'] for c in nrf), sum(c['red_model']['reduced'] for c in nrf)))
    nr = [c for c in cases if not c['reduced']]
    P("  non-reduced cases: %d ; dimMT(reduced model) == dimMT check ok: %d" % (len(nr), sum(c['red_dimMT_check'] for c in nr)))
    # failing reduced cases: full data
    rf = [c for c in cases if c['reduced'] and not c['passed']]
    if rf:
        P("  !!! REDUCED FAILURES:")
        for c in rf:
            G = groups[c['gid']]
            P("   gid=%d |G|=%d orbits=%s TI=%s Phi=%s dimMT=%d d=%d LU=%s rankW=%d LW=%s index=%s adm=%s" %
              (c['gid'], G['order'], [len(o) for o in G['orbits']], G['orbit_TI'], c['Phi'], c['dimMT'], c['d'],
               c['LU'], c['rankW'], c['LW'], c['index'], [(a['K'], a['J'], a['B'], a['size']) for a in c['adm']]))
            P("     gens=%s" % G['gens'])
    else:
        P("  REDUCED FAILURES: none")
    # reduced d>=1 list
    P("  Reduced (G,Phi) with d>=1:")
    P("   %-5s %-6s %-14s %-24s %-28s %-3s %-5s %-6s %s" % ("gid", "|G|", "orbit sizes", "orbit TI", "Phi", "d", "pass", "Ncls", "admissible (K,J,B,|B|)"))
    for c in cases:
        if c['reduced'] and c['d'] >= 1:
            G = groups[c['gid']]
            P("   %-5d %-6d %-14s %-24s %-28s %-3d %-5s %-6s %s" %
              (c['gid'], G['order'], [len(o) for o in G['orbits']], ",".join(G['orbit_TI']), c['Phi'], c['d'], c['passed'],
               c['Nclass'], [(a['K'], a['J'], a['B'], a['size']) for a in c['adm']]))
    # summary by orbit-size partition for reduced non-simple d>=1
    tab = collections.defaultdict(lambda: collections.Counter())
    for c in cases:
        if c['reduced'] and c['d'] >= 1:
            G = groups[c['gid']]
            key = tuple(sorted((len(o) // 2 for o in G['orbits']), reverse=True))
            tab[key]['cases'] += 1
            tab[key]['d=%d' % c['d']] += 1
            tab[key]['pass' if c['passed'] else 'FAIL'] += 1
            tab[key]['Ncls'] = tab[key]['Ncls']
    keys = sorted(tab, key=lambda k: (-len(k), k))
    P("  Reduced d>=1 by factor-dimension partition (dims of E_i / 2):")
    for k in sorted(tab, key=lambda k: (len(k), k)):
        ncl = len(set((c['gid'], c['Nclass']) for c in cases if c['reduced'] and c['d'] >= 1 and
                      tuple(sorted((len(o) // 2 for o in groups[c['gid']]['orbits']), reverse=True)) == k))
        P("   %-16s %s  (N_W-classes: %d)" % (k, dict(tab[k]), ncl))
    return D

if __name__ == '__main__':
    out = open(sys.argv[2], 'w') if len(sys.argv) > 2 else sys.stdout
    for g in map(int, sys.argv[1].split(',')):
        summarize(g, out)
