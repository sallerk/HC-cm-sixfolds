# tiers.py -- local python: tier classification of reduced d>=1 cases, from the enum_g*.json files of run_enum.py
#   usage:  python tiers.py 4,5,6   (a second argument also prints every non-T0 case)
#   tiers (labels of the output, not related to anything else): T0 = strict KEY TEST passes; otherwise, for the
#   imaginary-quadratic Weil characters with multiplicities (cmenum.mult_weil, Dmax = 12): T1 = they generate
#   Lambda_U at dimension D <= 6; T2 = first at some D with 6 < D <= 12; T3 = not up to D = 12.
import json, sys, os, collections
HERE = os.path.dirname(os.path.abspath(__file__))
for g in map(int, sys.argv[1].split(',')):
    D = json.load(open(os.path.join(HERE, 'enum_g%d.json' % g)))
    groups = {G['gid']: G for G in D['groups']}
    print("=" * 80); print("g =", g)
    tab = collections.Counter()
    rows = []
    for c in D['cases']:
        if not (c['reduced'] and c['d'] >= 1):
            continue
        G = groups[c['gid']]
        part = tuple(sorted((len(o) // 2 for o in G['orbits']), reverse=True))
        m = c['mult']
        tier = 'T0 strict' if c['passed'] else ('T1 mult D<=6' if m['Dmin'] is not None and m['Dmin'] <= 6 else
               ('T2 mult D=%s' % m['Dmin'] if m['Dmin'] is not None else 'T3 none<=%d (Qrank %d of d=%d)' % (m['Dmax'], m['qrank'], c['d'])))
        tab[(part, c['d'], tier)] += 1
        rows.append((part, c['gid'], G['order'], G['orbit_TI'], c['Phi'], c['d'], tier, c['Nclass'], c['LU'], m['examples'][:3]))
    for k in sorted(tab):
        print("  ", k, tab[k])
    if len(sys.argv) > 2:
        for r in rows:
            if not r[6].startswith('T0'):
                print("   ", r)
