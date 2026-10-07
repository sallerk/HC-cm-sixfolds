# verify_controls.py -- pure-Python re-check of the APFV rows: from APFV's own generators compute |G|, transitivity,
# primitivity, d, Lambda_U, balanced K-blocks, and compare with the stored case located by controls.py (controls_apfv.json).
import json, re, os
from verify_enum import closure_order, orbits_uf, set_orbit, rank, int_kernel, hnf, group_struct
HERE = os.path.dirname(os.path.abspath(__file__))
def parse(fn):
    t = open(fn).read()
    for m in re.finditer(r'<"([^"]+)",\s*\[(.*?)\]\s*,\s*(\d+),\s*(true|false)>', t, re.S):
        yield m.group(1), [list(map(int, x.split(','))) for x in re.findall(r'\[([^\[\]]+)\]', m.group(2))], int(m.group(3))
loc = json.load(open(os.path.join(HERE, 'controls_apfv.json')))
for g, fn in [(4, os.path.join(HERE, 'apfv_data', '4-0.000_0.000_0.000_0.000_1.00_1.00_1.00_1.00.m')),
              (6, os.path.join(HERE, 'apfv_data', '6-0.000_0.000_0.000_0.000_0.000_0.000_1.00_1.00_1.00_1.00_1.00_1.00.m'))]:
    D = json.load(open(os.path.join(HERE, 'enum_g%d.json' % g))); n = 2 * g
    groups = {G['gid']: G for G in D['groups']}
    cases = {(c['gid'], tuple(c['Phi'])): c for c in D['cases']}
    ok = 0; rows = list(parse(fn))
    for (lab, L, ar), lo in zip(rows, loc[str(g)]):
        idx = lambda s: s - 1 if s > 0 else g - s - 1
        gens = [[idx(l[i]) if i < g else idx(-l[i - g]) for i in range(n)] for l in L]
        gens = [[None] * n for _ in L]
        for k, l in enumerate(L):
            for i in range(1, g + 1):
                gens[k][idx(i)] = idx(l[i - 1]); gens[k][idx(-i)] = idx(l[g + i - 1])
        order = closure_order(gens, n)
        orbs = orbits_uf(gens, n)
        Phi = list(range(g))
        otypes = [sorted(T) for T in set_orbit(gens, Phi)]
        mu = [[1 if i in T else -1 for i in range(g)] for T in otypes]
        d = g - rank(mu)
        gs = group_struct(dict(gens=gens), g)
        prim = all(not all(set(b) <= set(Phi) or not (set(b) & set(Phi)) for b in part) for part in gs['sys'][0])
        bal = [B for K in gs['Ks'] for B in K.values() if 2 * len(set(Phi) & set(B)) == len(B)]
        st = cases[(lo['gid'], tuple(lo['Phi']))]
        G = groups[lo['gid']]
        good = (order == G['order'] and len(orbs) == 1 and G['transitive'] and d == st['d'] == 1 and prim and
                st['reduced'] and len(bal) == 1 and len(st['adm']) == 1 and st['passed'] and d == g - ar)
        ok += good
        print(g, lab, 'order', order, 'd', d, 'angle rank+1 = dimMT?', g + 1 - d == ar + 1, 'primitive', prim, '#balanced K', len(bal),
              '| stored gid', lo['gid'], 'order', G['order'], 'd', st['d'], 'passed', st['passed'], '->', 'OK' if good else 'MISMATCH')
    print("g=%d: %d/%d APFV rows consistent" % (g, ok, len(rows)))
