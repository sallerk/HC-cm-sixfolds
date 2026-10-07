# verify_mutation_test.py -- mutation test: corrupt stored outputs of the main computation (enum_g*.json) and
# check that verify_enum.verify_case detects every corruption
import json, random, copy, os
from verify_enum import *
HERE = os.path.dirname(os.path.abspath(__file__))
random.seed(1)
for g in [4, 5, 6]:
    D = json.load(open(os.path.join(HERE, 'enum_g%d.json' % g)))
    groups = {G['gid']: G for G in D['groups']}
    gsd = {}
    det = tot = 0
    cases = [c for c in D['cases'] if c['d'] >= 1]
    for c in random.sample(cases, min(150, len(cases))):
        c2 = copy.deepcopy(c)
        kind = random.choice(['LU', 'passed', 'adm', 'prim', 'd', 'Phi'])
        if kind == 'LU':
            c2['LU'][0][random.randrange(g)] += 1
        elif kind == 'passed':
            c2['passed'] = not c2['passed']
        elif kind == 'adm':
            if c2['adm']: c2['adm'].pop()
            else: c2['adm'].append(dict(B=[0]))
        elif kind == 'prim':
            c2['primitive'][0] = not c2['primitive'][0]
        elif kind == 'd':
            c2['d'] += 1
        elif kind == 'Phi':
            c2['Phi'] = sorted((x + g) % (2 * g) if x == c2['Phi'][-1] else x for x in c2['Phi'])
        gid = c['gid']
        if gid not in gsd:
            gsd[gid] = group_struct(groups[gid], g)
        tot += 1
        e = verify_case(c2, groups[gid], gsd[gid], g, check_mult=False)
        det += bool(e)
        if not e: print("  undetected", g, kind, c["Phi"], c2["Phi"], c["primitive"], c["adm"][:1], c["LU"], c2["LU"])
    print("g=%d mutations %d detected %d" % (g, tot, det))
