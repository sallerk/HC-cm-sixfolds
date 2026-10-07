# run_enum.py -- main computation: enumerate W(B_g)-classes of subgroups G containing rho, all CM types up to G,
# KEY TEST.   usage (Sage):  sage -python run_enum.py g      -> enum_g{g}.json in this folder
import sys, os, json, time, collections
from sage.all import libgap
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cmenum import *

g = int(sys.argv[1])
t0 = time.time()
n = 2 * g
W, rho = W_group(g)
cc = libgap.ConjugacyClassesSubgroups(W)
reps = [libgap.Representative(c) for c in cc]
reps = [H for H in reps if bool(libgap.IsSubset(H, [rho]))]
print("g=%d |W|=%d classes containing rho: %d  (%.1fs)" % (g, int(libgap.Size(W)), len(reps), time.time() - t0), flush=True)

groups_out, cases_out = [], []
for gid, G in enumerate(reps):
    gd, res = analyse_group(G, g, W, True)
    groups_out.append(dict(gid=gid, gens=gd['gens'], order=gd['order'], orbits=gd['orbits'],
                           orbit_TI=[o['TI'] for o in gd['oinfo']], orbit_img_order=[o['img_order'] for o in gd['oinfo']],
                           transitive=(len(gd['orbits']) == 1),
                           Ks=[dict(signs=K['signs'], blocks={str(k): v for k, v in K['blocks'].items()}) for K in gd['Ks']],
                           ntypeorbits=len(gd['torbs'])))
    for ot, ta in zip(gd['torbs'], res):
        c = dict(gid=gid, Phi=ta['Phi'], norb=ta['norb'], dimMT=ta['dimMT'], d=ta['d'], LU=ta['LU'],
                 adm=[dict(K=a['K'], J=a['J'], B=a['B'], size=a['size']) for a in ta['adm']],
                 rankW=ta['rankW'], LW=ta['LW'], index=ta['index'], satidx=ta['satidx'], passed=ta['passed'],
                 primitive=ta['primitive'], isog=ta['isog'], reduced=ta['reduced'], Nclass=ta.get('Nclass'))
        if ta['reduced'] and ta['d'] >= 1:
            c['mult'] = mult_weil(gd, ta, 12)
        if not ta['reduced']:
            r2, gp = reduced_dimMT(gd, ta, ot)
            c['red_g'] = gp
            c['red_dimMT_check'] = (r2 == ta['dimMT'])
            if not ta['passed']:
                gd2, ta2, Phi2 = analyse_reduced_model(gd, ta)
                c['red_model'] = dict(g=gp, gens=gd2['gens'], Phi=ta2['Phi'], dimMT=ta2['dimMT'], d=ta2['d'],
                                      passed=ta2['passed'], reduced=ta2['reduced'], index=ta2['index'],
                                      rankW=ta2['rankW'], orbit_sizes=[len(o) for o in gd2['orbits']],
                                      orbit_TI=[o['TI'] for o in gd2['oinfo']])
        cases_out.append(c)
    if gid % 200 == 0:
        print("  group %d/%d  cases so far %d  (%.1fs)" % (gid, len(reps), len(cases_out), time.time() - t0), flush=True)

json.dump(dict(g=g, groups=groups_out, cases=cases_out), open(os.path.join(HERE, 'enum_g%d.json' % g), 'w'), default=int)
print("done g=%d: groups %d, cases %d, time %.1fs" % (g, len(groups_out), len(cases_out), time.time() - t0), flush=True)
