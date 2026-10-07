# analyse_sage.py -- main computation (Sage + GAP, inside the sagemath docker image).
# Reads fields_out.json (explicit Galois actions from fields_pari.py) and runs the enumeration library of the parent
# folder (cmenum.type_analysis / mult_weil, general_weil2.analyse_general2, locate.Stored with ../enum_g6.json)
# on every CM type of every config.  Run from the parent folder:  sage -python open_families/analyse_sage.py
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from sage.all import libgap
from cmenum import group_data, type_analysis, mult_weil, W_group
from general_weil2 import analyse_general2
from locate import Stored

D = json.load(open(os.path.join(HERE, 'fields_out.json')))
st = Stored(6)
W, rho = W_group(6)
res = []
for cfg in D['configs']:
    t0 = time.time()
    G = libgap.Group([libgap.PermList([x + 1 for x in p]) for p in cfg['gens']])
    gd = group_data(G, 6, None, False)
    rec = dict(name=cfg['name'], desc=cfg['desc'], order=int(gd['order']),
               orbits=[len(o) for o in gd['orbits']], TI=[o['TI'] for o in gd['oinfo']], types=[])
    for ot in gd['torbs']:
        ta = type_analysis(gd, ot)
        if not ta['reduced'] or ta['d'] == 0:
            rec['types'].append(dict(Phi=ta['Phi'], d=ta['d'], reduced=ta['reduced'], primitive=ta['primitive']))
            continue
        mw = mult_weil(gd, ta, 16)
        r = dict(Phi=ta['Phi'], d=ta['d'], LU=ta['LU'], reduced=True, strict=ta['passed'], strict_rankW=ta['rankW'],
                 strict_index=ta['index'], IQmultDmin16=mw['Dmin'], IQqrank16=mw['qrank'])
        if not ta['passed']:
            gw = analyse_general2(cfg['gens'], 6, ta['Phi'], ta['d'], ta['LU'], Dmax=12)
            r.update(genDmin12=gw['Dmin'], genQrank12=gw['qrank'], gen_fields_used=gw['fields_used'])
        gid, case = st.locate(G, ta['Phi'])
        r.update(gid=gid, stored_Phi=case['Phi'] if case else None, stored_passed=case['passed'] if case else None,
                 stored_d=case['d'] if case else None, stored_multDmin=case['mult']['Dmin'] if case and 'mult' in case else None)
        rec['types'].append(r)
    rec['secs'] = round(time.time() - t0, 1)
    res.append(rec)
    fails = [t for t in rec['types'] if t.get('reduced') and t['d'] > 0 and not t['strict']]
    print('%-22s |G|=%-4d orbits=%s TI=%s  reduced d>0: %d  strict-fail: %s' % (
        rec['name'], rec['order'], rec['orbits'], rec['TI'],
        sum(1 for t in rec['types'] if t.get('reduced') and t['d'] > 0),
        [(t['Phi'], t['d'], t['IQmultDmin16'], t.get('genDmin12'), t['gid']) for t in fails]), flush=True)
json.dump(res, open(os.path.join(HERE, 'analyse_out.json'), 'w'), indent=1, default=lambda o: int(o) if hasattr(o, '__int__') else str(o))
