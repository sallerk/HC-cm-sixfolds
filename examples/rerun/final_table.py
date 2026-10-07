"""Assemble final table: main computation (produce_runs/produce_*.json, matched via examples_in.json) + independent
verifier (verify_A/B.json) + GAP (verify_gap.json) + exact K-certificate (verify_K.json).
Output final_table.txt / final_table.json (all paths relative to this script's directory)."""
import json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
EX = {e['name']: e for e in json.load(open(os.path.join(HERE, 'examples_in.json')))}
V = {r['name']: r for f in ('verify_A.json', 'verify_B.json') for r in json.load(open(os.path.join(HERE, f)))}
GP = {r['name']: r for r in json.load(open(os.path.join(HERE, 'verify_gap.json')))}
KC = {r['name']: r for r in json.load(open(os.path.join(HERE, 'verify_K.json')))}
PR = {}
import glob
for fn in glob.glob(os.path.join(HERE, 'produce_runs', 'produce_6T*.json')):
    for r in json.load(open(fn)):
        PR[(r['polE'], r['p'], tuple(r['S']))] = r
out, lines = [], []
for name, e in EX.items():
    v = V[name]; g = GP[name]; k = KC[name]
    pr = PR.get((e.get('polE'), e['p'], tuple(e.get('S', ())))) if 'polE' in e else None
    rec = dict(name=name, row=e['row'], p=e['p'], q=e['q'], q_exp=round(math.log(e['q'], e['p'])),
               polF=e.get('polF'), polK=e.get('polK'), polE=e.get('polE'), S=e.get('S'), weil_poly=e['weil_poly'],
               main_computation=None if pr is None else dict(angle_rank=pr['angle_rank'], rank_info=pr['rank_info'],
                    geom=pr['geom_simple'], ordinary=pr['ordinary'], irreducible=pr['irreducible'],
                    norm_relation=pr['norm_relation'], abs_ok=pr['abs_ok'], sym=pr['q_symmetric']),
               verifier=dict(ALL_OK=v['ALL_OK'], angle_rank_numerical=v['angle_rank_numerical'], relations=v['relations'],
                    no_further_relation_below=v['bound_no_further_relation'],
                    relation_certified=[c['certified'] for c in v['relation_certificates']]),
               gap_h=g['GaloisType_h'], gap_P=g['GaloisType_P'], K_d=k['d'], label_12T=k['label_12T'], label_ok=k['label_matches'])
    pa = rec['main_computation']['angle_rank'] if rec['main_computation'] else None
    rec['VERIFIED'] = bool(v['ALL_OK'] and k['label_matches'] and (pa is None or pa == v['angle_rank_numerical'])
                          and all(rec['verifier']['relation_certified']))
    out.append(rec)
    lines.append(f"{name:30s} {e['row']:7s} p={e['p']:<6d} q=p^{rec['q_exp']} K=Q(sqrt-{k['d']}) h:{g['GaloisType_h']:5s} "
                 f"main_ar={pa} verifier_ar={v['angle_rank_numerical']} cert={rec['verifier']['relation_certified']} "
                 f"geom={v['geom_simple']} ord={v['ordinary']} irr={v['irreducible']} weil={v['weil_sturm']} VERIFIED={rec['VERIFIED']}")
    lines.append(f"      F={e.get('polF')}  S={e.get('S')}")
    lines.append(f"      P(T) coeffs T^12..T^0 = {e['weil_poly']}")
json.dump(out, open(os.path.join(HERE, 'final_table.json'), 'w'), indent=1)
open(os.path.join(HERE, 'final_table.txt'), 'w').write('\n'.join(lines) + '\n')
print('\n'.join(l for l in lines if not l.startswith('      ')))
