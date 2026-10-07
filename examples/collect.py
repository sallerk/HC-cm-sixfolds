"""Collect, per APFV row (6Tj x C2), the example of the main computation (produce.py) with the smallest q among those
passing every check of the main computation (balanced, geometrically simple, ordinary, irreducible, q-symmetric,
|roots| ok, norm relation, angle rank 5), plus one non-degenerate control (unbalanced, geometrically simple,
angle rank 6) from the same field.
Reads produce_runs/produce_6T*.json and jacobi21.json; writes examples_in.json (verifier input) and examples_table.txt
(all paths relative to this script's directory)."""
import json, glob, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROW = {1: '12T2', 3: '12T10', 4: '12T7', 5: '12T18', 6: '12T25', 7: '12T23', 8: '12T24', 9: '12T37', 10: '12T40',
       11: '12T48', 12: '12T75', 13: '12T77', 14: '12T123', 15: '12T180', 16: '12T219'}
GORD = {1: 12, 3: 24, 4: 24, 5: 36, 6: 48, 7: 48, 8: 48, 9: 72, 10: 72, 11: 96, 12: 120, 13: 144, 14: 240, 15: 720, 16: 1440}
allr = []
for fn in glob.glob(os.path.join(HERE, 'produce_runs', 'produce_6T*.json')):
    for r in json.load(open(fn)):
        r['file'] = os.path.basename(fn); allr.append(r)
def ok(r):
    return (r['q_symmetric'] and r['abs_ok'] and r['irreducible'] and r['ordinary'] and r['geom_simple'])
best, ctrl = {}, {}
for r in allr:
    j = r['jF']
    if ok(r) and r['balanced'] and r['norm_relation'] and r['angle_rank'] == 5:
        if j not in best or r['q'] < best[j]['q']: best[j] = r
for j, b in best.items():
    c = [r for r in allr if r['polE'] == b['polE'] and ok(r) and not r['balanced'] and r['angle_rank'] == 6 and not r['norm_relation']]
    if c: ctrl[j] = min(c, key=lambda r: (r['p'] != b['p'], r['q']))
ex = []
lines = []
for j in sorted(best):
    b = best[j]
    ex.append(dict(name=f'{ROW[j]}_deg', row=ROW[j], jF=j, weil_poly=b['weil_poly'], p=b['p'], q=b['q'], G_order=GORD[j],
                   polF=b['polF'], polK=b['polK'], polE=b['polE'], S=b['S'], file=b['file'], pi=b['pi']))
    lines.append(f"{ROW[j]:7s} 6T{j:<2d} K={b['polK']:9s} p={b['p']:<7d} q={b['p']}^{round(__import__('math').log(b['q'], b['p']))} "
                 f"S={b['S']} angle_rank={b['angle_rank']} ({b['rank_info']}) file={b['file']}\n        P={b['weil_poly']}")
    if j in ctrl:
        c = ctrl[j]
        ex.append(dict(name=f'{ROW[j]}_ctrl', row=ROW[j], jF=j, weil_poly=c['weil_poly'], p=c['p'], q=c['q'], G_order=GORD[j],
                       polF=c['polF'], polK=c['polK'], polE=c['polE'], S=c['S'], file=c['file'], pi=c['pi']))
        lines.append(f"   ctrl: p={c['p']} q={c['q']} S={c['S']} angle_rank={c['angle_rank']}")
J = json.load(open(os.path.join(HERE, 'jacobi21.json')))
for r in J:
    if r['irreducible']:
        ex.append(dict(name=f"fermat21_jacobi_p{r['p']}_{'_'.join(map(str, r['abc']))}", row='12T2', jF=1, weil_poly=r['weil_poly'],
                       p=r['p'], q=r['p'], G_order=12))
json.dump(ex, open(os.path.join(HERE, 'examples_in.json'), 'w'), indent=1)
open(os.path.join(HERE, 'examples_table.txt'), 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines)); print(len(ex), 'verifier inputs')
