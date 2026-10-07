# compare.py -- case-by-case comparison of the main computation (path1) and the independent check (path2).
# Reads path1_out.json and path2_out.json, writes compare_out.txt (all next to this script).
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
p1 = {(r['disc'], r['ti'], tuple(r['q'])): r for r in json.load(open(os.path.join(HERE, 'path1_out.json')))}
p2 = {(r['disc'], r['ti'], tuple(r['q'])): r for r in json.load(open(os.path.join(HERE, 'path2_out.json')))}
common = set(p1) & set(p2)
mism_split = [k for k in common if p1[k]['split_hilbert'] != p2[k]['split_hyperbolic']]
mism_delta = [k for k in common if p1[k]['Delta'] != p2[k]['Delta']]
mism_b = [k for k in common if p1[k]['b'] != p2[k]['b']]
pred_bad2 = [k for k in common if p2[k]['pred_formula'] != p2[k]['split_hyperbolic']]
lines = [f'cases: path1 {len(p1)}, path2 {len(p2)}, common {len(common)}',
         f'split decision mismatches (path1 Hilbert-on-det vs path2 trace-form hyperbolicity): {len(mism_split)}',
         f'Delta mismatches (Q(sqrt a) vs Q(sqrt c) embedding): {len(mism_delta)}',
         f'b mismatches: {len(mism_b)}',
         f'criterion (-b,Delta)=D vs path2 split: mismatches {len(pred_bad2)}',
         f'split cases: {sum(p2[k]["split_hyperbolic"] for k in common)}']
print('\n'.join(lines)); open(os.path.join(HERE, 'compare_out.txt'), 'w').write('\n'.join(lines) + '\n')
