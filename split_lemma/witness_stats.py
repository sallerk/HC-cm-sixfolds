# witness_stats.py -- summary statistics of split_witnesses.json (output of verify_split.py) that are not
# printed by verify_split.py: composition of the test algebras, witnesses per search phase, witnesses that
# needed auxiliary primes, forced primes, and the split fraction among correctly-signed random elements a.
# Pure Python; reads the JSON only. The last line repeats the inline command that first produced it.
import json, os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'split_witnesses.json')))
pats = [p for c in D for p in c['patterns']]

print('cases', len(D))
print('cases by factor degrees', sorted(Counter(tuple(c['degs']) for c in D).items()))
fields = [c for c in D if len(c['degs']) == 1]
cyc = [c for c in fields if not c['label'].startswith('random')]
rnd = [c for c in fields if c['label'].startswith('random')]
print('sextic-field cases', len(fields), '= Q(zeta_m)^+ cases', len(cyc), sorted(set(c['label'] for c in cyc)),
      '+ random sextic cases', len(rnd), '(distinct fields: %d)' % len(set(tuple(c['polys']) for c in rnd)))
print('product cases', len(D) - len(fields))
print('d values', sorted(set(c['d'] for c in D)), ' max d_F', max(int(c['dF']) for c in D))
print('witnesses', sum(p['witness'] is not None for p in pats), 'of', len(pats), 'sign patterns')
print('witnesses by search phase', sorted(Counter(p['phase'] for p in pats).items()))
aux = [len(p['witness'].get('aux_primes', [])) for p in pats if p['witness'] is not None]
print('witnesses needing auxiliary primes', sum(1 for a in aux if a > 0),
      ' (number of auxiliary primes, count):', sorted(Counter(a for a in aux if a > 0).items()))
F = [f for c in D for f in c['forced'] if f['forced']]
print('forced primes', len(F), ' with t_p = 1:', sum(f['t_p'] == 1 for f in F),
      ' with t_p = -1:', sum(f['t_p'] == -1 for f in F))
su = [p['n_signed_tried'] - 1 for p in pats if p['phase'] == 'sunit']
print('S-unit-phase patterns', len(su), '; correctly-signed random a tried before (none split): min', min(su),
      'median', sorted(su)[len(su) // 2], 'max', max(su))
fr = []
for c in D:
    s = sum(p['n_signed_tried'] for p in c['patterns'] if p['phase'] == 'random')
    k = sum(p['n_split'] for p in c['patterns'] if p['phase'] == 'random')
    if s >= 50: fr.append(k / s)
print(len(fr), 'cases; split fraction among correctly-signed random a (random-phase patterns): min %.3f median %.3f max %.3f'
      % (min(fr), sorted(fr)[len(fr) // 2], max(fr)))
