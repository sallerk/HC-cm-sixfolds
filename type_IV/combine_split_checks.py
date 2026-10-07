# combine_split_checks.py -- writes split_checks.out: the sum of the STATS lines and FAILS counts of the five
# chunk outputs split_checks_0.out .. split_checks_4.out (from `python split_checks.py 5 i 5`, i = 0..4).
# The combined file was first written by an inline command; this file contains the same code.
# Note: "seconds" in the combined line is the sum of the five chunk run times.
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
tot = {}
for i in range(5):
    for l in open(os.path.join(HERE, 'split_checks_%d.out' % i)):
        if l.startswith('STATS'):
            st = json.loads(l[6:])
            for k, v in st.items(): tot[k] = round(tot.get(k, 0) + v, 1)
        if l.startswith('FAILS'):
            tot['FAILS'] = tot.get('FAILS', 0) + int(l.split()[1])
with open(os.path.join(HERE, 'split_checks.out'), 'w') as f:
    f.write('COMBINED over chunks 0..4 (python split_checks.py 5 i 5; chunk i: pairs[i::5], seed 51005+i)\n')
    f.write(json.dumps(tot) + '\n')
print(tot)
