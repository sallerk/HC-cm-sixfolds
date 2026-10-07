# all_delta.py -- check the expectation "every Delta < 0 occurs for every definite D" on a range:
# for each test algebra D and each squarefree m <= 60, find a diagonal T = <l1,l2,l3> (pure l_r, small
# coordinates) with Nrd(T) = m * square, i.e. Delta = -m mod squares. Also tabulate which classes are
# covered (Q(sqrt -m) embeds in D) vs not.
# Reads data.json, writes all_delta_out.txt (both next to this script).
import itertools, json, math, os
from fractions import Fraction as F
from exists import is_sq_Qp, disc_primes

HERE = os.path.dirname(os.path.abspath(__file__))

def sqfree(n):
    n = abs(n); r = 1; p = 2
    while p*p <= n:
        while n % (p*p) == 0: n //= p*p
        if n % p == 0: r *= p; n //= p
        p += 1
    return r*n

data = json.load(open(os.path.join(HERE, 'data.json')))
res = {}
lines = []
for A in data['algebras']:
    a, c, disc = A['a'], A['c'], A['disc']
    dp = disc_primes(a, c)
    norms = set()
    R = 6
    for v in itertools.product(range(-R, R+1), repeat=3):
        if any(v):
            norms.add(sqfree(-a*v[0]**2 - c*v[1]**2 + a*c*v[2]**2))
    norms = sorted(norms)
    reach = set()
    for x, y in itertools.combinations_with_replacement(norms[:60], 2):
        for z in norms[:60]:
            reach.add(sqfree(x*y*z))
    targets = [m for m in range(1, 61) if sqfree(m) == m]
    missing = [m for m in targets if m not in reach]
    covered = [m for m in targets if all(not is_sq_Qp(F(-m), p) for p in dp)]
    unc = [m for m in targets if m not in covered]
    line = (f'disc {disc}: Delta=-m realized for all squarefree m<=60 except {missing}; '
            f'uncovered (Q(sqrt-m) not in D) m: {unc}')
    print(line); lines.append(line)
open(os.path.join(HERE, 'all_delta_out.txt'), 'w').write('\n'.join(lines) + '\n')
