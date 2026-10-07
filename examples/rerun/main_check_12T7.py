"""Independent spot check of the 12T7 example (q = p = 73) with sympy and mpmath only; shares no code with
produce.py / produce_lib.gp or verify.py.
Checks: q-symmetry; irreducible over Q (sympy); |roots|^2 = q (numerically, 400 digits); ordinary;
geometric simplicity, numerically: (pi_i/pi_j)^N != 1 for i != j, N = lcm of n with phi(n) | 24 or phi(n) | 48;
multiplicative relations among the 6 arguments via mpmath.pslq (expect exactly one relation => angle rank 5),
then a search for a second relation with coefficients up to 1e12 (none expected)."""
import sympy as sp, mpmath as mp
from math import gcd
from functools import reduce
T = sp.symbols('T')
q = 73
c = [1, 36, 626, 5268, 79, -579096, -7290148, -42274008, 420991, 2049341556, 17777298866, 74630577348, 151334226289]
P = sp.Poly(c, T)
sym = all(c[12 - i] == q**(6 - i) * c[i] for i in range(7))
print('q-symmetric:', sym, ' constant = q^6:', c[12] == q**6)
print('irreducible over Q:', P.is_irreducible)
print('ordinary (q does not divide a6):', c[6] % q != 0)
mp.mp.dps = 400
roots = mp.polyroots(c, maxsteps=500, extraprec=2000)
print('max | |root|^2 - q |:', mp.nstr(max(abs(abs(r)**2 - q) for r in roots), 5))
# geometric simplicity: Q(pi^N) = Q(pi) for all N  <=>  pi_i^N distinct for all i != j (N up to roots of unity orders)
def phi(n): return sum(1 for k in range(1, n+1) if gcd(k, n) == 1)
Ns = [n for n in range(1, 2000) if 24 % phi(n) == 0 or 48 % phi(n) == 0]
N = reduce(lambda a, b: a*b//gcd(a, b), Ns)
logs = [mp.log(r) for r in roots]
mind = min(abs(mp.exp(N*(logs[i]-logs[j])) - 1) for i in range(12) for j in range(12) if i != j)
print('N =', N, ' min |(pi_i/pi_j)^N - 1| =', mp.nstr(mind, 5), '(nonzero => no ratio is a root of unity of order | N)')
# angle rank via PSLQ on arguments (one per conjugate pair) plus 2*pi/12
up = sorted([r for r in roots if mp.im(r) > 0], key=lambda z: mp.arg(z))
th = [mp.arg(z) for z in up]
rel = mp.pslq(th + [2*mp.pi/12], maxcoeff=10**6, maxsteps=10**6)
print('PSLQ relation (theta_1..6, 2pi/12):', rel)
if rel:
    print('  residual:', mp.nstr(sum(a*b for a, b in zip(rel, th + [2*mp.pi/12])), 5))
    # look for a second independent relation: drop a coordinate that the first relation uses
    for drop in range(6):
        if rel[drop] != 0:
            sub = [t for i, t in enumerate(th) if i != drop] + [2*mp.pi/12]
            rel2 = mp.pslq(sub, maxcoeff=10**12, maxsteps=10**6)
            print(f'  second relation search (drop theta_{drop+1}), coeff bound 1e12:', rel2)
            break
