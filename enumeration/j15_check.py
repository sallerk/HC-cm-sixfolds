# j15_check.py -- is the 6-dim sub-product A_15 x S_5 of J_15 = Jac(y^2=x^15+1) one of the failing (4,2) cases?
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from math import gcd
from sage.all import libgap
from cmenum import *
from locate import Stored
g = 6
U15 = [a for a in range(1, 15) if gcd(a, 15) == 1]          # Hom(Q(zeta15),C)
U5 = [1, 2, 3, 4]                                            # Hom(Q(zeta5),C)
pts = [('15', a) for a in U15] + [('5', b) for b in U5]
neg = lambda p: (p[0], (-p[1]) % int(p[0]))
reps = [('15', 1), ('15', 2), ('15', 4), ('15', 7), ('5', 1), ('5', 2)]   # = CM type of J_15 on these factors
lab = {}
for i, p in enumerate(reps):
    lab[p] = i; lab[neg(p)] = i + 6
assert len(lab) == 12
def mult(u):
    return libgap.PermList([lab[(p[0], (u * p[1]) % int(p[0]))] + 1 for p in sorted(lab, key=lambda q: lab[q])])
G = libgap.Group([mult(u) for u in [2, 7, 14]])
Phi = list(range(6))
gd = group_data(G, g, None, False)
for ot in gd['torbs']:
    if tuple(Phi) in ot:
        ta = type_analysis(gd, ot)
print("J15 sub-product A15 x S5: |G|=%d orbits=%s TI=%s d=%d LU=%s strict=%s primitive=%s reduced=%s" %
      (gd['order'], [len(o) for o in gd['orbits']], [o['TI'] for o in gd['oinfo']], ta['d'], ta['LU'], ta['passed'], ta['primitive'], ta['reduced']))
st = Stored(6)
gid, case = st.locate(G, Phi)
print("located gid=%s Phi=%s Nclass=%s passed=%s" % (gid, case['Phi'], case['Nclass'], case['passed']))
# also full J15 minus the elliptic part? (4,2,1) is g=7 - skip
