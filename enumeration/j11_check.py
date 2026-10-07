# j11_check.py -- concrete (5,1) instance: J_11 = Jac(y^2 = x^11 + 1) (CM by Q(zeta11), type {1..5}) x E (CM by Q(sqrt-11))
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sage.all import libgap
from cmenum import *
from locate import Stored
g = 6
U11 = list(range(1, 11))
sq = {1, 3, 4, 5, 9}                     # squares mod 11; Q(sqrt-11) <-> kernel of Legendre symbol
pts = [('11', a) for a in U11] + [('K', 1), ('K', -1)]
reps = [('11', a) for a in [1, 2, 3, 4, 5]] + [('K', 1)]
neg = lambda p: (p[0], (-p[1]) % 11) if p[0] == '11' else ('K', -p[1])
lab = {}
for i, p in enumerate(reps):
    lab[p] = i; lab[neg(p)] = i + 6
def act(u):
    leg = 1 if u % 11 in sq else -1
    def f(p):
        return ('11', (u * p[1]) % 11) if p[0] == '11' else ('K', p[1] * leg)
    return libgap.PermList([lab[f(p)] + 1 for p in sorted(lab, key=lambda q: lab[q])])
G = libgap.Group([act(2)])
Phi = [lab[('11', a)] for a in [1, 2, 3, 4, 5]] + [lab[('K', 1)]]
gd = group_data(G, g, None, False)
for ot in gd['torbs']:
    if tuple(sorted(Phi)) in ot:
        ta = type_analysis(gd, ot)
mw = mult_weil(gd, ta, 12)
print("J11 x E_{sqrt-11} (E with CM type sigma0 = the embedding with sigma0|K = embedding fixed by squares):")
print("  |G|=%d orbits=%s TI=%s Phi=%s d=%d LU=%s strict=%s primitive=%s reduced=%s IQmultDmin=%s ex=%s" %
      (gd['order'], [len(o) for o in gd['orbits']], [o['TI'] for o in gd['oinfo']], ta['Phi'], ta['d'], ta['LU'], ta['passed'],
       ta['primitive'], ta['reduced'], mw['Dmin'], mw['examples'][:1]))
for K in gd['Ks']:
    print("  K blocks:", {o: B for o, B in K['blocks'].items()}, "t:", {o: 2*len(set(Phi) & set(B)) - len(B) for o, B in K['blocks'].items()})
st = Stored(6)
gid, case = st.locate(G, sorted(Phi))
print("  located gid=%s Phi=%s Nclass=%s passed=%s" % (gid, case['Phi'], case['Nclass'], case['passed']))
# other choice of CM type on E (the conjugate one)
Phi2 = [lab[('11', a)] for a in [1, 2, 3, 4, 5]] + [lab[('K', -1)]]
gid2, case2 = st.locate(G, sorted(Phi2))
print("  with the other CM type on E: gid=%s Phi=%s d=%s passed=%s" % (gid2, case2['Phi'], case2['d'], case2['passed']))
