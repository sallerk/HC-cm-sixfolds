# nonsplit_loci.py -- side remark (note, Sec. 7): every (K, 3, delta) Weil family (any discriminant)
# contains 4-dimensional type-IV loci whose members ALSO carry a split L-compatible polarization.
# Construction (argued by hand, checked here): delta < 0 a rational representative of the class (Markman's
# normalization: det H mod Nm K^x; split <-> -1).  Pick y >= 1 with M = y^2 + delta > 0 not a square,
# F = Q(sqrt M), Delta = M + y sqrt M (N(Delta) = M delta < 0), H = diag(1, -1, -Delta) over L = K F
# (signatures (2,1) / (1,2) after labelling the places).  Then det(h_1) = N(Delta) (4M)^3 = delta mod squares,
# so the members with polarization E (a = 1) lie in the (K,3,delta)-family, while a0 = -sqrt(M) Delta gives a
# split polarization (checked with the independent checker).
import json, os
from fractions import Fraction as Fr
from cypari import pari
from biq import Ctx
import split_checks as S

out = []
for d in [1, 2, 3, 5, 7, 11]:
    for delta in [-1, -2, -3, -5, -6, -7, -10, -11, -13, -14, -15, -17, -19, -21, -22, -23]:
        y = 1
        while True:
            M = y*y + delta
            if M > 1 and not pari('issquare(%d)' % M):
                break
            y += 1
        C = Ctx(d, M)
        Dl = C.Fel(M, y)
        Z = C.L()
        H = [[C.Fel(1), Z, Z], [Z, C.Fel(-1), Z], [Z, Z, C.neg(Dl)]]
        if C.signature(H, 1) != (2, 1):
            H = [[C.rho(H[i][j]) for j in range(3)] for i in range(3)]
            Dl = C.rho(Dl)
        assert C.signature(H, 1) == (2, 1) and C.signature(H, 2) == (1, 2)
        det1 = C.detK(C.gram_K(H, C.Fel(1)))
        same_class = C.is_norm_K(Fr(delta) * det1)          # det(h_1) == delta mod Nm K^x
        nonsplit = not C.is_norm_K(Fr(-delta))               # is the class delta non-split?
        a0 = C.neg(C.mul(C.Fel(0, 1), C.detL(H)))
        ch = S.checker(C, H, a0)
        out.append(dict(d=d, delta=delta, M=M, y=y, delta_nonsplit=nonsplit, h1_in_class=same_class,
                        a0_totpos=C.totpos(a0), a0_split=(ch['det_ok'] and ch['hilbert_split'] and ch['witt'] == 6)))
ok = all(r['h1_in_class'] and r['a0_totpos'] and r['a0_split'] for r in out)
print('cases', len(out), 'nonsplit classes among them', sum(r['delta_nonsplit'] for r in out), 'ALL OK' if ok else 'FAIL')
for r in out[:12]:
    print(r)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nonsplit_loci.json'), 'w'), indent=1)
