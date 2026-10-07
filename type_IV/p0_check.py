# p0_check.py -- closed-formula proof of the split lemma (note, Sec. 7: the element a0):
#   a0 = -sqrt(m) * Delta,  Delta = det_L(H) in F,  sqrt(m) taken positive at the place f1 where H has
#   signature (2,1).  Then a0 >> 0 and det(h_a0) = N(a0)^3 N(Delta) (4m)^3 = -m^3 N(Delta)^4 (4m)^3,
#   which is -1 times a rational square, so h_a0 is hyperbolic by Landherr's theorem.
# Checked here on fresh random general Hermitian H (all 156 (d,m) pairs incl. d = m) with the independent
# checker (explicit Gram determinant; Hilbert symbols; Witt index via qfsolve), plus the alternative
# a0' = -sqrt(m)/Delta.
import sys, json, time, random
from fractions import Fraction as Fr
from biq import Ctx
import split_checks as S

S.rng.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 424242)
per = int(sys.argv[1]) if len(sys.argv) > 1 else 3
tot = dict(cases=0, totpos=0, det_is_minus_square=0, checker_ok=0, checker_ok_inv=0)
bad = []
t0 = time.time()
for d in [1, 2, 3, 5, 6, 7, 11, 15, 19, 23, 31, 35]:
    for m in [2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 21, 30]:
        C = Ctx(d, m)
        for _ in range(per):
            H = S.good_herm(C)
            if H is None:
                continue
            tot['cases'] += 1
            D = C.detL(H)
            a0 = C.neg(C.mul(C.Fel(0, 1), D))
            a1 = C.neg(C.div(C.Fel(0, 1), D))
            tp = C.totpos(a0) and C.totpos(a1)
            tot['totpos'] += int(tp)
            G = C.gram_K(H, a0)
            det = C.detK(G)
            q = -det
            from cypari import pari
            sq = q > 0 and bool(pari("issquare(%d)" % (q.numerator * q.denominator)))
            tot['det_is_minus_square'] += int(sq)
            c0 = S.checker(C, H, a0); c1 = S.checker(C, H, a1)
            ok0 = c0['det_ok'] and c0['hilbert_split'] and c0['witt'] == 6
            ok1 = c1['det_ok'] and c1['hilbert_split'] and c1['witt'] == 6
            tot['checker_ok'] += int(ok0); tot['checker_ok_inv'] += int(ok1)
            if not (tp and sq and ok0 and ok1):
                bad.append((d, m, tp, sq, c0, c1))
tot['seconds'] = round(time.time() - t0, 1)
print(json.dumps(tot))
print('bad', len(bad))
for b in bad[:20]:
    print(b)
