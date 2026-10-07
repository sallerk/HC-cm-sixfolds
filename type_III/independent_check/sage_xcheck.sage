# Sage cross-check of indep_check.py: (1) Hilbert symbols, (2) norm test -detH in Nm(Q(sqrt -b)) via Sage/PARI is_norm,
# (3) criterion of Proposition 6.2 of the note ("Thm3") via Sage QuaternionAlgebra ramification.
# Inputs (next to this script): sage_hilb_input.txt (rows a b p h, h = indep_check.hilbert(a,b,p), p = 0 the real place),
# sage_norm_input.txt (rows b detH split Delta al be thm3, taken from run_*.json).
import os
HERE = os.path.dirname(os.path.abspath(__file__))
bad=0; n=0
for line in open(os.path.join(HERE, 'sage_hilb_input.txt')):
    a,b,p,h = [ZZ(t) for t in line.split()]
    hs = hilbert_symbol(a,b,-1 if p==0 else p)
    n+=1
    if hs != h: bad+=1; print("HILB MISMATCH", a,b,p,h,hs)
print("hilbert symbols checked", n, "mismatches", bad)
bad=0; n=0; bad3=0
for line in open(os.path.join(HERE, 'sage_norm_input.txt')):
    b,d,sp,Dl,al,be,t3 = line.split()
    b=ZZ(b); d=QQ(d); sp=bool(int(sp)); Dl=QQ(Dl); al=ZZ(al); be=ZZ(be); t3=bool(int(t3))
    K = QuadraticField(-b, 'w')
    isn = (-d).is_norm(K)
    n+=1
    if isn != sp: bad+=1; print("NORM MISMATCH", b, d, sp, isn)
    # Proposition 6.2 via Sage quaternion algebras: (-b, Delta) ~ (al, be) ?
    Q1 = QuaternionAlgebra(QQ, -b, Dl.numerator()*Dl.denominator()); Q2 = QuaternionAlgebra(QQ, al, be)
    same = (Q1.ramified_primes() == Q2.ramified_primes())   # both definite here
    if same != sp: bad3+=1; print("THM3 MISMATCH", b, d, Dl, al, be, sp, same)
print("norm rows checked", n, "split-decision mismatches", bad, "Thm3(Sage quaternion) mismatches", bad3)
