# cor4_check.py -- constructive test of the second statement of Proposition 6.2 of the note, direction (<=)
# (printed as "Cor4": if Q(sqrt Delta) embeds in D, the K = Q(x) of the proof is split), and hand examples
# E1/E2/Abdulali-type cases (every K = Q(q), q primitive pure with |coords| <= 3: split == criterion).
# Uses only indep_check.py (no PARI). Output on stdout (saved as cor4_check_out.txt).
import random, sys, math, itertools
from fractions import Fraction as Fr
from indep_check import *

def is_square_rat(x):
    x = Fr(x)
    if x <= 0: return False
    return math.isqrt(x.numerator)**2 == x.numerator and math.isqrt(x.denominator)**2 == x.denominator

def find_y(al, be, N, bound=12):
    # pure y with Nrd(y) in N * Q^x2  (then y^2 = -Nrd(y) in Delta Q^x2)
    for B in range(1, bound+1):
        for y in itertools.product(range(-B, B+1), repeat=3):
            if max(map(abs, y)) != B: continue
            yy = (0,)+y
            n = qnrd(yy, al, be)
            if n and is_square_rat(Fr(n) / N): return yy
    return None

def split_for(S, q):
    al, be = S.al, S.be
    b = qnrd(q, al, be)
    dH, _, _ = hermitian_det(S, q, b)
    assert dH[1] == 0
    return is_norm_from(-b, -dH[0]), dH[0], b

def constructive(seed, ncases):
    random.seed(seed)
    stats = dict(cases=0, embedded=0, constructed=0, split_ok=0, split_fail=0, y_not_found=0, notemb=0, notemb_anysplit=0)
    for n in range(ncases):
        al, be = random.choice(ALGS) if n % 2 else (-random.randint(1,40), -random.randint(1,40))
        M = rand_M(random.choice(['diag','sparse','dense']))
        S = Setup(al, be, M)
        if det_Q(S.G) == 0: continue
        N = S.nrdT(); Delta = -N
        emb, ram = embeds(Delta, al, be)
        stats['cases'] += 1
        if emb:
            stats['embedded'] += 1
            y = find_y(al, be, N)
            if y is None: stats['y_not_found'] += 1; continue
            x = anticommuting_pure(y, al, be)
            g = math.gcd(math.gcd(x[1], x[2]), x[3]); x = (0, x[1]//g, x[2]//g, x[3]//g)
            sp, dH, b = split_for(S, x)
            stats['constructed'] += 1
            if sp: stats['split_ok'] += 1
            else:
                stats['split_fail'] += 1; print("FAIL", al, be, M, Delta, y, x, dH)
        else:
            stats['notemb'] += 1
            anysplit = False
            for _ in range(12):
                q = rand_pure(-4, 4)
                g = math.gcd(math.gcd(q[1], q[2]), q[3]); q = (0, q[1]//g, q[2]//g, q[3]//g)
                if split_for(S, q)[0]: anysplit = True
            stats['notemb_anysplit'] += anysplit
    print("constructive Cor4 (seed %d):" % seed, stats)

def all_K(S, bound):
    out = []
    seen = set()
    for q in itertools.product(range(-bound, bound+1), repeat=3):
        if q == (0,0,0): continue
        g = math.gcd(math.gcd(q[0], q[1]), q[2])
        if g != 1: continue
        if tuple(-c for c in q) in seen: continue
        seen.add(q)
        qq = (0,)+q
        sp, dH, b = split_for(S, qq)
        out.append((qq, b, sp))
    return out

def diagM(l1, l2, l3):
    M = [[Z4]*3 for _ in range(3)]; M[0][0], M[1][1], M[2][2] = l1, l2, l3; return M

def examples():
    I, J, K_ = (0,1,0,0), (0,0,1,0), (0,0,0,1)
    cases = [
        ("E1  D=(-2,-5) disc5, T=<i,j,k>", -2, -5, diagM(I, J, K_)),
        ("E2  D=(-1,-7) disc7, T=<i,j,k>", -1, -7, diagM(I, J, K_)),
        ("D=(-1,-3) disc3, T=<i,j,k>", -1, -3, diagM(I, J, K_)),
        ("D=(-1,-1) disc2, T=<i,j,k>", -1, -1, diagM(I, J, K_)),
        ("D=(-3,-5) disc5, T=<i,j,k>", -3, -5, diagM(I, J, K_)),
        ("D=(-1,-3), T=<i,i,j> (Abdulali Ex.3.1.3, m=3)", -1, -3, diagM(I, I, J)),
    ]
    for name, al, be, M in cases:
        S = Setup(al, be, M); N = S.nrdT(); Delta = -N
        emb, ram = embeds(Delta, al, be)
        res = all_K(S, 3)
        spl = sorted(set(b for (q, b, sp) in res if sp))
        # square-free parts of -b for the split K, and which classical K (Q(i), Q(sqrt-3)) occur, split or not
        def sqf(n):
            f = factor(n); r = 1
            for p, e in f.items():
                if e % 2: r *= p
            return r
        fields_all = sorted(set(sqf(b) for (q, b, sp) in res))
        fields_split = sorted(set(sqf(b) for (q, b, sp) in res if sp))
        print(f"{name}: ramified {ram}, Nrd(T)={N}, Delta={Delta}, Q(sqrt Delta) embeds: {emb}; "
              f"#K tested {len(res)}; d with Q(sqrt-d) tested: {fields_all[:15]}...; split d: {fields_split}")
        # consistency: split iff (-b,Delta) = D for each
        for (q, b, sp) in res:
            t3 = all(hilbert(-b, Delta, v) == hilbert(al, be, v) for v in places(-b, Delta, al, be))
            assert t3 == sp, (name, q, b, sp, t3)
    print("hand examples: split == Thm3 for every K tested")

if __name__ == '__main__':
    examples()
    for s in (11, 12):
        constructive(s, 150)
