# check_prop41.py -- pure Python 3, standard library only; shares no code with the enumeration programs.
# Checks the steps of the proof of Proposition 4.1 of the note (Lemmas 4.2-4.5) by direct computation:
#   A. every transitive pair (G, Phi) with g = 6 in ../enum_g6.json (the 142 W(B_6)-classes of transitive G):
#      Lemma 4.2 (2g <= 2^s for primitive Phi, s = dim of the span of the Galois conjugates of mu_Phi; and if two of
#      the functionals +-x_i agree on that span, Phi is imprimitive); d <= 1 for primitive Phi; Lemma 4.3 for sextic
#      CM subfields (Phi with one element in each fibre is imprimitive) and the parity statement (no Phi has
#      signature (3,3) with respect to two imaginary quadratic subfields); counts of primitive types with d = 1.
#   B. Lemma 4.4: E Galois with group S_3 x C_2 (12 points, regular action): the formula for d in terms of the
#      2x2 matrix A = sum eps_sigma std(sigma), and "rank A <= 1 => Phi imprimitive".
#   C. Lemma 4.5 and Table 1: the 16 transitive groups H_0 of degree 6, built from their definitions; block systems,
#      the criterion "P is a block of size 3 or a transversal of a system of blocks of size 2", the number n of
#      primitive 3-subsets, orbits, and the comparison with ../../examples/combi_rows.json and ../enum_g6.json.
#   usage:  python check_prop41.py > check_prop41_out.txt
import json, os, itertools
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ENUM = os.path.dirname(HERE)
REPO = os.path.dirname(ENUM)
FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print('FAIL:', msg)


# ---------------------------------------------------------------- permutation groups (perm = tuple, p[i] = image)
def mul(a, b):                      # (a*b)(i) = a(b(i))
    return tuple(a[b[i]] for i in range(len(b)))


def closure(gens, n):
    e = tuple(range(n))
    seen, frontier = {e}, [e]
    gens = [tuple(g) for g in gens]
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = mul(g, x)
                if y not in seen:
                    seen.add(y)
                    new.append(y)
        frontier = new
    return seen


def set_orbit(gens, S):
    S = frozenset(S)
    seen, frontier = {S}, [S]
    while frontier:
        new = []
        for T in frontier:
            for g in gens:
                U = frozenset(g[i] for i in T)
                if U not in seen:
                    seen.add(U)
                    new.append(U)
        frontier = new
    return seen


def is_transitive(gens, n):
    return len(set_orbit(gens, [0])) == n


def block_systems(gens, n):
    """All non-trivial block systems of a transitive group, as lists of frozensets."""
    out = []
    for k in range(2, n):
        if n % k:
            continue
        for rest in itertools.combinations(range(1, n), k - 1):
            B = frozenset((0,) + rest)
            orb, ok, frontier = {B}, True, [B]
            while frontier and ok:
                new = []
                for T in frontier:
                    for g in gens:
                        U = frozenset(g[i] for i in T)
                        if U not in orb:
                            if any(U & V for V in orb):
                                ok = False
                                break
                            orb.add(U)
                            new.append(U)
                    if not ok:
                        break
                frontier = new
            if ok:
                out.append(sorted(orb, key=min))
    return out


def rank(rows):
    M = [[Fraction(x) for x in r] for r in rows]
    r = 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r


# ---------------------------------------------------------------- CM types on 2g points, rho(i) = i + g mod 2g
def cm_types(g):
    for bits in itertools.product((0, 1), repeat=g):
        yield frozenset(i + g * b for i, b in enumerate(bits))


def mu(Phi, g):
    return tuple(1 if i in Phi else -1 for i in range(g))


def analyse(gens, g, Phi, systems):
    """s = dim span of the conjugates of mu_Phi; matrix rows; primitivity (via block systems with rho free)."""
    orbit = set_orbit(gens, Phi)
    rows = [mu(P, g) for P in orbit]
    s = rank(rows)
    prim = True
    for sysm in systems:
        B0 = sysm[0]
        rhoB0 = frozenset((i + g) % (2 * g) for i in B0)
        if rhoB0 == B0:
            continue                               # totally real subfield
        if all(B <= Phi or not (B & Phi) for B in sysm):
            prim = False
            break
    return s, rows, prim


def iq_blocks(systems, g):
    """Imaginary quadratic subfields K: systems of two blocks of size g interchanged by rho (block of index-0 side)."""
    out = []
    for sysm in systems:
        if len(sysm) == 2:
            B = sysm[0]
            if frozenset((i + g) % (2 * g) for i in B) == sysm[1]:
                out.append(B)
    return out


# ================================================================ A. the 142 transitive classes with g = 6
g = 6
D = json.load(open(os.path.join(ENUM, 'enum_g6.json')))
trans = [G for G in D['groups'] if G['transitive']]
print('A. transitive W(B_6)-classes of subgroups containing rho:', len(trans))
check(len(trans) == 142, 'expected 142 transitive classes')
stats = dict(pairs=0, prim=0, prim_d1=0, prim_d1_orbits=0, max_d_prim=0, min_s_prim=99, func_coincide_prim=0,
             func_coincide_imprim_ok=0, sextic_systems=0, sextic_balanced_types=0, sextic_balanced_prim=0,
             two_balanced_K=0)
stored = {}
for c in D['cases']:
    stored.setdefault(c['gid'], []).append(c)
classes_d1 = 0
for G in trans:
    gens = [tuple(x) for x in G['gens']]
    rho = tuple((i + g) % (2 * g) for i in range(2 * g))
    for x in gens:
        check(mul(x, rho) == mul(rho, x), 'gid %d: generator does not commute with rho' % G['gid'])
    check(is_transitive(gens, 2 * g), 'gid %d: stored as transitive, but G is not transitive' % G['gid'])
    check(rho in closure(gens, 2 * g), 'gid %d: rho is not in G' % G['gid'])
    systems = block_systems(gens, 2 * g)
    Ks = iq_blocks(systems, g)
    sextic = [s_ for s_ in systems if len(s_) == 6 and frozenset((i + g) % 12 for i in s_[0]) != s_[0]]
    stats['sextic_systems'] += len(sextic)
    seen_orbits = set()
    for Phi in cm_types(g):
        stats['pairs'] += 1
        s, rows, prim = analyse(gens, g, Phi, systems)
        d = g - s
        cols = list(zip(*rows))
        coincide = any(cols[i] == cols[j] or cols[i] == tuple(-x for x in cols[j])
                       for i in range(g) for j in range(i + 1, g))
        if prim:
            stats['prim'] += 1
            stats['max_d_prim'] = max(stats['max_d_prim'], d)
            stats['min_s_prim'] = min(stats['min_s_prim'], s)
            check(2 * g <= 2 ** s, 'Lemma 4.2 bound fails: gid %d Phi %s' % (G['gid'], sorted(Phi)))
            check(not coincide, 'Lemma 4.2: functionals coincide for primitive Phi, gid %d' % G['gid'])
            check(d <= 1, 'd >= 2 for a primitive type: gid %d Phi %s' % (G['gid'], sorted(Phi)))
            if coincide:
                stats['func_coincide_prim'] += 1
            if d == 1:
                stats['prim_d1'] += 1
                key = frozenset(set_orbit(gens, Phi))
                if key not in seen_orbits:
                    seen_orbits.add(key)
                    stats['prim_d1_orbits'] += 1
        elif coincide:
            stats['func_coincide_imprim_ok'] += 1
        bal = [K for K in Ks if len(K & Phi) == len(K) // 2]
        check(len(bal) <= 1, 'parity: two balanced K, gid %d' % G['gid'])
        stats['two_balanced_K'] += (len(bal) >= 2)
        for sysm in sextic:
            if all(len(B & Phi) == 1 for B in sysm):
                stats['sextic_balanced_types'] += 1
                if prim:
                    stats['sextic_balanced_prim'] += 1
                check(not prim, 'Lemma 4.3 (sextic): primitive Phi with one element per fibre, gid %d' % G['gid'])
    # comparison with the stored records (G-orbits of CM types)
    st = [c for c in stored.get(G['gid'], []) if all(c['primitive']) and c['d'] >= 1]
    check(len(st) == len(seen_orbits), 'gid %d: stored primitive d>=1 orbits %d vs %d'
          % (G['gid'], len(st), len(seen_orbits)))
print('   CM types examined (64 per class):', stats['pairs'])
print('   primitive types:', stats['prim'], '; min s over primitive:', stats['min_s_prim'],
      '; max d over primitive:', stats['max_d_prim'])
print('   Lemma 4.2: 2g <= 2^s for all primitive types; functionals +-x_i coincide for a primitive type:',
      stats['func_coincide_prim'], '; types where they coincide (all imprimitive):', stats['func_coincide_imprim_ok'])
print('   primitive types with d = 1:', stats['prim_d1'], 'in', stats['prim_d1_orbits'], 'G-orbits (stored: 21)')
check(stats['prim_d1_orbits'] == 21, 'expected 21 G-orbits of primitive types with d = 1')
print('   Lemma 4.3: sextic CM subfields (systems of 6 blocks of size 2, rho free):', stats['sextic_systems'],
      '; types with one element in each fibre:', stats['sextic_balanced_types'], ', of which primitive:',
      stats['sextic_balanced_prim'])
print('   parity: types balanced for two imaginary quadratic subfields:', stats['two_balanced_K'])


# ================================================================ B. E Galois with group S_3 x C_2
print()
print('B. Lemma 4.4: G = S_3 x C_2 acting on itself (regular), iota = (1, -1)')
S3 = sorted(itertools.permutations(range(3)))          # 6 elements, S3[0] = identity
idx = {s_: j for j, s_ in enumerate(S3)}


def s3mul(a, b):
    return tuple(a[b[i]] for i in range(3))


def lmul(h):                                            # left multiplication by h = (sigma, e) on the 12 points
    sig, e = h
    p = [0] * 12
    for j, t in enumerate(S3):
        for b in (0, 1):                                # point j + 6b = (t, (-1)^b)
            u = s3mul(sig, t)
            bb = b ^ (1 if e == -1 else 0)
            p[j + 6 * b] = idx[u] + 6 * bb
    return tuple(p)


gensB = [lmul(((1, 0, 2), 1)), lmul(((1, 2, 0), 1)), lmul(((0, 1, 2), -1))]
GB = closure(gensB, 12)
check(len(GB) == 12 and is_transitive(gensB, 12), 'S3 x C2 regular: order/transitivity')
check(lmul(((0, 1, 2), -1)) == tuple((i + 6) % 12 for i in range(12)), 'iota is rho')
sysB = block_systems(gensB, 12)


def std(t):                                             # S_3 on the plane x0 + x1 + x2 = 0, basis f1 = e0 - e1, f2 = e1 - e2
    def img(v):
        w = [0, 0, 0]
        for i in range(3):
            w[t[i]] += v[i]
        return w
    cols = []
    for v in ([1, -1, 0], [0, 1, -1]):
        w = img(v)                                      # w = a f1 + b f2 with a = w0, b = w0 + w1
        cols.append((w[0], w[0] + w[1]))
    return ((cols[0][0], cols[1][0]), (cols[0][1], cols[1][1]))


sgn = {t: (1 if sum(1 for i in range(3) for j in range(i + 1, 3) if t[i] > t[j]) % 2 == 0 else -1) for t in S3}
cntB = dict(types=0, prim=0, rank_le1=0, rank_le1_prim=0, formula_ok=0, prim_bal=0)
for Phi in cm_types(6):
    cntB['types'] += 1
    s, rows, prim = analyse(gensB, 6, Phi, sysB)
    d = 6 - s
    eps = {t: (1 if idx[t] in Phi else -1) for t in S3}          # Phi contains (t, eps_t)
    A = [[sum(eps[t] * std(t)[r][c_] for t in S3) for c_ in range(2)] for r in range(2)]
    rkA = rank(A)
    chi1 = sum(eps.values())                                     # triv x sgn: K-signature defect
    chi2 = sum(sgn[t] * eps[t] for t in S3)                      # sgn x sgn: K'-signature defect
    dform = (chi1 == 0) + (chi2 == 0) + 2 * (2 - rkA)
    cntB['formula_ok'] += (dform == d)
    check(dform == d, 'Lemma 4.4: d formula fails for Phi %s' % sorted(Phi))
    if rkA <= 1:
        cntB['rank_le1'] += 1
        cntB['rank_le1_prim'] += prim
        check(not prim, 'Lemma 4.4: rank A <= 1 but Phi primitive')
    if prim:
        cntB['prim'] += 1
        check(d == 0, 'Lemma 4.4: primitive Phi with d > 0')
        cntB['prim_bal'] += (chi1 == 0 or chi2 == 0)
print('   CM types: %(types)d ; d-formula holds for %(formula_ok)d ; rank A <= 1 for %(rank_le1)d (primitive among them: '
      '%(rank_le1_prim)d) ; primitive types: %(prim)d, all with d = 0 ; primitive and balanced for K or K\': %(prim_bal)d'
      % cntB)


# ================================================================ C. transitive groups of degree 6 (Lemma 4.5, Table 1)
print()
print('C. Lemma 4.5 and Table 1: the 16 transitive groups H_0 of degree 6 (built from their definitions)')


def perm_from_map(f, pts):
    ix = {p: i for i, p in enumerate(pts)}
    return tuple(ix[f(p)] for p in pts)


def on_points(pts, fs):
    return [perm_from_map(f, pts) for f in fs]


def s_n_action(gens_tuples, pts, act):
    return [perm_from_map(lambda p, gt=gt: act(gt, p), pts) for gt in gens_tuples]


Z6 = list(range(6))
# 6T1 C6, 6T3 D6 (hexagon)
H = {}
H[1] = ('C6', on_points(Z6, [lambda i: (i + 1) % 6]))
# 6T2 S3 acting on itself
H[2] = ('S3 (regular)', [perm_from_map(lambda t, a=a: s3mul(a, t), S3) for a in [(1, 0, 2), (1, 2, 0)]])
H[3] = ('D6 (hexagon)', on_points(Z6, [lambda i: (i + 1) % 6, lambda i: (-i) % 6]))
# 6T4 A4 and 6T7 S4 on the 2-subsets of {0,1,2,3}
pairs4 = [frozenset(p) for p in itertools.combinations(range(4), 2)]
act2 = lambda gt, p: frozenset(gt[i] for i in p)
H[4] = ('A4 on 2-subsets of {0..3}', s_n_action([(1, 2, 0, 3), (1, 0, 3, 2)], pairs4, act2))
H[7] = ('S4 on 2-subsets of {0..3}', s_n_action([(1, 0, 2, 3), (1, 2, 3, 0)], pairs4, act2))
# 6T8 S4 on the 6 left cosets of C4 = <(0 1 2 3)>
S4 = list(itertools.permutations(range(4)))
c4 = (1, 2, 3, 0)
C4 = [tuple(range(4)), c4, mul(c4, c4), mul(c4, mul(c4, c4))]
cosets = sorted({frozenset(mul(a, z) for z in C4) for a in S4}, key=lambda s_: sorted(s_))
H[8] = ('S4 on cosets of C4', [perm_from_map(lambda C, a=a: frozenset(mul(a, z) for z in C), cosets)
                               for a in [(1, 0, 2, 3), (1, 2, 3, 0)]])
# wreath-type groups on two blocks {0,1,2}, {3,4,5} (point 3b + x) or three blocks {0,1}, {2,3}, {4,5} (point 2b + e)
P2 = [(b, x) for b in range(2) for x in range(3)]
c0 = lambda p: (p[0], (p[1] + 1) % 3) if p[0] == 0 else p
c1 = lambda p: (p[0], (p[1] + 1) % 3) if p[0] == 1 else p
t0 = lambda p: (p[0], {0: 0, 1: 2, 2: 1}[p[1]]) if p[0] == 0 else p
tt = lambda p: (p[0], {0: 0, 1: 2, 2: 1}[p[1]])
sw = lambda p: (1 - p[0], p[1])
t0sw = lambda p: t0(sw(p))
H[5] = ('C3 wr C2 (= C3 x S3)', on_points(P2, [c0, sw]))
H[9] = ('S3 x S3 (index 2 in S3 wr C2, contains the block swap)', on_points(P2, [c0, c1, tt, sw]))
H[10] = ('3^2:4 (index 2 in S3 wr C2)', on_points(P2, [c0, c1, t0sw]))
H[13] = ('S3 wr C2', on_points(P2, [c0, t0, sw]))
P3 = [(b, e) for b in range(3) for e in range(2)]
f0 = lambda p: (p[0], 1 - p[1]) if p[0] == 0 else p
r3 = lambda p: ((p[0] + 1) % 3, p[1])
s01 = lambda p: ({0: 1, 1: 0, 2: 2}[p[0]], p[1])
H[6] = ('C2 wr C3 (= C2 x A4)', on_points(P3, [f0, r3]))
H[11] = ('C2 wr S3 (= C2 x S4)', on_points(P3, [f0, r3, s01]))
# 6T12 PSL(2,5), 6T14 PGL(2,5) on P^1(F_5) = {0,...,4, oo}
P1 = [0, 1, 2, 3, 4, 'oo']
tr = lambda x: 'oo' if x == 'oo' else (x + 1) % 5
sq = lambda x: 'oo' if x == 'oo' else (4 * x) % 5
inv = lambda x: 0 if x == 'oo' else ('oo' if x == 0 else (-pow(x, 3, 5)) % 5)     # x -> -1/x
db = lambda x: 'oo' if x == 'oo' else (2 * x) % 5
H[12] = ('PSL(2,5) = A5 on P^1(F_5)', on_points(P1, [tr, sq, inv]))
H[14] = ('PGL(2,5) = S5 on P^1(F_5)', on_points(P1, [tr, db, inv]))
H[15] = ('A6', [tuple({0: 1, 1: 2, 2: 0}.get(i, i) for i in range(6))]
         + [tuple({0: 1, 1: k, k: 0}.get(i, i) for i in range(6)) for k in range(3, 6)])
H[16] = ('S6', [tuple((i + 1) % 6 for i in range(6)), (1, 0, 2, 3, 4, 5)])

orders = {1: 6, 2: 6, 3: 12, 4: 12, 5: 18, 6: 24, 7: 24, 8: 24, 9: 36, 10: 36, 11: 48, 12: 60, 13: 72, 14: 120,
          15: 360, 16: 720}
combi = json.load(open(os.path.join(REPO, 'examples', 'combi_rows.json')))


def even(p):
    seen, par = set(), 0
    for i in range(len(p)):
        if i in seen:
            continue
        L, j = 0, i
        while j not in seen:
            seen.add(j)
            j = p[j]
            L += 1
        par += L - 1
    return par % 2 == 0


tot_orbits = tot_classes = 0
rows_out = []
for j in range(1, 17):
    name, gens6 = H[j]
    els = closure(gens6, 6)
    check(len(els) == orders[j] and is_transitive(gens6, 6), '6T%d: order %d (expected %d)' % (j, len(els), orders[j]))
    sys6 = block_systems(gens6, 6)
    b2 = [s_ for s_ in sys6 if len(s_) == 3]
    b3 = [s_ for s_ in sys6 if len(s_) == 2]
    in_A6 = all(even(x) for x in gens6)
    # standard tables: the transitive subgroups of A6 of degree 6 are 6T4, 6T7, 6T10, 6T12 and 6T15
    check(in_A6 == (j in (4, 7, 10, 12, 15)), '6T%d: contained in A6 is %s, unlike the standard tables' % (j, in_A6))
    # G = H_0 x <iota> on the 12 points (i, +) = i, (i, -) = i + 6
    gensG = [tuple(list(x) + [6 + x[i] for i in range(6)]) for x in gens6] + [tuple((i + 6) % 12 for i in range(12))]
    sysG = block_systems(gensG, 12)
    triples = [frozenset(P) for P in itertools.combinations(range(6), 3)]
    crit_prim, direct_prim, d1 = [], [], 0
    for P in triples:
        is_block3 = any(P in s_ for s_ in b3)
        is_transv = any(all(len(P & B) == 1 for B in s_) for s_ in b2)
        crit = not (is_block3 or is_transv)
        Phi = frozenset(list(P) + [i + 6 for i in range(6) if i not in P])
        s, rows, prim = analyse(gensG, 6, Phi, sysG)
        check(crit == prim, '6T%d: Lemma 4.5 criterion disagrees for P = %s' % (j, sorted(P)))
        if prim:
            direct_prim.append(P)
            check(6 - s == 1, '6T%d: primitive balanced type with d = %d' % (j, 6 - s))
    # orbits of the primitive P under H_0 and complementation
    comp = lambda P: frozenset(range(6)) - P
    rem, norb = set(direct_prim), 0
    while rem:
        P = rem.pop()
        orb = {frozenset(x[i] for i in P) for x in els}
        orb |= {comp(Q) for Q in orb}
        rem -= orb
        norb += 1
    # primitive types of G with d = 1 for the other imaginary quadratic subfield (if any), up to G
    Ks = iq_blocks(sysG, 6)
    other = 0
    if len(Ks) > 1:
        seen = set()
        for Phi in cm_types(6):
            s, rows, prim = analyse(gensG, 6, Phi, sysG)
            KK = [K for K in Ks if len(K & Phi) == 3]
            if prim and 6 - s == 1 and KK and KK[0] != frozenset(range(6)) and frozenset(range(6)) not in KK:
                key = frozenset(set_orbit(gensG, Phi))
                if key not in seen:
                    seen.add(key)
                    other += 1
    n = len(direct_prim)
    cr = combi[str(j)]['summary']
    check(n == cr['n_prim_balanced'], '6T%d: n = %d vs stored n_prim_balanced %d' % (j, n, cr['n_prim_balanced']))
    nprim_all = sum(1 for Phi in cm_types(6) if analyse(gensG, 6, Phi, sysG)[2])
    check(nprim_all == cr['n_prim'], '6T%d: primitive types %d vs stored n_prim %d' % (j, nprim_all, cr['n_prim']))
    if n:
        tot_orbits += norb + other
        tot_classes += 1
    rows_out.append((j, name, len(els), len(b2), len(b3), in_A6, n, norb, other, cr['TI'], len(Ks)))
print('   6Tj  |H_0|  #sys(3x2)  #sys(2x3)  H_0<A6  n  orbits(K)  orbits(K\')  #IQ  stored 12T  name')
for r in rows_out:
    print('   %-4s %5d  %9d  %9d  %6s %3d  %9d  %10d  %3d  %-10s %s'
          % ('6T%d' % r[0], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[10], r[9], r[1]))
print('   groups with primitive types of signature (3,3): %d ; G-orbits of primitive types with d = 1: %d'
      % (tot_classes, tot_orbits))
check(tot_classes == 15 and tot_orbits == 21, 'expected 15 groups and 21 orbits')
# 6T7 versus 6T8 in the stored enumeration: project the stored generators of 12T23 / 12T24 to the 6 conjugate pairs
for gid, lab in ((877, '12T23'), (879, '12T24')):
    G = next(G for G in D['groups'] if G['gid'] == gid)
    proj = [tuple(x[i] % 6 for i in range(6)) for x in G['gens']]
    print('   stored gid %d (%s): image on the 6 conjugate pairs has order %d, contained in A6: %s'
          % (gid, lab, len(closure(proj, 6)), all(even(x) for x in proj)))
    check(len(closure(proj, 6)) == 24 and all(even(x) for x in proj) == (lab == '12T23'),
          'gid %d: the image on the 6 pairs is not %s' % (gid, '6T7' if lab == '12T23' else '6T8'))

print()
print('RESULT:', 'all checks passed' if not FAILS else '%d FAILURES' % len(FAILS))
