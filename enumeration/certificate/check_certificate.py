# check_certificate.py -- pure Python 3, standard library only; reads only certificate.json.
# Verifies, for every record (a reduced pair (G, Phi) with d >= 1 and g <= 6, up to G):
#   0. the entries are integers of the right shape: gens, orbits, Phi (strictly increasing, in [0, 2g)), d, LU, and
#      for "kinds a-d" the signs chi (+1 or -1), the multiplicities c_i (non-zero) and T;
#   1. G = <gens> lies in W(B_g) (commutes with rho(i) = i + g mod 2g) and contains rho; the stored orbits are the
#      G-orbits; Phi is a CM type;
#   2. the pair is reduced: on each orbit the CM type is primitive (no block system of G on the orbit, with rho
#      permuting the blocks without fixed points, has Phi as a union of blocks), and no G-set isomorphism between two
#      orbits maps Phi onto Phi;
#   3. d = g - rank{mu_h : h in G}; the rows of LU are orthogonal to every mu_h and form a saturated lattice of rank d,
#      so they are a basis of Lambda_U;
#   4. "kinds a-d": each listed character is the Weil character of a balanced sub-product B of powers of A of one of
#      the kinds (a) dim B = 2, (b) dim B = 4, (c) dim B = 6 with a factor of odd dimension, (d) B = a simple factor
#      of dimension 6: the signs chi define a homomorphism G -> {+1,-1} with chi(rho) = -1, the field K it defines
#      lies in each factor used, and sum c_i t_i = 0; and LU = T * (characters) with T an integer matrix.  Hence
#      Lambda_W = Lambda_U.
#   5. "table2": the Weil characters of imaginary quadratic fields on balanced sub-products of powers of A of
#      dimension <= 6 (all fields, all multiplicities), which are checked to lie in Lambda_U, do not generate Lambda_U.
#      For these records the program also prints the rank of the lattice generated with dimension <= 12 and whether
#      it equals Lambda_U (Table 2, "IQ").
# and, over all records:
#   6. no two records with the same group G have CM types in the same G-orbit;
#   7. for the "table2" records (by family) and for the transitive records with g = 6: the classes up to conjugacy in
#      W(B_g), computed by brute force over W(B_g) ((G_1, Phi_1) and (G_2, Phi_2) are conjugate if some w in W(B_g)
#      conjugates G_1 onto G_2 and maps Phi_1 into the G_2-orbit of Phi_2), are the classes given by the stored labels
#      (gid, Nclass); the program prints their numbers.
# That the records include every reduced pair with d >= 1 rests on the enumeration and cannot be checked from the file.
#   usage:  python check_certificate.py [certificate.json] > check_certificate_out.txt
import json, os, sys, itertools
from fractions import Fraction


def mul(a, b):
    return tuple(a[b[i]] for i in range(len(b)))


def is_int(x):
    return type(x) is int


def int_matrix(M, ncols):
    return isinstance(M, list) and all(isinstance(row, list) and len(row) == ncols and all(is_int(x) for x in row)
                                       for row in M)


def elements_with(gens, n, signs=None):
    """BFS over G = <gens>; returns dict element -> value of the sign character (or 1), or None if the signs
    do not define a homomorphism."""
    e = tuple(range(n))
    val = {e: 1}
    frontier = [e]
    while frontier:
        new = []
        for x in frontier:
            for k, s in enumerate(gens):
                y = mul(s, x)
                v = val[x] * (signs[k] if signs else 1)
                if y in val:
                    if val[y] != v:
                        return None
                else:
                    val[y] = v
                    new.append(y)
        frontier = new
    return val


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


def rank(rows):
    M = [[Fraction(x) for x in r] for r in rows]
    r = 0
    for c in range(len(M[0]) if M else 0):
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


def hnf(rows):
    """Row Hermite normal form of the integer lattice spanned by rows (zero rows dropped)."""
    A = [list(r) for r in rows if any(r)]
    out = []
    if not A:
        return out
    ncol = len(A[0])
    for c in range(ncol):
        nz = [r for r in A if r[c] != 0]
        A = [r for r in A if r[c] == 0]
        while len(nz) > 1:
            nz.sort(key=lambda r: abs(r[c]))
            p = nz[0]
            rest = []
            for r in nz[1:]:
                q = r[c] // p[c]
                v = [a - q * b for a, b in zip(r, p)]
                (rest if v[c] != 0 else A).append(v)
            nz = [p] + rest
        if nz:
            p = nz[0]
            if p[c] < 0:
                p = [-x for x in p]
            for i, o in enumerate(out):
                q = o[c] // p[c]
                out[i] = [a - q * b for a, b in zip(o, p)]
            out.append(p)
            A = [r for r in A if any(r)]
    return out


def det(M):
    M = [[Fraction(x) for x in r] for r in M]
    n, s = len(M), Fraction(1)
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c] != 0), None)
        if p is None:
            return 0
        if p != c:
            M[c], M[p] = M[p], M[c]
            s = -s
        s *= M[c][c]
        for i in range(c + 1, n):
            f = M[i][c] / M[c][c]
            M[i] = [a - f * b for a, b in zip(M[i], M[c])]
    return int(s)


def gcd_minors(B):
    from math import gcd
    d, g = len(B), len(B[0])
    out = 0
    for cols in itertools.combinations(range(g), d):
        out = gcd(out, abs(det([[r[c] for c in cols] for r in B])))
    return out


def blocks_on_orbit(gens, O):
    """Non-trivial block systems of G on the orbit O (O is a G-orbit), as lists of frozensets."""
    O = sorted(O)
    n, x0, out = len(O), O[0], []
    for k in range(2, n):
        if n % k:
            continue
        for rest in itertools.combinations(O[1:], k - 1):
            B = frozenset((x0,) + rest)
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
                out.append(list(orb))
    return out


def check_record(r):
    errs = []
    g = r.get('g')
    if not is_int(g) or not 1 <= g <= 8:
        return ['g is not an integer between 1 and 8']
    n = 2 * g
    if not r['gens'] or not int_matrix(r['gens'], n):
        return ['gens are not integer lists of length 2g']
    if not isinstance(r['orbits'], list) or not all(isinstance(o, list) and all(is_int(x) for x in o)
                                                    for o in r['orbits']):
        return ['orbits are not lists of integers']
    if not (isinstance(r['Phi'], list) and all(is_int(x) and 0 <= x < n for x in r['Phi'])
            and r['Phi'] == sorted(set(r['Phi']))):
        return ['Phi is not a strictly increasing list of integers in [0, 2g)']
    if not is_int(r['d']) or not int_matrix(r['LU'], g):
        return ['d or LU is not integral']
    gens = [tuple(x) for x in r['gens']]
    rho = tuple((i + g) % n for i in range(n))
    for x in gens:
        if sorted(x) != list(range(n)) or mul(x, rho) != mul(rho, x):
            errs.append('generator not in W(B_g)')
            return errs
    els = elements_with(gens, n)
    if rho not in els:
        errs.append('rho not in G')
        return errs
    orbs = sorted((sorted(o) for o in {frozenset(e[i] for e in els) for i in range(n)}), key=min)
    if orbs != sorted((sorted(o) for o in r['orbits']), key=min):
        errs.append('orbits differ')
        return errs
    Phi = frozenset(r['Phi'])
    if len(Phi) != g or any(((i + g) % n) in Phi for i in Phi):
        errs.append('Phi is not a CM type')
        return errs
    # reduced: primitive on each orbit, no isomorphic orbits with CM types
    for O in orbs:
        for sysm in blocks_on_orbit(gens, O):
            B0 = sysm[0]
            if frozenset((i + g) % n for i in B0) == B0:
                continue
            if all(B <= Phi or not (B & Phi) for B in sysm):
                errs.append('Phi imprimitive on an orbit')
                break
    for Oa, Ob in itertools.combinations(orbs, 2):
        if len(Oa) != len(Ob):
            continue
        a0 = Oa[0]
        stab = [x for x in els if x[a0] == a0]
        for b in Ob:
            if all(x[b] == b for x in stab):          # a0 -> b extends to a G-set isomorphism Oa -> Ob
                f = {}
                for x in els:
                    f[x[a0]] = x[b]
                if all((f[i] in Phi) == (i in Phi) for i in Oa):
                    errs.append('two orbits with isomorphic CM types (not reduced)')
                    break
    # d and Lambda_U
    mus = [[1 if i in P else -1 for i in range(g)] for P in set_orbit(gens, Phi)]
    d = g - rank(mus)
    if d != r['d'] or d < 1:
        errs.append('d differs: %d' % d)
    LU = r['LU']
    if len(LU) != d or rank(LU) != d or any(sum(a * m for a, m in zip(v, mu)) for v in LU for mu in mus):
        errs.append('LU is not in Lambda_U or has the wrong rank')
    elif gcd_minors(LU) != 1:
        errs.append('LU not saturated')
    if r['status'] == 'kinds a-d':
        W = []
        for wc in r['weil']:
            chi = wc['chi']
            if not (isinstance(chi, list) and len(chi) == len(gens) and all(is_int(s) and s in (1, -1) for s in chi)):
                errs.append('chi is not a list of signs, one for each generator')
                break
            val = elements_with(gens, n, chi)
            if val is None or val[rho] != -1:
                errs.append('chi is not a homomorphism with chi(rho) = -1')
                break
            w, t, D, odd, copies = [0] * g, 0, 0, False, 0
            for o_s, c in wc['c'].items():
                if not (o_s.isdigit() and int(o_s) < len(orbs)) or not is_int(c) or c == 0:
                    errs.append('multiplicities must be non-zero integers indexed by orbits')
                    break
                O = orbs[int(o_s)]
                B = frozenset(x[O[0]] for x, v in val.items() if v == 1)
                if len(B) * 2 != len(O):
                    errs.append('K is not contained in the factor of orbit %s' % o_s)
                    break
                for x in B:
                    if x < g:
                        w[x] += c
                    else:
                        w[x - g] -= c
                t += c * (2 * len(B & Phi) - len(B))
                D += abs(c) * len(O) // 2
                odd = odd or (len(O) // 2) % 2 == 1
                copies += abs(c)
            if t != 0:
                errs.append('sub-product not balanced')
            kind = ('a' if D == 2 else 'b' if D == 4 else 'c' if D == 6 and odd else
                    'd' if D == 6 and copies == 1 else None)
            if kind is None or kind != wc['kind']:
                errs.append('kind %s does not match (D=%d)' % (wc['kind'], D))
            if any(sum(a * m for a, m in zip(w, mu)) for mu in mus):
                errs.append('Weil character not orthogonal to the Hodge cocharacters')
            W.append(w)
        T = r['T']
        if len(T) != len(LU) or not int_matrix(T, len(W)):
            errs.append('T is not an integer matrix of the right shape')
        elif any([sum(T[i][k] * W[k][j] for k in range(len(W))) for j in range(g)] != LU[i] for i in range(len(LU))):
            errs.append('LU != T * W')
    elif r['status'] == 'table2':
        res = iq_lattice(r, gens, els, orbs, Phi, 6)
        if any(sum(a * m for a, m in zip(w, mu)) for w in res for mu in mus):
            errs.append('an imaginary-quadratic Weil character is not in Lambda_U')
        elif hnf(res) == hnf(LU):
            errs.append('table2 record is generated at dimension <= 6')
    else:
        errs.append('unknown status')
    return errs


def iq_lattice(r, gens, els, orbs, Phi, DMAX):
    g, n = r['g'], 2 * r['g']
    rho = tuple((i + g) % n for i in range(n))
    out = []
    for signs in itertools.product((1, -1), repeat=len(gens)):
        val = elements_with(gens, n, list(signs))
        if val is None or val[rho] != -1:
            continue
        blocks = []
        for oi, O in enumerate(orbs):
            B = frozenset(x[O[0]] for x, v in val.items() if v == 1)
            if 2 * len(B) == len(O):
                blocks.append((O, B))
        ranges = [range(-(DMAX // (len(O) // 2)), DMAX // (len(O) // 2) + 1) for O, B in blocks]
        for cv in itertools.product(*ranges):
            D = sum(abs(c) * len(O) // 2 for c, (O, B) in zip(cv, blocks))
            if D == 0 or D > DMAX:
                continue
            if sum(c * (2 * len(B & Phi) - len(B)) for c, (O, B) in zip(cv, blocks)) != 0:
                continue
            w = [0] * g
            for c, (O, B) in zip(cv, blocks):
                for x in B:
                    if x < g:
                        w[x] += c
                    else:
                        w[x - g] -= c
            if any(w):
                out.append(w)
    return out


# ---- checks over all records: duplicates, and classes up to conjugacy in W(B_g) by brute force ----

_W, _CONJ = {}, {}


def wbg(g):
    """The elements of W(B_g), as pairs (w, w^-1) of permutations of {0, ..., 2g-1} that commute with rho."""
    if g not in _W:
        out = []
        for p in itertools.permutations(range(g)):
            for s in itertools.product((0, g), repeat=g):
                w = [0] * (2 * g)
                for i in range(g):
                    w[i], w[i + g] = p[i] + s[i], p[i] + g - s[i]
                wi = [0] * (2 * g)
                for i, x in enumerate(w):
                    wi[x] = i
                out.append((tuple(w), tuple(wi)))
        _W[g] = out
    return _W[g]


def prepare(r):
    gens = [tuple(x) for x in r['gens']]
    els = elements_with(gens, 2 * r['g'])
    return dict(g=r['g'], gens=gens, els=els, key=frozenset(els), orbs=[frozenset(o) for o in r['orbits']],
                Phi=frozenset(r['Phi']), Phiorb=set_orbit(gens, r['Phi']), label=(r['gid'], r['Nclass']))


def conjugators(A, B):
    """All w in W(B_g) with w G_A w^-1 = G_B."""
    k = (A['key'], B['key'])
    if k not in _CONJ:
        out = []
        if A['g'] == B['g'] and len(A['els']) == len(B['els']):
            orbsB, n = set(B['orbs']), 2 * A['g']
            for w, wi in wbg(A['g']):
                if set(frozenset(w[x] for x in O) for O in A['orbs']) != orbsB:
                    continue
                if all(tuple(w[s[wi[i]]] for i in range(n)) in B['els'] for s in A['gens']):
                    out.append(w)
        _CONJ[k] = out
    return _CONJ[k]


def partition(items, related):
    parent = list(range(len(items)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, j in itertools.combinations(range(len(items)), 2):
        if find(i) != find(j) and related(items[i], items[j]):
            parent[find(i)] = find(j)
    out = {}
    for i in range(len(items)):
        out.setdefault(find(i), set()).add(i)
    return set(frozenset(v) for v in out.values())


def conjugate_pairs(A, B):
    if sorted(map(len, A['orbs'])) != sorted(map(len, B['orbs'])):
        return False
    return any(frozenset(w[x] for x in A['Phi']) in B['Phiorb'] for w in conjugators(A, B))


def global_checks(recs):
    """Returns (list of errors, {name: (number of classes, number of records)})."""
    errs, counts = [], {}
    items = [prepare(r) for r in recs]
    for i, j in itertools.combinations(range(len(items)), 2):
        A, B = items[i], items[j]
        if A['g'] == B['g'] and A['key'] == B['key'] and B['Phi'] in A['Phiorb']:
            errs.append('records %d and %d: same group and CM types in the same G-orbit' % (i, j))
    subsets = {}
    for it, r in zip(items, recs):
        if r['status'] == 'table2':
            subsets.setdefault(('table2', tuple(sorted((len(o) // 2 for o in r['orbits']), reverse=True))),
                               []).append(it)
        elif r['g'] == 6 and len(r['orbits']) == 1:
            subsets.setdefault(('transitive g=6', (6,)), []).append(it)
    for name, sub in subsets.items():
        brute = partition(sub, conjugate_pairs)
        labels = partition(sub, lambda A, B: A['label'] == B['label'])
        if brute != labels:
            errs.append('%s %s: the stored labels (gid, Nclass) do not give the classes up to W(B_g)' % name)
        counts[name] = (len(brute), len(sub))
    return errs, counts


def main(path):
    data = json.load(open(path))
    recs = data['records']
    bad = 0
    summary = {}
    t2 = {}
    for r in recs:
        errs = check_record(r)
        if errs:
            bad += 1
            print('ERROR g=%s gid=%s Phi=%s: %s' % (r.get('g'), r.get('gid'), r.get('Phi'), '; '.join(errs)))
        key = (r['g'], r['status'])
        summary[key] = summary.get(key, 0) + 1
        if r['status'] == 'table2' and not errs:
            gens = [tuple(x) for x in r['gens']]
            els = elements_with(gens, 2 * r['g'])
            orbs = sorted((sorted(o) for o in r['orbits']), key=min)
            L12 = hnf(iq_lattice(r, gens, els, orbs, frozenset(r['Phi']), 12))
            part = tuple(sorted((len(o) // 2 for o in orbs), reverse=True))
            e = t2.setdefault(part, dict(ranks=set(), equal=set()))
            e['ranks'].add(len(L12))
            e['equal'].add(L12 == hnf(r['LU']))
    print('records:', len(recs), '; with errors:', bad)
    for k in sorted(summary):
        print('  g=%d %-10s %d' % (k[0], k[1], summary[k]))
    gerrs, counts = global_checks(recs) if bad == 0 else (['not run (record errors)'], {})
    for e in gerrs:
        print('ERROR', e)
    print('Table 2 (factor dimensions: classes up to W(B_6) (G-orbits); IQ rank at dimension <= 12; equal to '
          'Lambda_U):')
    for part in sorted(t2, key=lambda p: (-p[0], p)):
        e = t2[part]
        c = counts.get(('table2', part), (None, None))
        print('  %-10s %s (%s) ; rank %s ; equal %s' % (part, c[0], c[1], sorted(e['ranks']), sorted(e['equal'])))
    c = counts.get(('transitive g=6', (6,)), (None, None))
    print('transitive pairs with g = 6 (simple sixfolds, d = 1): %s classes up to W(B_6) (%s G-orbits)' % c)
    print('no two records with the same group and CM types in the same G-orbit:', not any('same G-orbit' in e
                                                                                          for e in gerrs))
    ok = bad == 0 and not gerrs
    print('RESULT:', 'all records verified' if ok else 'FAILED (%d records, %d other errors)' % (bad, len(gerrs)))
    return 0 if ok else 1


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, 'certificate.json')))
