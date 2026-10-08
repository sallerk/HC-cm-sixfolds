# mutation_test.py -- corrupts certificate.json in memory in several ways and counts how many corrupted records
# check_certificate.py rejects.  Mutations that must always be rejected are marked "must"; for the others the
# corrupted record can be valid by accident (e.g. a changed CM type in the same G-orbit), and the count is reported.
# The last two mutations concern the checks over all records (duplicates and the class labels).
#   usage:  python mutation_test.py > mutation_test_out.txt
import copy, itertools, json, os
from fractions import Fraction
from check_certificate import check_record, set_orbit, elements_with, global_checks, rank

HERE = os.path.dirname(os.path.abspath(__file__))
recs = json.load(open(os.path.join(HERE, 'certificate.json')))['records']


def m_lu(r):            # change one entry of the Lambda_U basis
    r['LU'][0][0] += 1


def m_d(r):             # wrong defect
    r['d'] += 1


def m_T(r):             # wrong transformation matrix
    if r['status'] != 'kinds a-d':
        return False
    r['T'][0][0] += 1


def m_Tfloat(r):        # same value, but not an integer
    if r['status'] != 'kinds a-d':
        return False
    r['T'][0][0] = float(r['T'][0][0])


def m_kind(r):          # wrong kind label
    if r['status'] != 'kinds a-d':
        return False
    w = r['weil'][0]
    w['kind'] = {'a': 'b', 'b': 'c', 'c': 'd', 'd': 'a'}[w['kind']]


def m_c(r):             # wrong multiplicity
    if r['status'] != 'kinds a-d':
        return False
    w = r['weil'][0]
    k = sorted(w['c'])[0]
    w['c'][k] += 1 if w['c'][k] > 0 else -1


def m_c0(r):            # an extra factor with multiplicity 0
    if r['status'] != 'kinds a-d':
        return False
    w = r['weil'][0]
    free = [str(i) for i in range(len(r['orbits'])) if str(i) not in w['c']]
    if not free:
        return False
    w['c'][free[0]] = 0


def m_status(r):        # swap the status
    if r['status'] == 'kinds a-d':
        r['status'] = 'table2'
    else:
        r['status'] = 'kinds a-d'
        r['weil'], r['T'] = [], []


def m_chi(r):           # flip the sign of one generator in the first character
    if r['status'] != 'kinds a-d':
        return False
    r['weil'][0]['chi'][0] *= -1


def m_alias(r):         # replace the largest element x of Phi by x - 2g (the same element modulo 2g)
    r['Phi'][-1] -= 2 * r['g']
    r['Phi'].sort()


def m_phi(r):           # replace one element of Phi by its conjugate
    i = r['Phi'][0]
    r['Phi'][0] = (i + r['g']) % (2 * r['g'])
    r['Phi'].sort()


def solve(LU, W):
    """Rational T with LU = T * W (W of full row rank)."""
    d = len(W)
    T = []
    for v in LU:
        # least-squares-free solve: x W = v, using Gauss-Jordan on the augmented transpose
        M = [[Fraction(W[k][j]) for k in range(d)] + [Fraction(v[j])] for j in range(len(v))]
        piv, row = [], 0
        for c in range(d):
            p = next((i for i in range(row, len(M)) if M[i][c] != 0), None)
            if p is None:
                continue
            M[row], M[p] = M[p], M[row]
            M[row] = [x / M[row][c] for x in M[row]]
            for i in range(len(M)):
                if i != row and M[i][c] != 0:
                    M[i] = [a - M[i][c] * b for a, b in zip(M[i], M[row])]
            piv.append(c)
            row += 1
        x = [Fraction(0)] * d
        for i, c in enumerate(piv):
            x[c] = M[i][d]
        if [sum(x[k] * W[k][j] for k in range(d)) for j in range(len(v))] != v:
            return None
        T.append(x)
    return T


def rational_fakes():
    """For each "table2" record whose imaginary-quadratic Weil characters of kinds a-d span Lambda_U (x) Q but do not
    generate Lambda_U: a false record with status "kinds a-d" and a non-integral T (the lattice-versus-span error)."""
    out = []
    for r in recs:
        if r['status'] != 'table2':
            continue
        g, n = r['g'], 2 * r['g']
        gens = [tuple(x) for x in r['gens']]
        rho = tuple((i + g) % n for i in range(n))
        orbs = sorted((sorted(o) for o in r['orbits']), key=min)
        Phi = frozenset(r['Phi'])
        entries = []
        for signs in itertools.product((1, -1), repeat=len(gens)):
            val = elements_with(gens, n, list(signs))
            if val is None or val[rho] != -1:
                continue
            blocks = []
            for oi, O in enumerate(orbs):
                B = frozenset(x[O[0]] for x, v in val.items() if v == 1)
                if 2 * len(B) == len(O):
                    blocks.append((oi, O, B))
            ranges = [range(-(6 // (len(O) // 2)), 6 // (len(O) // 2) + 1) for oi, O, B in blocks]
            for cv in itertools.product(*ranges):
                D = sum(abs(c) * len(O) // 2 for c, (oi, O, B) in zip(cv, blocks))
                if D == 0 or D > 6 or sum(c * (2 * len(B & Phi) - len(B)) for c, (oi, O, B) in zip(cv, blocks)):
                    continue
                odd = any(c and (len(O) // 2) % 2 for c, (oi, O, B) in zip(cv, blocks))
                copies = sum(abs(c) for c in cv)
                kind = ('a' if D == 2 else 'b' if D == 4 else 'c' if D == 6 and odd else
                        'd' if D == 6 and copies == 1 else None)
                if kind is None:
                    continue
                w = [0] * g
                for c, (oi, O, B) in zip(cv, blocks):
                    for x in B:
                        if x < g:
                            w[x] += c
                        else:
                            w[x - g] -= c
                entries.append(({'chi': list(signs), 'c': {str(oi): c for c, (oi, O, B) in zip(cv, blocks) if c},
                                 'kind': kind}, w))
        basis = []
        for e in entries:
            if rank([b[1] for b in basis] + [e[1]]) > len(basis):
                basis.append(e)
        if len(basis) < r['d']:
            continue
        T = solve(r['LU'], [b[1] for b in basis])
        if T is None or all(x.denominator == 1 for row in T for x in row):
            continue
        q = copy.deepcopy(r)
        q['status'], q['weil'], q['T'] = 'kinds a-d', [b[0] for b in basis], T
        out.append(q)
    return out


MUT = [('LU entry', m_lu, True), ('defect', m_d, True), ('T entry', m_T, True), ('T not int', m_Tfloat, True),
       ('kind label', m_kind, True), ('multiplicity', m_c, True), ('zero mult.', m_c0, True),
       ('status', m_status, True), ('Phi alias', m_alias, True), ('chi sign', m_chi, False), ('Phi', m_phi, False)]
total_must = rejected_must = 0
for name, f, must in MUT:
    tried = rej = 0
    for r in recs:
        q = copy.deepcopy(r)
        if f(q) is False:
            continue
        tried += 1
        rej += bool(check_record(q))
    print('%-13s %-5s corrupted %3d, rejected %3d' % (name, 'must' if must else '', tried, rej))
    if name == 'Phi':
        same = 0
        for r in recs:
            q = copy.deepcopy(r)
            f(q)
            if not check_record(q):
                same += frozenset(q['Phi']) in set_orbit([tuple(x) for x in r['gens']], frozenset(r['Phi']))
        print('              accepted Phi mutations: %d, of which in the G-orbit of the original CM type (the same '
              'pair): %d' % (tried - rej, same))
    if must:
        total_must += tried
        rejected_must += rej

fakes = rational_fakes()
reasons = [check_record(q) for q in fakes]
rej = sum(bool(e) for e in reasons)
only_T = sum(e == ['T is not an integer matrix of the right shape'] for e in reasons)
print('%-13s %-5s corrupted %3d, rejected %3d (only for the non-integral T: %d)' % ('rational T', 'must', len(fakes),
                                                                                 rej, only_T))
total_must += len(fakes)
rejected_must += rej

# duplicates: a second record for the same pair, with a CM type in the same G-orbit
tried = rej = 0
for r in recs:
    q = copy.deepcopy(r)
    gens = [tuple(x) for x in r['gens']]
    other = sorted(sorted(P) for P in set_orbit(gens, r['Phi']) if sorted(P) != r['Phi'])
    if not other:
        continue
    q['Phi'] = other[0]
    tried += 1
    errs, counts = global_checks([r, q])
    rej += any('same G-orbit' in e for e in errs)
print('%-13s %-5s corrupted %3d, rejected %3d' % ('duplicate', 'must', tried, rej))
total_must += tried
rejected_must += rej

# class labels: give a record the label (gid, Nclass) of another record of the same group in a different class
errs0, counts0 = global_checks(recs)
tried = rej = 0
for i, r in enumerate(recs):
    if not (r['status'] == 'table2' or (r['g'] == 6 and len(r['orbits']) == 1)):
        continue
    others = sorted({s['Nclass'] for s in recs if s['gid'] == r['gid'] and s['Nclass'] != r['Nclass']})
    if not others:
        continue
    mod = copy.deepcopy(recs)
    mod[i]['Nclass'] = others[0]
    tried += 1
    errs, counts = global_checks(mod)
    rej += any('labels' in e for e in errs)
print('%-13s %-5s corrupted %3d, rejected %3d   (unmodified certificate: %d errors)' % ('class label', 'must', tried,
                                                                                     rej, len(errs0)))
total_must += tried
rejected_must += rej
print('must-reject mutations: %d corrupted, %d rejected' % (total_must, rejected_must))
print('RESULT:', 'all must-reject mutations detected' if total_must == rejected_must else 'SOME NOT DETECTED')
