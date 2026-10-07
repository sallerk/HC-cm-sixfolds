# verify_split.py -- main computation: split-polarization lemma (note, Sec. 3, Theorem A), numerical evidence.
# Lemma: E = K (x) F a CM algebra, F = prod F_i totally real of total degree 2n, K = Q(sqrt(-d)).
# For every sign vector s on the 2n real places with exactly n entries -1 (= every CM type of
# signature (n,n) w.r.t. K) there is a in F^x with sign(a) = s and
#   delta(a) := (-1)^n N_{F/Q}(a) d_F  in  Nm(K^x)      (=> the K-hermitian form Tr_{E/K}(a x ybar) is split).
# Witness certificate: rational (X,Y,Z), Z != 0, with X^2 + d Y^2 = delta Z^2 (from qfsolve).
# Also checks the "forced prime" identity: if -d is a square in F_v for every v | p then
# t_p := ((-1)^n d_F, -d)_p = 1.
import json, os, random, itertools, sys, time
from cypari import pari

HERE = os.path.dirname(os.path.abspath(__file__))   # output is written next to this script
random.seed(20261004)
pari.allocatemem(2 * 10**9)

def fields_cases():
    C = []
    cyc = lambda m, k: str(pari(f'polsubcyclo({m},{k})')) if k else None
    # F_i given as defining polynomials (strings); 'x-1' is Q
    # exact minimal polynomial of zeta_m + zeta_m^{-1}
    plus = lambda m: str(pari(f'polredabs(factor(polresultant(polcyclo({m},y), y^2 - x*y + 1, y))[1,1])'))
    q21 = plus(21); q13 = plus(13); q28 = plus(28); q36 = plus(36); q9 = plus(9); q7 = plus(7)
    q5 = plus(5); q8 = plus(8); q12 = plus(12); q11 = plus(11); q15 = plus(15); q16 = plus(16)
    for d in (3, 7):
        C.append(('Q(z21)+', [q21], d))
    for d in (1, 2, 3, 7, 11, 13):
        C.append(('Q(z13)+', [q13], d))
    for d in (1, 7, 15):
        C.append(('Q(z28)+', [q28], d))
    for d in (1, 3, 5):
        C.append(('Q(z36)+', [q36], d))
    for d in (1, 3, 23):
        C.append(('Q(z9)+ x Q(z7)+', [q9, q7], d))
    for d in (2, 5, 31):
        C.append(('Q(z5)+ x Q(z16)+', [q5, q16], d))
    for d in (1, 6, 35):
        C.append(('Q(z5)+ x Q(z8)+ x Q(z12)+', [q5, q8, q12], d))
    for d in (3, 19, 39):
        C.append(('Q x Q(z11)+', ['x-1', q11], d))
    for d in (1, 47):
        C.append(('Q x Q x Q(z15)+', ['x-1', 'x-1', q15], d))
    # random totally real sextics
    found = 0; tries = 0
    while found < 60 and tries < 400000:
        tries += 1
        co = [random.randint(-9, 9) for _ in range(6)]
        f = pari(f'x^6 + {co[0]}*x^5 + {co[1]}*x^4 + {co[2]}*x^3 + {co[3]}*x^2 + {co[4]}*x + {co[5]}')
        if co[5] == 0 or not pari.polisirreducible(f) or pari.polsturm(f) != 6:
            continue
        f = pari.polredabs(f)
        d = random.choice([1, 2, 3, 5, 6, 7, 11, 15, 23, 30, 31, 35, 39, 47, 55, 57, 59, 71])
        C.append((f'random sextic #{found}', [str(f)], d))
        found += 1
    return C


def f2_solve(rows, rhs):
    # Gaussian elimination over F2; rows: list of int bitmasks (bit j = coefficient of e_j)
    m = len(rows); piv = []
    R = [(r, b) for r, b in zip(rows, rhs)]
    sol_rows = []
    for (r, b) in R:
        for (pr, pb, pbit) in sol_rows:
            if r >> pbit & 1:
                r ^= pr; b ^= pb
        if r == 0:
            if b: return None
            continue
        pbit = r.bit_length() - 1
        # reduce existing rows
        new_rows = []
        for (pr, pb, pb2) in sol_rows:
            if pr >> pbit & 1:
                pr ^= r; pb ^= b
            new_rows.append((pr, pb, pb2))
        sol_rows = new_rows + [(r, b, pbit)]
    e = 0
    for (r, b, pbit) in sol_rows:
        if b: e |= 1 << pbit
    return e

def gens_for(nf, bnf, k, primes):
    if k == 1:
        return [pari(-1)] + [pari(p) for p in primes], ['-1'] + [str(p) for p in primes]
    S = [pr for p in primes for pr in pari.idealprimedec(nf, p)]
    su = pari.bnfsunit(bnf, S)
    G = [pari.nfalgtobasis(nf, -1)] + [pari.nfalgtobasis(nf, u) for u in pari('(b)->b.fu')(bnf)]         + [pari.nfalgtobasis(nf, u) for u in su[0]]
    return G, None

def sunit_witness(nfs, bnfs, degs, target, d, n, dF, P0):
    extra = []
    q = 2
    for rnd in range(12):
        primes = sorted(set(P0) | set(extra))
        gens = []  # (factor index, element)
        for i, (nf, bnf, k) in enumerate(zip(nfs, bnfs, degs)):
            G, _ = gens_for(nf, bnf, k, primes)
            gens += [(i, g) for g in G]
        rows_sign = [0] * sum(degs); rows_h = [0] * len(primes)
        for j, (i, g) in enumerate(gens):
            nf = nfs[i]; off = sum(degs[:i])
            sg = [int(t) for t in pari.nfeltsign(nf, g)] if degs[i] > 1 else [1 if g > 0 else -1]
            for t, sv in enumerate(sg):
                if sv == -1: rows_sign[off + t] |= 1 << j
            Ng = pari.nfeltnorm(nf, g) if degs[i] > 1 else g
            for t, p in enumerate(primes):
                if int(pari.hilbert(Ng, -d, p)) == -1: rows_h[t] |= 1 << j
        flat_target = [sv for tg in target for sv in tg]
        rhs_sign = [1 if sv == -1 else 0 for sv in flat_target]
        rhs_h = [1 if int(pari.hilbert((-1)**n * dF, -d, p)) == -1 else 0 for p in primes]
        e = f2_solve(rows_sign + rows_h, rhs_sign + rhs_h)
        if e is not None:
            parts = []
            for i, (nf, k) in enumerate(zip(nfs, degs)):
                a = pari.nfalgtobasis(nf, 1) if k > 1 else pari(1)
                for j, (ii, g) in enumerate(gens):
                    if ii == i and (e >> j) & 1:
                        a = pari.nfeltmul(nf, a, g) if k > 1 else a * g
                a = pari([a]).Col() if k == 1 else pari.nfalgtobasis(nf, a)
                parts.append(a)
            N = 1
            for nf, v, k in zip(nfs, parts, degs):
                N *= pari.nfeltnorm(nf, v) if k > 1 else v[0]
            delta = (-1)**n * N * dF
            sol = pari.qfsolve(pari.matdiagonal(pari([1, d, -delta])))
            ok_signs = all(([int(t) for t in pari.nfeltsign(nf, v)] if k > 1 else [1 if v[0] > 0 else -1]) == tg
                           for nf, v, k, tg in zip(nfs, parts, degs, target))
            if str(pari.type(sol)) == 't_COL' and sol[2] != 0 and ok_signs:
                return ({'a': [[str(c) for c in v] for v in parts],
                         'a_alg': [str(pari.nfbasistoalg(nf, v).lift()) if k > 1 else str(v[0])
                                   for nf, v, k in zip(nfs, parts, degs)],
                         'N_a': str(N), 'delta': str(delta), 'XYZ': [str(c) for c in sol]}, list(extra))
            raise RuntimeError('F2 solution did not certify: logic error')
        # enlarge S by the next prime not yet used (Chebotarev step)
        q = int(pari.nextprime(q + 1))
        while q in primes: q = int(pari.nextprime(q + 1))
        extra.append(q)
    return None, list(extra)

def run_case(label, polys, d):
    nfs = [pari.nfinit(pari(p)) for p in polys]
    bnfs = [pari.bnfinit(pari(p), 1) if int(pari.poldegree(pari(p))) > 1 else None for p in polys]
    degs = [int(pari.poldegree(nf[0])) for nf in nfs]
    two_n = sum(degs); n = two_n // 2
    dF = 1
    for nf in nfs:
        dF *= int(nf[2])
    assert dF > 0
    # forced primes
    P = set(int(p) for p in pari.factor(2 * d * dF)[0])
    forced_rows = []
    for p in sorted(P):
        forced = True
        for nf in nfs:
            for pr in pari.idealprimedec(nf, p):
                if int(pari.nfislocalpower(nf, pr, -d, 2)) != 1:
                    forced = False
        tp = int(pari.hilbert((-1)**n * dF, -d, p))
        forced_rows.append({'p': p, 'forced': forced, 't_p': tp})
    forced_ok = all(r['t_p'] == 1 for r in forced_rows if r['forced'])
    # sign patterns with exactly n negatives over the 2n real places (concatenated over factors)
    pats = [s for s in itertools.product((1, -1), repeat=two_n) if s.count(-1) == n]
    results = []
    for s in pats:
        target = []
        pos = 0
        for k in degs:
            target.append(list(s[pos:pos+k])); pos += k
        wit = None; n_signed = 0; n_split = 0; phase = 'random'
        for attempt in range(400):
            parts = []
            ok = True
            for nf, k, tg in zip(nfs, degs, target):
                # random element in the integral basis, small coefficients; rejection on signs
                v = pari([random.randint(-4, 4) for _ in range(k)]).Col()
                if all(c == 0 for c in v):
                    ok = False; break
                sg = [int(t) for t in pari.nfeltsign(nf, v)] if k > 1 else [1 if int(v[0]) > 0 else -1]
                if sg != tg:
                    ok = False; break
                parts.append(v)
            if not ok:
                continue
            n_signed += 1
            N = 1
            for nf, v in zip(nfs, parts):
                N *= pari.nfeltnorm(nf, v)
            delta = (-1)**n * N * dF
            sol = pari.qfsolve(pari.matdiagonal(pari([1, d, -delta])))
            if str(pari.type(sol)) == 't_COL' and sol[2] != 0:
                n_split += 1
                if wit is None:
                    wit = {'a': [[str(c) for c in v] for v in parts],
                           'a_alg': [str(pari.nfbasistoalg(nf, v).lift()) for nf, v in zip(nfs, parts)],
                           'N_a': str(N), 'delta': str(delta), 'XYZ': [str(c) for c in sol]}
            if wit is not None and n_signed >= 40:
                break
        if wit is None:
            # S-unit phase, mirroring the proof: a = prod g^e over generators g of the S-unit groups
            # (S = primes over 2 d dF, enlarged by auxiliary primes = the Chebotarev step). For P-units
            # the conditions are F2-linear in e: signs at the 2n real places, and
            # sum_g e_g [ (N g, -d)_p = -1 ] = [ t_p = -1 ] for p in S (other p automatic).
            phase = 'sunit'
            wit, extra = sunit_witness(nfs, bnfs, degs, target, d, n, dF, sorted(P))
            n_signed += 1
            if wit is not None:
                n_split += 1
                wit['aux_primes'] = extra
        results.append({'signs': list(s), 'witness': wit, 'n_signed_tried': n_signed, 'n_split': n_split,
                        'phase': phase})
    return {'label': label, 'polys': [str(pari(p)) for p in polys], 'degs': degs, 'd': d, 'n': n,
            'dF': str(dF), 'forced': forced_rows, 'forced_identity_ok': forced_ok, 'patterns': results,
            'all_patterns_have_witness': all(r['witness'] is not None for r in results)}

if __name__ == '__main__':
    t0 = time.time()
    out = []
    for lab, polys, d in fields_cases():
        r = run_case(lab, polys, d)
        nw = sum(x['witness'] is not None for x in r['patterns'])
        ns = sum(x['n_signed_tried'] for x in r['patterns']); nsp = sum(x['n_split'] for x in r['patterns'])
        nforced = sum(f['forced'] for f in r['forced']); nbad = sum(f['t_p'] == -1 for f in r['forced'])
        print(f"{lab:28s} d={d:3d} dF={r['dF']:>12s} witnesses {nw}/{len(r['patterns'])} "
              f"split-fraction {nsp}/{ns}  forced primes {nforced} (identity ok: {r['forced_identity_ok']})  "
              f"#p with t_p=-1: {nbad}", flush=True)
        out.append(r)
    json.dump(out, open(os.path.join(HERE, 'split_witnesses.json'), 'w'), indent=1)
    print('cases', len(out), ' all witnesses found:', all(r['all_patterns_have_witness'] for r in out),
          ' forced identity holds everywhere:', all(r['forced_identity_ok'] for r in out),
          f' time {time.time()-t0:.1f}s')
