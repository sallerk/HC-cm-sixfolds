"""Main computation: explicit ordinary geometrically simple abelian sixfolds with exceptional Tate classes.
Construction (cf. APFV 2505.09589 Alg. 6.2, Prop. 6.7): E = F.K, F totally real sextic with Galois group 6Tj,
K imaginary quadratic; then Gal(E) = 6Tj x C2 = the APFV degenerate row (combi_rows.out). p splits completely
in E; primes of E above p <-> pairs (j, eta) (root of polF mod p, root of polK mod p). A CM type = choice
eta_j for each j; it is BALANCED (signature (3,3) w.r.t. K) iff exactly 3 eta_j equal the first root of polK.
a = prod_j P_(j,eta_j); h = class order; alpha generates a^h; u = alpha*tau(alpha)/p^h in O_F^{x,+};
pi = alpha/v with v*tau(v) = u if such a unit exists (q = p^h), else pi = alpha^2/u (q = p^(2h)).
Checks (exact unless stated): pi*tau(pi) = q; Weil polynomial integral with q-symmetry; |roots| = sqrt q
(numerical, PARI polroots); irreducible; ordinary (p does not divide a_6); geometric simplicity
(charpoly(pi^N) squarefree mod a large prime, N = lcm{n : phi(n) | |G|}); exceptional relation
N_{E/K}(pi)^12 = q^36 (exact, in K); angle rank = rank(M) - 1 where M[(g,eps),(j,e)] = v(pi at P_(g j, e+eps)),
with the action of Gal(F~/Q) on the p-adic roots of polF computed from nfsplitting (|G_F| <= 72), or the
full Sym/Alt (6T16/6T15); for 6T12/6T14 the rank is labelling-independent (combi_rows: all 20 balanced types
primitive with rank 6, all unbalanced non-induced types rank 7), and is computed for every labelling.
Usage: python produce.py <tag> <polF> <polK> <j> [nprimes] [pmax] [qmax]
  -> produce_runs/produce_<tag>.json / .out (next to this script)"""
import sys, json, itertools, time, os
from cypari import pari
HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, 'produce_runs')
os.makedirs(OUTDIR, exist_ok=True)
pari.allocatemem(4 * 10**9)
pari('default(new_galois_format,1)')
pari('default(realprecision,60)')
pari('read("%s")' % os.path.join(HERE, 'produce_lib.gp').replace('\\', '/'))
gp = pari

tag, polF, polK, jF = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
nprimes = int(sys.argv[5]) if len(sys.argv) > 5 else 4
pmax = int(sys.argv[6]) if len(sys.argv) > 6 else 3000
LOG = open(os.path.join(OUTDIR, f'produce_{tag}.out'), 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); LOG.write(s + '\n'); LOG.flush()

t0 = time.time()
polF, polK = gp(polF), gp(polK)
assert int(gp.polgalois(polF)[2]) == jF and int(gp.polsturm(polF)) == 6
assert gp.poldegree(polK) == 2 and int(gp.polsturm(polK)) == 0
GF_ORDER = int(gp.polgalois(polF)[0]); GORDER = 2 * GF_ORDER
R, a, b, k, tauimg = gp("buildE")(polF, polK)
red = gp.polredbest(R, 1)
R2, m = red[0], red[1]
a2 = gp.subst(a, 'x', m).lift(); b2 = gp.subst(b, 'x', m).lift()
trb = -polK.polcoef(1) / polK.polcoef(2)
# tau in new coordinates: the automorphism fixing a2 and sending b2 -> trb - b2
auts = gp.nfgaloisconj(R2)
tau = None
for s in auts:
    if gp(f'Mod(subst({a2}, x, {s}), {R2}) == Mod({a2}, {R2})') and \
       gp(f'Mod(subst({b2}, x, {s}), {R2}) == Mod({trb} - ({b2}), {R2})'):
        tau = s
assert tau is not None
R = R2; a = a2; b = b2; tauimg = tau
P(f'[{tag}] F = {polF}  (6T{jF}, |G_F| = {GF_ORDER}),  K = {polK},  E: R = {R}')
bnf = gp.bnfinit(R, 1)
nf = bnf.bnf_get_nf() if hasattr(bnf, 'bnf_get_nf') else gp('(b)->b.nf')(bnf)
cyc = gp('(b)->b.cyc')(bnf)
P(f'  disc(E) = {gp.factor(gp("(n)->n.disc")(nf))}, class group cyc = {cyc}, time {time.time()-t0:.1f}s')
# basis {a^i b^e} for relative coordinates
basisAB = gp.matconcat(gp(f'vector(12, n, my(i = (n-1)%6, e = (n-1)\\6); nfalgtobasis({gp("(n)->n")(nf)}, Mod({a}, {R})^i * Mod({b}, {R})^e))'))
NROOT = int(gp("rootsofunitybound")(GORDER))
P(f'  N (multiple of #mu(L)) = {NROOT}')
ELL = [int(gp.nextprime(10**40)), int(gp.nextprime(10**41))]

# Galois action of Gal(F~/Q) on the p-adic roots of polF (indexed by sorted residues mod p)
def gf_action(p, rF):
    rF = [int(r) for r in rF]
    if jF == 16:
        return [list(s) for s in itertools.permutations(range(6))], 'Sym(6) (6T16)'
    if jF == 15:
        def sign(s):
            return (-1) ** sum(1 for i in range(6) for j in range(i + 1, 6) if s[i] > s[j])
        return [list(s) for s in itertools.permutations(range(6)) if sign(s) == 1], 'Alt(6) (6T15)'
    if jF in (12, 14):
        return None, 'labelling-independent (combi_rows)'
    S = gp.subst(gp.nfsplitting(polF), 'x', 'y')
    assert gp.poldegree(S) == GF_ORDER
    rts = gp.nfroots(S, polF)
    assert len(rts) == 6
    th = gp.polrootspadic(S, p, 60)
    assert len(th) == GF_ORDER
    perms = []
    def cls(m):
        out = []
        for r in rts:
            v = gp.subst(r.lift(), 'y', th[m])
            res = int(gp('(v,p)->if(valuation(v,p)>=0, truncate(v)%p, -1)')(v, p))
            out.append(rF.index(res))
        return out
    c0 = cls(0)
    for mm in range(GF_ORDER):
        cm = cls(mm)
        perm = [None] * 6
        for i in range(6):
            perm[c0[i]] = cm[i]
        perms.append(perm)
    assert len(set(map(tuple, perms))) == GF_ORDER
    return perms, 'nfsplitting + p-adic roots'

def rank_M(perms, wv):
    rows = []
    for s in perms:
        for eps in (0, 1):
            rows.append([wv[s[j]][(e + eps) % 2] for j in range(6) for e in (0, 1)])
    return int(gp.matrank(gp.matrix(len(rows), 12, [x for r in rows for x in r])))

def all_labellings_rank(wv):
    # for 6T12/6T14: rank over every conjugate of the abstract group (labelling-independence check)
    G6 = [list(s) for s in json.load(open(os.path.join(HERE, 'combi_rows.json')))[str(jF)]['gens12']]
    gens = [g[:6] for g in G6[:-1]]
    # closure
    elts = {tuple(range(6))}; frontier = [tuple(range(6))]
    while frontier:
        nxt = []
        for e in frontier:
            for g in gens:
                c = tuple(g[e[i]] for i in range(6))
                if c not in elts: elts.add(c); nxt.append(c)
        frontier = nxt
    elts = [list(e) for e in elts]
    ranks = set()
    for x in itertools.permutations(range(6)):
        xi = [x.index(i) for i in range(6)]
        conj = [[x[e[xi[i]]] for i in range(6)] for e in elts]
        ranks.add(rank_M(conj, wv))
    return sorted(ranks)

results = []
pcount = 0
DISC = int(gp.poldisc(R)) * int(gp.poldisc(polF))
DK = gp.poldisc(polK)
QMAX = int(sys.argv[7]) if len(sys.argv) > 7 else 10**40
for p in gp.primes([2, pmax]):
    p = int(p)
    if pcount >= nprimes: break
    if DISC % p == 0: continue
    if int(gp.kronecker(DK, p)) != 1: continue
    if len(gp.polrootsmod(polF, p)) != 6: continue
    pp = gp("primepairs")(nf, polF, polK, a, b, p)
    if pp == 0: continue
    dec, tab, rF, rK = pp
    pcount += 1
    perms, how = gf_action(p, rF)
    P(f'  p = {p}: split completely; G_F action: {how}')
    for S in itertools.product([1, 2], repeat=6):
        nS = sum(1 for s in S if s == 1)
        if S[0] != 1: continue                       # complement = conjugate type (same Weil polynomial)
        if nS not in (2, 3): continue                # balanced (3) and one unbalanced family (2) as control
        try:
            piw, q, h, sq = gp("weilnumber")(bnf, R, tauimg, dec, tab, p, list(S), QMAX)
        except Exception as ex:
            P('   S', S, 'ERROR', ex); continue
        q = int(q)
        if q > QMAX:
            P(f'   S={S} bal={nS==3} q={p}^{(2*h if sq else h)} > QMAX: skipped'); continue
        Pw = gp.charpoly(piw)
        co = [int(Pw.polcoef(i)) for i in range(12, -1, -1)]
        sym = all(co[12 - i] == q ** (6 - i) * co[i] for i in range(7))
        absok = bool(gp('(P,q)->my(r=polroots(P)); vecmax(apply(z->abs(norm(z)/q-1), r)) < 1e-20')(Pw, q))
        irr = bool(gp.polisirreducible(Pw))
        ordinary = co[6] % p != 0
        gs = [bool(gp("geomsimple_modl")(Pw, NROOT, l)) for l in ELL]
        geom = any(gs)
        # exceptional relation N_{E/K}(pi) / q^3 root of unity
        cvec = gp.matsolve(basisAB, gp.nfalgtobasis(nf, piw))
        Nrel = gp(f'my(t = varhigher("t"), c = {cvec}, A = sum(n = 1, 12, c[n] * t^((n-1)%6) * Mod(x, {polK})^((n-1)\\6))); norm(Mod(A, subst({polF}, x, t)))')
        rel = bool(gp(f'({Nrel} / {q}^3)^12 == 1'))
        wvg = gp('(nf,pw,dec,tab,d)->vector(6,j,vector(2,e,nfeltval(nf,pw,dec[tab[j,e]])/d))')(nf, piw, dec, tab, (2 * h if sq else h))
        wv = [[int(wvg[j][e]) for e in (0, 1)] for j in range(6)]
        assert all(sorted(r) == [0, 1] for r in wv)
        if perms is not None:
            rk = rank_M(perms, wv); rkinfo = 'exact'
        else:
            rks = all_labellings_rank(wv); rk = rks[0] if len(rks) == 1 else None; rkinfo = f'all labellings {rks}'
        ar = None if rk is None else rk - 1
        rec = dict(tag=tag, jF=jF, polF=str(polF), polK=str(polK), polE=str(R), p=p, q=q, h=int(h), squared=bool(sq),
                   S=list(S), balanced=(nS == 3), weil_poly=co, q_symmetric=sym, abs_ok=absok, irreducible=irr,
                   ordinary=ordinary, geom_simple=geom, geom_simple_tests=gs, N_rootsofunity=NROOT,
                   norm_relation=rel, valuations=wv, rank_M=rk, rank_info=rkinfo, angle_rank=ar,
                   pi=str(piw.lift()), G_F_perms=perms if (perms is not None and len(perms) <= 200) else None)
        results.append(rec)
        P(f'   S={S} bal={nS==3} q={p}^{(2*h if sq else h)} sym={sym} |r|ok={absok} irr={irr} ord={ordinary} '
          f'geom={geom} normrel={rel} rankM={rk} ({rkinfo}) angle_rank={ar}')
json.dump(results, open(os.path.join(OUTDIR, f'produce_{tag}.json'), 'w'), indent=0)
P(f'done {time.time()-t0:.1f}s')
