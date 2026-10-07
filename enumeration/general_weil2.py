# general_weil2.py -- main computation (Sage+GAP): fully general CM-field Weil characters on sub-products of powers.
#   usage (Sage):  sage -python general_weil2.py 4,5,6 12     -> general_weil2_out_D12.json in this folder
# For K' <-> M <= G (rho not in M, [G:M]=2k), a "unit" is a G-orbit Omega on O_i x G/M: it describes K' inside
# End^0(A_i^m) = M_m(E_i) (m = |Omega|/|O_i|; K' commutes with the centre E_i) -- the tau-eigenspace of
# H^1(A_i^m) over sigma has dimension #{(sigma,tau) in Omega}.  Unit contributes, for tau in G/M,
#   beta_tau = sum_{sigma:(sigma,tau) in Omega} [sigma],   s_tau = 2|Phi cap F_tau| - |F_tau|,  dimension m*g_i.
# A sub-product of a power of A with a K'-action is a non-negative combination of units; Weil type iff
# sum c*s_tau = 0 for all tau; Weil characters chi_tau = sum c*beta_tau.  ("direct" embeddings = graphs, m=1.)
import sys, os, json, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sage.all import libgap, matrix, ZZ, prod
from cmenum import *
from cmenum import _cvecs


def units_for(G, gl, ggens, orbs, M, g, rho):
    n = 2 * g
    idx = int(libgap.Index(G, M))
    cos = libgap.RightCosets(G, M)
    hom = libgap.ActionHomomorphism(G, cos, libgap.OnRight)
    imgs = [[int(libgap.OnPoints(p, libgap.Image(hom, x))) - 1 for p in range(1, idx + 1)] for x in ggens]
    rimg = [int(libgap.OnPoints(p, libgap.Image(hom, rho))) - 1 for p in range(1, idx + 1)]
    taus = []
    for t in range(idx):
        if t not in taus and rimg[t] not in taus:
            taus.append(t)
    units = []
    for i, O in enumerate(orbs):
        # orbits of G on O x G/M
        seen = set()
        oms = []
        for x in O:
            for t in range(idx):
                if (x, t) in seen:
                    continue
                om = {(x, t)}
                st = [(x, t)]
                while st:
                    a, b = st.pop()
                    for q, qc in zip(gl, imgs):
                        y = (q[a], qc[b])
                        if y not in om:
                            om.add(y)
                            st.append(y)
                seen |= om
                oms.append(frozenset(om))
        done = set()
        for om in oms:
            if om in done:
                continue
            conj = frozenset((a, rimg[b]) for a, b in om)
            done.add(om); done.add(conj)
            m = len(om) // len(O)
            beta, sig = [], []
            for t in taus:
                F = [a for a, b in om if b == t]
                w = [0] * g
                for a in F:
                    if a < g:
                        w[a] += 1
                    else:
                        w[a - g] -= 1
                beta.append(w)
                sig.append(F)
            units.append(dict(orbit=i, m=m, weight=m * len(O) // 2, beta=beta, F=sig))
    return idx, taus, units


def analyse_general2(gens, g, Phi, d, LU, Dmax=12, maxindex=24):
    n = 2 * g
    G = libgap.Group([libgap.PermList([x + 1 for x in p]) for p in gens])
    rho = rho_perm(g)
    ggens = list(libgap.GeneratorsOfGroup(G))
    gl = [perm_list(x, n) for x in ggens]
    orbs = sorted([sorted(int(x) - 1 for x in o) for o in libgap.Orbits(G, list(range(1, n + 1)))], key=min)
    Ps = set(Phi)
    allv = []
    per_field = []
    for cc in libgap.ConjugacyClassesSubgroups(G):
        M = libgap.Representative(cc)
        if bool(libgap.IsSubset(M, [rho])):
            continue
        if int(libgap.Index(G, M)) > maxindex:
            continue
        idx, taus, units = units_for(G, gl, ggens, orbs, M, g, rho)
        k = idx // 2
        for u in units:
            u['s'] = [2 * len(Ps & set(F)) - len(F) for F in u['F']]
        units = [u for u in units if any(any(b) for b in u['beta'])]   # drop units with only even characters
        if not units:
            continue
        best = None
        for c in _cvecs([u['weight'] for u in units], Dmax):
            if any(sum(ci * units[j]['s'][ti] for j, ci in enumerate(c)) for ti in range(len(taus))):
                continue
            D = sum(abs(ci) * units[j]['weight'] for j, ci in enumerate(c))
            chis = [[sum(ci * units[j]['beta'][ti][x] for j, ci in enumerate(c)) for x in range(g)] for ti in range(len(taus))]
            chis = [w for w in chis if any(w)]
            if not chis:
                continue
            if d == 0:
                raise RuntimeError("nonzero Weil character with d=0")
            assert matrix(ZZ, LU + chis).rank() == d, "Weil character not in Lambda_U"
            assert D % k == 0
            for w in chis:
                allv.append((D, idx, w))
            best = D if best is None else min(best, D)
        if best is not None:
            per_field.append((idx, str(libgap.StructureDescription(M)), best))
    allv.sort(key=lambda z: z[0])
    Dmin = 0 if d == 0 else None
    if d:
        for D in sorted(set(z[0] for z in allv)):
            Mx = matrix(ZZ, [z[2] for z in allv if z[0] <= D])
            if Mx.rank() == d and prod([int(x) for x in Mx.elementary_divisors() if x != 0]) == 1:
                Dmin = D
                break
    qr = matrix(ZZ, [z[2] for z in allv]).rank() if allv else 0
    sat = None
    if allv and qr == d:
        Mx = matrix(ZZ, [z[2] for z in allv])
        sat = prod([int(x) for x in Mx.elementary_divisors() if x != 0])
    used = sorted(set((z[1], z[0]) for z in allv if Dmin is not None and z[0] <= Dmin))
    return dict(Dmin=Dmin, qrank=qr, satindex_all=sat, fields_used=used, per_field=per_field)


if __name__ == '__main__':
    t0 = time.time()
    Dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    outall = {}
    for g in map(int, sys.argv[1].split(',')):
        D = json.load(open(os.path.join(HERE, 'enum_g%d.json' % g)))
        groups = {G['gid']: G for G in D['groups']}
        res = []
        for c in D['cases']:
            if not (c['reduced'] and c['d'] >= 1):
                continue
            G = groups[c['gid']]
            r = analyse_general2(G['gens'], g, c['Phi'], c['d'], c['LU'], Dmax)
            part = tuple(sorted((len(o) // 2 for o in G['orbits']), reverse=True))
            res.append(dict(gid=c['gid'], Phi=c['Phi'], Nclass=c['Nclass'], part=part, d=c['d'], strict=c['passed'],
                            multDmin=c['mult']['Dmin'], **r, TI=G['orbit_TI'], order=G['order']))
            print(g, c['gid'], c['Phi'], part, G['orbit_TI'], 'd=%d' % c['d'], 'strict=%s' % c['passed'],
                  'IQmultDmin=%s' % c['mult']['Dmin'], 'genDmin=%s' % r['Dmin'], 'genQrank=%d' % r['qrank'],
                  'sat(all)=%s' % r['satindex_all'], 'used(index,D)=%s' % r['fields_used'], 'perfield=%s' % r['per_field'], flush=True)
        outall[g] = res
    json.dump(outall, open(os.path.join(HERE, 'general_weil2_out_D%d.json' % Dmax), 'w'), default=int)
    print("time %.1fs" % (time.time() - t0))
