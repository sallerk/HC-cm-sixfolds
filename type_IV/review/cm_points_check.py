# cm_points_check.py -- referee check (Theorem D / Remark rem:flapan), from scratch.
# CM members of the family: M = L*C, L = Q(sqrt(-d), sqrt(m)), C a real cyclic cubic field, so
# Gal(M/Q) = Gal(L/Q) x Z/3 = {1, gs, gw, gsgw} x Z/3 (gs: sqrt m -> -sqrt m, gw = complex conjugation).
# Embeddings of M <-> group elements (after fixing one embedding phi0 restricting to sigma1 on L).
# For every CM type Phi of M whose multiplicities on L are n_sigma1 = 2, n_sigma2 = 1 (sigma1, sigma2 over the
# same embedding tau of K), compute:
#   dim Hg(A) = rank of {mu_{h Phi} : h in G} in X_*(U_M) (x) Q;
#   pairing of the character chi_tau = det_K (N_{M/K} at tau) with all mu_{h Phi}  (expected 0: Hg in ker det_K);
#   pairing of chi_sigma1 = det_L at sigma1 (N_{M/L} at sigma1) with mu_Phi       (expected n_s1 - n_s1bar = 1:
#     Hg is not contained in R_{F/Q} SU(H));
#   whether Hg = U_{M/K} (dim 5), i.e. (A,K) almost-neat (Remark rem:IV-converse), and whether Phi is primitive.
import itertools
import sympy as sp

GL = [(0, 0), (1, 0), (0, 1), (1, 1)]          # (gs exponent, gw exponent); gw = complex conjugation
G = [(l, t) for l in GL for t in range(3)]
def mul(g, h): return ((g[0][0] + h[0][0]) % 2, (g[0][1] + h[0][1]) % 2), (g[1] + h[1]) % 3
IOTA = ((0, 1), 0)
REPS = [g for g in G if g[0][1] == 0]           # one of each pair {g, g*iota}: those over tau (gw-exponent 0)
def restr(g):  # restriction to L as a label
    return {(0, 0): 's1', (1, 0): 's2', (0, 1): 's1b', (1, 1): 's2b'}[g[0]]

def mu_vec(Phi):
    return [1 if r in Phi else -1 for r in REPS]

def main():
    rows = []
    allok = True
    cnt = 0
    for choice in itertools.product([0, 1], repeat=6):
        Phi = set()
        for r, c in zip(REPS, choice):
            Phi.add(r if c == 0 else mul(r, IOTA))
        n = {lab: sum(1 for g in Phi if restr(g) == lab) for lab in ('s1', 's1b', 's2', 's2b')}
        if not (n['s1'] == 2 and n['s2'] == 1):
            continue
        cnt += 1
        vecs = [mu_vec({mul(h, g) for g in Phi}) for h in G]
        dimHg = sp.Matrix(vecs).rank()
        chi_tau = [1]*6                                   # all reps lie over tau
        chi_s1 = [1 if restr(r) == 's1' else 0 for r in REPS]
        pair_tau = [sum(a*b for a, b in zip(chi_tau, v)) for v in vecs]
        pair_s1 = sum(a*b for a, b in zip(chi_s1, mu_vec(Phi)))
        # primitivity: Phi is induced from a proper CM subfield iff its stabilizer in G (h Phi = Phi) is nontrivial
        stab = [h for h in G if {mul(h, g) for g in Phi} == Phi]
        ok = all(p == 0 for p in pair_tau) and pair_s1 == n['s1'] - n['s1b'] == 1 and dimHg <= 5
        allok &= ok
        rows.append((sorted(Phi), n, dimHg, pair_tau == [0]*12, pair_s1, len(stab)))
    print('CM types of M = L*C3 with (n_s1, n_s2) = (2, 1): %d' % cnt)
    from collections import Counter
    print('dim Hg distribution:', dict(Counter(r[2] for r in rows)))
    print('primitive (trivial stabilizer):', sum(1 for r in rows if r[5] == 1), ' of', len(rows))
    print('dim Hg = 5 (= U_{M/K}, almost-neat) among primitive:', sum(1 for r in rows if r[5] == 1 and r[2] == 5))
    print('all: <det_K, h mu> = 0 for all h, <det_L at s1, mu> = 1, dim Hg <= 5:', allok)
    for r in rows[:6]:
        print('  Phi', r[0], 'n', r[1], 'dimHg', r[2], 'detK-trivial', r[3], '<chi_s1,mu>', r[4], '|stab|', r[5])

if __name__ == '__main__':
    main()
