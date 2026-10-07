"""crit_nonembed.py: negative direction of Theorem C(i) / Prop. typeIII-crit(ii).
For many definite D and DIAGONAL T chosen so that Q(sqrt Delta) does NOT embed in D, test many
imaginary quadratic K = Q(q) (q pure, |coords| <= 5) by the direct determinant computation: none may
be split.  For comparison (positive control) the same number of K is tested for a T of the same D
for which Q(sqrt Delta) embeds; there split K should be found."""
import random, json, sys
from crit_check import Quat, analyse, embeds, random_pure, ram_set, nrd_T, psi_matrix, sqfree_part


def diag_T(l1, l2, l3):
    z = (0, 0, 0, 0)
    return [[l1, z, z], [z, l2, z], [z, z, l3]]


def main():
    rnd = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 77)
    nK = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    algs = [(-1, -1), (-1, -3), (-2, -5), (-3, -5), (-1, -7), (-1, -11), (-2, -13), (-3, -17), (-1, -19),
            (-7, -10), (-3, -10), (-6, -11), (-5, -14), (-1, -30), (-13, -15)]
    tot = dict(nonembed_T=0, nonembed_K_tested=0, nonembed_K_split=0, embed_T=0, embed_K_tested=0,
               embed_K_split=0, embed_T_with_split=0, crit_mismatch=0)
    for (a, c) in algs:
        Q = Quat(a, c)
        R = sorted(ram_set(a, c))
        found = {False: None, True: None}
        for _ in range(400):
            T = diag_T(random_pure(rnd, 3), random_pure(rnd, 3), random_pure(rnd, 3))
            N = nrd_T(Q, T)
            if N == 0:
                continue
            e = embeds(-N, a, c)
            if found[e] is None:
                found[e] = T
            if found[False] is not None and found[True] is not None:
                break
        line = f"D=({a},{c}) ram={R}:"
        for e in (False, True):
            T = found[e]
            if T is None:
                line += f" [no T with embeds={e} found]"
                continue
            Psi = psi_matrix(Q, T)
            ns = 0
            for _ in range(nK):
                q = random_pure(rnd, 5)
                r, _, _, Delta = analyse(Q, T, q, Psi)
                ns += r["split"]
                tot["crit_mismatch"] += (r["split"] != r["crit"])
            key = "embed" if e else "nonembed"
            tot[f"{key}_T"] += 1; tot[f"{key}_K_tested"] += nK; tot[f"{key}_K_split"] += ns
            if e:
                tot["embed_T_with_split"] += (ns > 0)
            line += f" | embeds={e}: Delta={int(sqfree_part(Delta))}, split K {ns}/{nK}"
        print(line); sys.stdout.flush()
    print("TOTAL:", json.dumps(tot))


if __name__ == "__main__":
    main()
