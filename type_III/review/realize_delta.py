"""
realize_delta.py  (review of Remark rem:III-examples: "every Delta < 0 occurs for every D")

Claim (proof in REVIEW_thmC.txt, finding M5): for every definite D and every squarefree m > 0 there is
a diagonal T = <l1, l2, l3> (l_r pure) with Delta = -Nrd(T) = -m mod squares.
Local fact used: n > 0 is the reduced norm of a pure quaternion of D iff -n is not a square in Q_p
for every finite p | disc(D)  (Hasse-Minkowski for the ternary form of pure quaternions).

This script: for the 18 discriminants of the project's all_delta.py, takes D = (a,c) (first found),
enumerates pure quaternions with coordinates in [-B, B], records the square classes of their reduced
norms, checks that every recorded class satisfies the local condition above (consistency check), and
searches for every squarefree m <= MMAX a triple of classes with product m.  Each triple found is
re-checked by computing Nrd(T) independently (crit_check.nrd_T: determinant of x -> xT on Q^12).
"""
import sys, itertools, json
from crit_check import Quat, ram_set, sqfree_part, nrd_T, is_square_Qp

DISCS = [2, 3, 5, 7, 11, 13, 17, 19, 23, 30, 42, 66, 70, 78, 102, 105, 110, 130]


def primes_of(n):
    out, p = [], 2
    while p * p <= n:
        while n % p == 0:
            out.append(p); n //= p
        p += 1
    if n > 1:
        out.append(n)
    return sorted(set(out))


def find_algebra(N):
    target = {0} | set(primes_of(N))
    for a in range(-1, -80, -1):
        for c in range(-1, -400, -1):
            if set(ram_set(a, c)) == target:
                return a, c
    raise ValueError(N)


def main():
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    MMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    ms = [m for m in range(1, MMAX + 1) if sqfree_part(m) == m]
    tot = dict(algebras=0, targets=0, realized=0, nrd_check_fail=0, local_condition_violations=0)
    out = []
    for N in DISCS:
        a, c = find_algebra(N)
        Q = Quat(a, c)
        ps = primes_of(N)
        classes = {}
        for x, y, z in itertools.product(range(-B, B + 1), repeat=3):
            if (x, y, z) == (0, 0, 0):
                continue
            lam = (0, x, y, z)
            n = Q.nrd(lam)
            assert n > 0
            s = int(sqfree_part(n))
            # local condition: -s not a square in Q_p for p | disc D
            if any(is_square_Qp(-s, p) for p in ps):
                tot["local_condition_violations"] += 1
            if s not in classes or sum(map(abs, lam)) < sum(map(abs, classes[s])):
                classes[s] = lam
        S = sorted(classes)
        missing = []
        for m in ms:
            found = None
            for s1, s2 in itertools.combinations_with_replacement(S, 2):
                s3 = int(sqfree_part(m * s1 * s2))
                if s3 in classes:
                    found = (s1, s2, s3)
                    break
            tot["targets"] += 1
            if found is None:
                missing.append(m)
                continue
            z = (0, 0, 0, 0)
            T = [[classes[found[0]], z, z], [z, classes[found[1]], z], [z, z, classes[found[2]]]]
            Nrd = nrd_T(Q, T)
            if int(sqfree_part(Nrd)) != m:
                tot["nrd_check_fail"] += 1
                missing.append(m)
                continue
            tot["realized"] += 1
        tot["algebras"] += 1
        line = (f"disc {N}: D=({a},{c}) ram={sorted(ram_set(a, c))} B={B}: {len(S)} norm classes; "
                f"Delta=-m realized for all squarefree m<={MMAX} except {missing}")
        print(line); sys.stdout.flush()
        out.append(dict(disc=N, a=a, c=c, n_classes=len(S), missing=missing))
    print("TOTAL:", json.dumps(tot))
    json.dump(out, open("realize_delta_out.json", "w"), indent=0)


if __name__ == "__main__":
    main()
