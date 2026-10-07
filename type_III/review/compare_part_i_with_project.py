"""compare_part_i_with_project.py: final comparison (done after the independent checks) of this review's
split decisions with the project's stored outputs for Proposition 6.2(i).
Input: ../data.json and ../path1_out.json (folder type_III/) (21168 records (D,T,q)).
For every record (or a random sample of size N), recompute with crit_check.analyse (own code):
split decision (local norm test on det phi, phi built from psi = Trd o T), Delta (sqfree part), b;
compare with the project's split_hilbert, Delta, b."""
import json, sys, random, time
from fractions import Fraction as Fr
from crit_check import Quat, analyse, psi_matrix, sqfree_part

BASE = "../"


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    data = json.load(open(BASE + "data.json"))
    recs = json.load(open(BASE + "path1_out.json"))
    algs = {(A["a"], A["c"]): A for A in data["algebras"]}
    idx = list(range(len(recs)))
    if N:
        random.Random(31).shuffle(idx)
        idx = sorted(idx[:N])
    stats = dict(compared=0, split_agree=0, split_disagree=0, delta_agree=0, b_agree=0, proj_split=0)
    cache = {}
    t0 = time.time()
    bad = []
    for n, i in enumerate(idx):
        r = recs[i]
        a, c = r["a"], r["c"]
        Q = Quat(a, c)
        T = [[tuple(int(x) for x in ent) for ent in row] for row in algs[(a, c)]["Ts"][r["ti"]]["T"]]
        key = (a, c, r["ti"])
        if key not in cache:
            cache.clear(); cache[key] = psi_matrix(Q, T)
        q = tuple(int(Fr(x)) for x in r["q"])
        res, _, _, Delta = analyse(Q, T, q, cache[key])
        stats["compared"] += 1
        stats["proj_split"] += bool(r["split_hilbert"])
        if res["split"] == bool(r["split_hilbert"]):
            stats["split_agree"] += 1
        else:
            stats["split_disagree"] += 1; bad.append(i)
        stats["delta_agree"] += (res["Delta"] == int(sqfree_part(Fr(r["Delta"]))))
        stats["b_agree"] += (res["b"] == int(sqfree_part(Fr(r["b"]))))
        if n % 2000 == 0:
            print(f"  {n}/{len(idx)} {json.dumps(stats)} ({time.time() - t0:.0f}s)"); sys.stdout.flush()
    print("FINAL:", json.dumps(stats), "first disagreements:", bad[:10])


if __name__ == "__main__":
    main()
