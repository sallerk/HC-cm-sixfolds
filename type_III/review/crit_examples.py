"""crit_examples.py: the hand examples of Remark rem:III-examples, recomputed with crit_check.analyse
(direct determinant of the hermitian form; no use of the criterion in the split decision)."""
from crit_check import Quat, analyse, embeds, random_pure, sqfree_part
import random

def T_diag(l1, l2, l3):
    z = (0, 0, 0, 0)
    return [[l1, z, z], [z, l2, z], [z, z, l3]]

i, j, k = (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)
print("Example 1: D=(-1,-1), T=<i,i,i>")
Q = Quat(-1, -1); T = T_diag(i, i, i)
for q, lab in [((0, 1, 0, 0), "b=1"), ((0, 1, 1, 0), "b=2"), ((0, 1, 2, 0), "b=5"), ((0, 1, 1, 1), "b=3"), ((0, 1, 1, 2), "b=6")]:
    r, _, det, Delta = analyse(Q, T, q)
    print(f"  q={q} {lab}: Nrd(q)={Q.nrd(q)} Delta={r['Delta']} det(phi) sqf={r['detphi_sqf']} split={r['split']} criterion={r['crit']}")
print("  embeds(Q(sqrt Delta) in D):", embeds(-1, -1, -1))
print("Example 2: D=(-2,-5) (disc 5), T=<i,j,k>")
Q = Quat(-2, -5); T = T_diag(i, j, k)
rnd = random.Random(3); nsplit = 0; ntest = 0; fields = set()
for _ in range(60):
    q = random_pure(rnd, 4)
    r, _, det, Delta = analyse(Q, T, q)
    ntest += 1; nsplit += r["split"]; fields.add(r["b"])
    assert r["split"] == r["crit"]
print(f"  Delta={r['Delta']}; embeds={embeds(Delta, -2, -5)}; {ntest} random K (fields -b, b in {sorted(fields)[:12]}...): split {nsplit}")
print("Example 3: D=(-3,-5) (disc 5, contains sqrt-3), T=<i,j,k>")
Q = Quat(-3, -5); T = T_diag(i, j, k)
r, _, det, Delta = analyse(Q, T, (0, 1, 0, 0))
print(f"  K=Q(sqrt-3): Delta={r['Delta']} split={r['split']} crit={r['crit']}; embeds={embeds(Delta, -3, -5)}")
print("Example 4: D=(-1,-7) (disc 7), T=<i,j,k>")
Q = Quat(-1, -7); T = T_diag(i, j, k)
for q, lab in [((0, 1, 0, 0), "Q(i)"), ((0, 0, 1, 0), "Q(sqrt-7)"), ((0, 0, 0, 1), "Q(sqrt-7) via k")]:
    r, _, det, Delta = analyse(Q, T, q)
    print(f"  {lab}: Delta={r['Delta']} split={r['split']} crit={r['crit']}")
print("Example 5: D=(-1,-3) (disc 3), T=<i,j,k>, K=Q(sqrt-3)")
Q = Quat(-1, -3); T = T_diag(i, j, k)
r, _, det, Delta = analyse(Q, T, (0, 0, 1, 0))
print(f"  Delta={r['Delta']} split={r['split']} crit={r['crit']}")
