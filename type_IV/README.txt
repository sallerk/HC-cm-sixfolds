type_IV -- computations for type-IV abelian sixfolds with L = K F (note, Sec. 7, Theorem D)
==========================================================================================

Purpose
-------
Setting (note, Sec. 7; the note writes b for d and A^k for the powers A^n below):
K = Q(sqrt(-d)), F = Q(sqrt(m)) real quadratic, L = K F, A an abelian sixfold
with L acting, Hodge multiplicities (2,1) and (1,2) at the two embeddings of L over one embedding of K;
H a 3x3 L-hermitian matrix (signature (2,1) at one real place f1 of F, (1,2) at the other), and for
a in F the K-hermitian form h_a = Tr_{L/K}(a H) on V = L^3, computed in the K-basis {e_i, sqrt(m) e_i}.
The scripts check (in exact arithmetic, except for the ranks mod p in hodge_ring_gen.py and the
floating-point ranks in hodge_pullback_span.py):
  - Lemma 7.1: det h_a = N(a)^3 N(det_L H) (4m)^3;
  - that the explicit element a0 = -sqrt(m) det_L H (sqrt(m) > 0 at f1) is totally positive and h_a0
    is split (also a0' = -sqrt(m)/det_L H), and two further constructions of a split h_a;
  - the side remark that type-IV loci of this kind lie inside non-split (K,3,delta) Weil families;
  - dimensions and generators of the Hodge ring of A^n for the Hodge group of Proposition 7.2
    (these scripts take that group as input and do not verify the proposition; this part is not
    needed for the Hodge conjecture statement).
"Split" is tested in two ways: -det h_a is a norm from K (Hilbert symbols), and the 12-dimensional
rational quadratic form x -> h_a(x,x) is hyperbolic (Witt index 6, or signature, discriminant and Hasse
invariants in the independent check). The computations support the statements; they do not prove them.

Files
-----
biq.py                   Library used by split_checks.py, p0_check.py, degenerate_cases.py and
                         nonsplit_loci.py: exact arithmetic in L (Fractions), 6x6 Gram matrix of h_a over K
                         and its determinant (PARI matdet over K), Hilbert-symbol norm test, forced primes,
                         relative norm equations (rnfisnorm), Witt index of a rational quadratic form by
                         repeated qfsolve (find an isotropic vector, split off a hyperbolic plane, take the
                         complement with matkerint).
split_checks.py          Supports Lemma 7.1 and the existence of a split polarization (Theorem D).
                         Random L-hermitian H with the required signatures for all 156 pairs (d,m),
                         d in {1,2,3,5,6,7,11,15,19,23,31,35}, m in {2,3,5,6,7,10,11,13,14,15,17,21,30}
                         (d = m included), 5 per pair = 780 cases. For each H: computation A (local-global:
                         Hilbert symbols at forced primes, then qfsolve on u^2 + d v^2 - c x^2 + c m y^2,
                         c = -N(det_L H) d_F), computation B (explicit: an H-isotropic vector from a relative
                         norm equation in a small search box, a = sqrt(m)/delta and an explicit
                         3-dimensional h_a-isotropic K-subspace W), and the checker (exact Gram determinant
                         vs Lemma 7.1, Hilbert test, Witt index) on both elements and on the control a = 1.
                         Usage: python split_checks.py PER I K  (chunk I of K: pairs[I::K], seed 51005+I).
                         Writes split_checks_I.json; stdout is split_checks_I.out. In the "case" lines the
                         last number is elapsed seconds.
combine_split_checks.py  Sums the STATS/FAILS lines of the five chunk outputs into split_checks.out
                         ("seconds" there is the sum of the five chunk times). First done by an inline
                         command; this file (added when the folder was prepared) has the same code.
p0_check.py              Supports Lemma 7.1 and the explicit element a0. Fresh random H (3 per pair,
                         468 cases, seed 424242): a0 and a0' totally positive, -det h_a0 a rational square
                         (so det h_a0 = -1 modulo squares), checker on a0 and a0'.
degenerate_cases.py      Supports the local-global argument for a split h_a and its product case; not
                         needed for a0. (1) F = Q x Q (L = K x K): the construction of computation B on
                         each factor (isotropic vectors by qfsolve), 15 cases for each d in
                         {1,2,3,5,7,15,23} (seed 77); (2) census of forced primes (-dm a square in Q_p,
                         -d not) by local type over the (d,m) grid, 3 H per pair (generator of
                         split_checks.py, seed 51005), with the symbols (N det_L H, -d)_p and (-d_F, -d)_p;
                         (3) obstruction control: replacing c by c q (q inert in K, -dm a square mod q)
                         makes qfsolve report a local obstruction, while c itself is solvable.
nonsplit_loci.py         Supports the side remark on type-IV loci inside non-split Weil families. For
                         d in {1,2,3,5,7,11} and 16 values delta < 0 (96 pairs): F = Q(sqrt(M)),
                         M = y^2 + delta, det_L H = M + y sqrt(M), H = diag(1, -1, -det_L H); checks that
                         h_1 is in the class of delta (so the family lies in the (K,3,delta) Weil family),
                         that a0 is totally positive and that h_a0 is split.
independent_split_check.py  Supports Lemma 7.1 and the element a0 by an independent implementation (does
                         not import biq.py or split_checks.py). 242 random H (121 pairs (d,m), 2 rounds,
                         seed from the command line, 20261005 for the stored run); for a0: total
                         positivity, Lemma 7.1, Hilbert-norm test, hyperbolicity of the 12-dim trace form
                         via signature, discriminant and Hasse invariants; controls a = 1 and a random
                         totally positive a (the two split tests must agree). It first prints a CM example
                         (Kubota rank of a CM type of M = L*C3 with the same multiplicities: MT rank 6,
                         dim Hg 5), which concerns members with End^0 larger than L (almost-neatness); it
                         is not needed for Lemma 7.1 or a0. The script stops adding cases after 420 s
                         (BUDGET); a run takes about 2 s, so all 242 cases are done.
                         Usage: python independent_split_check.py SEED [JSON_PATH] (default JSON file:
                         independent_split_check.json next to the script); the summary goes to stdout.
                         Fix applied before the stored run: the hyperbolicity test first called
                         pari.issquare, which cypari does not provide (AttributeError); it was replaced
                         by the integer test math.isqrt(prod)**2 == prod.
hodge_dims.py            Hodge ring for the group of Proposition 7.2 (over C: {(g1,g2) in GL3 x GL3 :
                         det g1 det g2 = 1}): dim B^k(A^n) of Hodge classes and dim D^k(A^n) of Lefschetz
                         classes, n = 1, 2, 3, by the Weyl character formula. Not needed for the Hodge
                         conjecture statement.
hodge_ring_gen.py        Independent check of those numbers and of generation: ranks (mod p = 2^31-1) of
                         the ring generated by divisor classes and the pullbacks f^*W(A,K), in an exterior
                         algebra model. Usage: python hodge_ring_gen.py N KMAX.
hodge_pullback_span.py   Numerical rank of span{f^* Omega : f in L^n} for n = 1..4 (part1_pullback_span of
                         hodge_ring_gen.py with 20/60/250/600 random trials). Part 1 of
                         hodge_ring_gen_n3.out reports "60 of 100" only because its default of 60 trials
                         caps the rank; hodge_pullback_span.out has the full rank 100. First run as an
                         inline command; this file (added when the folder was prepared) has the same code.
*.out, *.json            Stored outputs (see below).
Third-party files: none (all files were written for this project).

How to run
----------
Python 3.12.9, cypari 2.5.6 (PARI 2.15.4), numpy 2.1.3 (versions used for the re-run). From this folder
(these are the shell lines that produced the stored outputs; "DONE_..." lines are appended by the shell):
  python -u p0_check.py 3 > p0_check.out 2>&1; echo DONE_P0 >> p0_check.out
  for i in 0 1 2 3 4; do python -u split_checks.py 5 $i 5 > split_checks_$i.out 2>&1; echo "DONE_CHUNK rc=$?" >> split_checks_$i.out; done
  python combine_split_checks.py
  python degenerate_cases.py > degenerate_cases.out 2>&1
  python nonsplit_loci.py > nonsplit_loci.out 2>&1
  python independent_split_check.py 20261005 independent_split_check.json > independent_split_check.out 2>&1
  python hodge_dims.py > hodge_dims.out 2>&1
  python hodge_ring_gen.py 2 12 > hodge_ring_gen_n2.out 2>&1
  python hodge_ring_gen.py 3 4 > hodge_ring_gen_n3.out 2>&1
  python hodge_pullback_span.py > hodge_pullback_span.out
The five chunks were first run in parallel; they are independent and give the same results when run
one after another. biq.py asks PARI for a 2 GB stack.

Expected output (key lines)
---------------------------
split_checks.out (sum of the five chunks): cases 780, A_ok 780, C_A_ok 780, B_ok 592, C_B_ok 592,
  B_notfound 188 (computation B found no isotropic vector in its search box; not a failure),
  det_formula_ok 780, h1_split 84 (a = 1 is split in 84 of 780 cases), h1_consistent 780 (the two split
  tests agree for a = 1), forced_cases 203, forced2_cases 95 (p = 2 forced), FAILS 0.
p0_check.out: {"cases": 468, "totpos": 468, "det_is_minus_square": 468, "checker_ok": 468,
  "checker_ok_inv": 468, ...}, bad 0.
degenerate_cases.out: product case 15/15 OK for each of the 7 values of d (105/105); forced primes by
  type: inert-both 15, ramified-both 51, p=2 ramified-both 21, p=2 inert-both 36 (123 in all), all
  symbols 1; obstruction control: qfsolve reports an obstruction for c q in 8/8 tests (at q in five of
  them, at 7, 5 and 3 in the other three) and finds a solution for c in 8/8.
nonsplit_loci.out: "cases 96 nonsplit classes among them 58 ALL OK".
independent_split_check.out: CM example n_tau 3, (n_s1, n_s2) = (2, 1), Kubota rank 6, dim Hg 5;
  cases 242, a0_totpos = a0_formula = a0_hilbert = a0_witt = 242, c1_split 33, cR_split 31,
  c1_agree = cR_agree = 242, fails [].
hodge_dims.out: A: dim B^k = 1,2,3,6,3,2,1 (2 classes beyond the Lefschetz ones, in k = 3);
  A^2: 1,8,36,152,376,656,844,656,376,152,36,8,1; A^3: 1,18,171,1340,7335,... (full list in the file).
hodge_ring_gen_n2.out: generated ring = Hodge classes in all degrees k = 1..12, "ALL OK".
hodge_ring_gen_n3.out: the same for k <= 4, "ALL OK".
hodge_pullback_span.out: ranks 1, 16, 100, 400 for n = 1..4 (= C(n+2,3)^2).

Runtime (one core)
------------------
split_checks.py about 12.5 min for the five chunks run one after another (one case, chunk 1 case 140,
d = 31, m = 11, takes about 8 min by itself); p0_check.py 13 s; each other script under 5 s.
Details and the comparison with the stored outputs: RERUN_LOG.txt.

Review programs
---------------
review/ contains the programs written for the second internal review of Theorem D (Lemma 7.1, a0, consequences of
Proposition 7.2 at explicit period points, Hodge rings for n <= 4, the example of Remark 7.3); see review/README.txt.
