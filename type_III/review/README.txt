type_III/review/ -- programs written for the second internal review of Theorem C (Section 6 of the note)
=========================================================================================================

These programs were written separately from the main code in type_III/ and share no code with it. They support
Lemma 6.1, Proposition 6.2, Lemmas 6.3 and 6.4, Proposition 6.5 and Remark 6.6. All were run on 6 October 2026
(Python 3, python-flint, SymPy; cypari only where stated). Thread count: environment variable NPROC (<= 7 used).

Proposition 6.2 and Lemma 6.1 (splitness criterion)
  crit_check.py        -> crit_check_out.txt         3120 random triples (D,T,K), mostly non-diagonal T; phi built
                                                     from psi = Trd o T; own Hilbert symbols (cross-checked by the
                                                     product formula, SymPy and explicit global solutions).
                                                     Expected: FINAL stats with "mism_crit": 0, "lemma_fail": 0;
                                                     the controls (wrong formulas) disagree, as they should.
  crit_nonembed.py     -> crit_nonembed_out.txt      600 K for forms with Q(sqrt Delta) not in D: 0 split;
                                                     positive control 57/600 split.
  crit_examples.py     -> crit_examples_out.txt      hand examples (e.g. D = (-3,-5), D = (-1,-7), T = <i,j,k>).
  compare_part_i_with_project.py [N] -> compare_part_i_with_project_out.txt
                                                     compares with the stored outputs ../data.json and
                                                     ../path1_out.json (all 21168 triples, or N random ones).
                                                     Expected: FINAL "split_disagree": 0.
  realize_delta.py     -> realize_delta_out.txt/.json  every Delta = -m (m squarefree, m <= 60) occurs for each of
                                                     the 18 algebras: 666/666 (Remark 6.6 now has a proof).

Lemma 6.3 (Hodge group) and Lemma 6.4 (the character rho)
  so6_subgroups.py     -> so6_subgroups_out.txt      faithful irreducible 6-dimensional representations of semisimple
                                                     Lie algebras of rank <= 3; exactly one is orthogonal (A3).
  lemma_rho_check.py   -> lemma_rho_check_out.txt    an element of L(A)(Q) outside the identity component acts by -1
                                                     on W(A,K) (10/10); controls act trivially (20/20).

Proposition 6.5 (generation) -- checks only; the proof in the note uses no computation
  inv_dims.py K DMAX   -> inv_dims_k*_d*.json        dimensions of SO(6)-invariants per dominant weight of GL(U)
                                                     (Weyl character formula).
  gen_check.py K DMAX MODE [DMIN] -> gen_check_k{K}_d{DMIN}-{DMAX}_{MODE}.txt/.json
                                                     ranks modulo primes of random projections of products of the
                                                     classes theta(u,u') and omega(u), compared with inv_dims.
                                                     MODE normal: k=2 all degrees, k=3 degrees <= 16, k=4 degrees
                                                     <= 10, all weights OK (k=3, degree 16: 55 weights,
                                                     3111 = 3111, per-weight times sum to 3.2 h; for k=3 the
                                                     degrees >= 20 follow by hard Lefschetz, and the middle
                                                     degree 18 was not run).
                                                     MODE theta_only (control): divisor classes alone fall short
                                                     exactly in degrees 6..18 for k=2. MODE bad (control): a
                                                     planted non-invariant generator is detected.
  compare_with_project.txt                           the 58 dominant weights for k=2 also computed by
                                                     ../hodge_ring.sage agree.
  pullback_span.py, pullback_span_big.py, pullback_span_k5.py -> pullback_span*_out.txt/.json
                                                     rank modulo a prime near 5e7 of the span of the pull-backs f^*w,
                                                     f in D^k: 7, 84, 462, 1716 for k = 1..4 (five pairs (D,K)) and
                                                     5005 for k = 5 (two pairs), i.e. dim Sym^6(C^{2k}). Controls:
                                                     projections alone give 7k; a random 6-dimensional subspace gives
                                                     more than dim Sym^6 (96 > 84, 474 > 462). pullback_span_k5.py
                                                     takes about 15 minutes per pair and several GB of memory.

Usage examples:  python crit_check.py ;  python gen_check.py 2 12 normal ;  python pullback_span_k5.py
