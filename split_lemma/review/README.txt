split_lemma/review/ -- programs written for the second internal review of Theorem A (Section 3 of the note)
===========================================================================================================

These programs were written separately from the main code in split_lemma/ and share no code with it. No proof in
the note uses them. Run on 6 October 2026 (Python 3; cypari only where stated, and every PARI result re-verified).

  qsym.py                 own Hilbert symbols over Q_p and R, square classes (pure Python).
  test_qsym.py -> test_qsym.out
                          24853 comparisons with PARI's hilbert, 4000 product-formula checks, 9000 comparisons with
                          a brute-force search for primitive zeros mod p^k. Expected: "disagreements: 0".
  alg.py                  own exact arithmetic in F = prod Q[y]/(f_i) and E = K (x) F.
  test_psia.py -> test_psia.out
                          Lemma 3.1 on 120 random instances (fields of degree 2, 4, 6 and products of fields):
                          phi = 2d * conj(h_a), det h_a = N(a) disc, the Riemann-form sign criterion and the
                          signature. Expected: "failures: 0".
  inst.py, witness.py, herm.py
                          random instances; search for split polarizations psi_a (PARI's bnfinit and qfsolve used
                          only as search tools); hermitian forms over K with an explicit Witt decomposition.
  smoke_tests_console.txt console output of interactive runs of witness.py and herm.py: certified witnesses for
                          products of shapes (3,3), (2,2,2), (1,1) and for a sextic field; non-split samples
                          certified. A planned larger batch was not run.
