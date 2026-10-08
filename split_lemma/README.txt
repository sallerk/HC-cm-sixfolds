split_lemma -- numerical evidence for the split-polarization lemma (note, Sec. 3, Theorem A)
==========================================================================================

Purpose
-------
Notation: the integer d with K = Q(sqrt(-d)) is called b in the note.
Theorem A: every CM abelian variety of Weil type (E = K (x) F a CM field or a product of CM fields,
K = Q(sqrt(-d)) embedded diagonally, signature (n,n)) has an E-compatible polarization whose K-hermitian
form h_a(x,y) = Tr_{E/K}(a x ybar), a in F^x, has determinant (-1)^n modulo norms from K (is split).
Since det h_a = N_{F/Q}(a) d_F, the condition is delta(a) := (-1)^n N(a) d_F in Nm(K^x), with the signs
of a at the 2n real places prescribed by the CM type (exactly n minus signs).
The scripts test this for n = 3 on 80 test algebras and all C(6,3) = 20 sign patterns of signature (3,3),
and check the forced-prime identity used in the proof. The computations support the theorem; they do
not prove it.

Files
-----
verify_split.py        Main computation (PARI via cypari); supports Theorem A (instances with n = 3)
                       and the forced-prime step of its proof. Builds the 80 test cases:
                         14 cases F = Q(zeta_m)^+ for m = 13, 21, 28, 36 (with several d each);
                         52 cases F = a random totally real sextic field (48 distinct fields; the
                            generator asks for up to 60 but stops after 400000 random polynomials,
                            having found 52);
                         14 product cases with factor degrees (3,3), (2,4), (2,2,2), (1,5), (1,1,4).
                       d takes 20 values between 1 and 71; max d_F = 2699728625.
                       For every case and sign pattern it searches a in F^x with these signs and
                       delta(a) in Nm(K^x): first up to 400 random elements with small coefficients
                       (wrong signs rejected), then, if none was split, an S-unit phase that solves
                       the sign and Hilbert-symbol conditions F_2-linearly, enlarging S by auxiliary
                       primes when necessary (the Chebotarev step of the proof). Each witness
                       carries a certificate (X, Y, Z), Z != 0, with X^2 + d Y^2 = delta Z^2 (PARI
                       qfsolve). It also computes, for every p | 2 d d_F, whether p is forced (-d a
                       square in F_v for all v | p) and t_p = ((-1)^n d_F, -d)_p, and checks t_p = 1
                       at forced primes (proof of Theorem A, forced-prime step). Fixed seed:
                       random.seed(20261004). Writes split_witnesses.json (next to the script);
                       prints verify_split.out.
recheck_split.py       Independent check (sympy + mpmath, no PARI) of the same; for every witness in
                       split_witnesses.json: signs at the real roots (60 digits), exactly n minus
                       signs and coverage of all 20 patterns per case; N(a) as an exact resultant;
                       delta modulo squares via polynomial discriminants; the norm equation
                       X^2 + d Y^2 = delta Z^2 exactly; det h_a = N(a) disc on the power basis (first
                       witness of each case); t_p recomputed with its own Hilbert-symbol code.
forced_nontrivial.py   Counts the forced primes in split_witnesses.json and those among them where
                       -d is not a square in Q_p (K_p a field inside every F_v: the non-trivial case
                       of the forced-prime step). The stored forced_nontrivial.out was first produced
                       by an inline command; this file (added when the folder was prepared) contains
                       the same code and reproduces it byte for byte.
witness_stats.py       Descriptive statistics of split_witnesses.json that verify_split.py does not print
                       (case composition, witnesses per search phase, auxiliary primes, forced
                       primes, split fraction among correctly-signed random a). Added when the
                       folder was prepared; its last line repeats the inline command that produced
                       the split-fraction numbers. Not needed for Theorem A.
split_witnesses.json   80 cases x 20 patterns: signs, witness (a in the integral basis and as a
                       polynomial, N(a), delta, (X,Y,Z), auxiliary primes), counters, search phase;
                       forced-prime table per case.
verify_split.out, recheck_split.out, forced_nontrivial.out, witness_stats.out   stored outputs.
Third-party files: none (all files were written for this project).

How to run
----------
Python 3.12.9, cypari 2.5.6 (PARI 2.15.4), sympy 1.13.1, mpmath 1.3.0 (versions used for the re-run).
From this folder (the shell lines are the ones used to make the stored outputs; "EXIT" is appended by
the shell):
  python verify_split.py > verify_split.out 2>&1; echo EXIT $? >> verify_split.out
  python recheck_split.py > recheck_split.out 2>&1; echo EXIT $? >> recheck_split.out
  python forced_nontrivial.py > forced_nontrivial.out
  python witness_stats.py > witness_stats.out
verify_split.py overwrites split_witnesses.json; the other three scripts read it. verify_split.py
asks PARI for a 2 GB stack (pari.allocatemem); the first line of verify_split.out is PARI's message.

Expected output (key lines)
---------------------------
verify_split.out:      80 lines "... witnesses 20/20 ... (identity ok: True) ...", last line
                       "cases 80  all witnesses found: True  forced identity holds everywhere: True  time ..."
                       (the "forced primes" column sums to 118).
recheck_split.out:     "cases 80 witnesses rechecked 1600 problems 0"
forced_nontrivial.out: "forced primes 118  non-trivial (K_p a field inside every F_v): 17"
witness_stats.out:     witnesses 1600 of 1600; by search phase random 555, sunit 1045; 14 witnesses
                       needing auxiliary primes (6 with one, 8 with two); forced primes 118, all with
                       t_p = 1; "37 cases; split fraction among correctly-signed random a
                       (random-phase patterns): min 0.066 median 0.293 max 0.505".
Remarks on the statistics: the per-case "split-fraction" column of verify_split.out counts all
patterns (an S-unit witness counts as one split try), so it differs from the witness_stats figure,
which uses only patterns solved in the random phase (37 cases with >= 50 such tries). Patterns sent to
the S-unit phase had no split element among their correctly-signed random tries (0 to 38 tries, median
6, out of 400 attempts); these patterns are not included in the 6.6%-50.5% range.

Runtime (one core)
------------------
verify_split.py about 25 s; recheck_split.py about 7 s; forced_nontrivial.py and witness_stats.py < 1 s.
See RERUN_LOG.txt for the re-run and the comparison with the stored outputs.

Review programs
---------------
review/ contains the programs written for the second internal review of Theorem A (own Hilbert symbols, checks of
Lemma 3.1 including products of fields, further certified witnesses); see review/README.txt.
