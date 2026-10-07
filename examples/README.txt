examples/ -- explicit Weil polynomials: ordinary, geometrically simple abelian sixfolds of angle rank 5

PURPOSE
For each of the 15 degenerate Galois classes of [APFV] for ordinary abelian sixfolds (Galois group 6Tj x C2 in the
product action, j = 1, 3, 4, ..., 16) this folder gives an explicit Weil q-polynomial P(T) of degree 12 whose isogeny
class of abelian sixfolds over F_q is ordinary, geometrically simple, of angle rank 5 and has the labelled Galois
group, together with the scripts that construct the examples and verify them independently.
Supports:
  - note Sec. 8, Corollary 8.1 and Table 3 (the 15 rows *_deg of final_table.txt);
  - note Sec. 4, Example: row 12T7, F = Q[x]/(x^6 - 12x^4 + 35x^2 - 25), K = Q(i), E = F.K, q = p = 73.
[APFV] S. Arango-Pineros, S. Frengley, S. Vemulapalli, Galois groups of simple abelian varieties over finite fields
and exceptional Tate classes, arXiv:2505.09589.

RESULT (final_table.txt / final_table.json, rows *_deg)
  row     Gal(F)  F (totally real sextic)                    disc(F)          K          p      q
  12T2    6T1     x^6-x^5-7x^4+2x^3+7x^2-2x-1                300125           Q(sqrt-3)  139    139
  12T7    6T4     x^6-12x^4+35x^2-25                         45697600         Q(i)       73     73
  12T10   6T3     x^6-x^5-10x^4+x^3+12x^2-3x-1               2738000          Q(i)       269    269
  12T18   6T5     x^6-x^5-6x^4+7x^3+4x^2-5x+1                722000           Q(sqrt-3)  331    331
  12T23   6T7     x^6-6x^4+8x^2-1                            3356224          Q(i)       617    617
  12T24   6T8     x^6-11x^4+37x^2-37                         51868672         Q(i)       881    881
  12T25   6T6     x^6-8x^4+19x^2-13                          1997632          Q(sqrt-3)  211    211
  12T37   6T9     x^6-x^5-13x^4+10x^3+34x^2-39x+9            46360125         Q(sqrt-3)  331    331
  12T40   6T10    x^6-21x^4-21x^3+99x^2+198x+99              55130625         Q(sqrt-3)  991    991
  12T48   6T11    x^6-3x^5-2x^4+9x^3-x^2-4x+1                1387029          Q(sqrt-3)  907    907
  12T75   6T12    x^6-x^5-91x^4-123x^3+840x^2+1310x+446      37949591784976   Q(sqrt-3)  8317   8317^2
  12T77   6T13    x^6-2x^5-4x^4+8x^3+2x^2-5x+1               485125           Q(sqrt-3)  619    619
  12T123  6T14    x^6-3x^5-35x^4+30x^3+90x^2-74x-13          12200571365625   Q(sqrt-3)  12409  12409
  12T180  6T15    x^6-3x^5-9x^4+14x^3+27x^2-3x-10            170067681        Q(sqrt-3)  2437   2437
  12T219  6T16    x^6-2x^5-4x^4+6x^3+4x^2-3x-1               1134389          Q(i)       2069   2069
  The CM type S, E = Q[x]/(polE), pi (as an element of E), the Weil polynomial (coefficients of T^12, ..., T^0)
  and all check results are in final_table.txt / final_table.json (pi is in examples_in.json).
  Example, row 12T7: P(T) = T^12 + 36T^11 + 626T^10 + 5268T^9 + 79T^8 - 579096T^7 - 7290148T^6 - 42274008T^5
  + 420991T^4 + 2049341556T^3 + 17777298866T^2 + 74630577348T + 151334226289 (= 73^6); with pi_1, ..., pi_6 the
  roots in the upper half plane ordered by argument, pi_1 pi_2 pi_3 conj(pi_4) conj(pi_5) conj(pi_6) = -73^3.
  The 12T2 row has abelian E (E is a subfield of Q(zeta_105)); the other 14 rows have non-abelian Galois closure.
  Other rows: *_ctrl = control, an unbalanced CM type on the same E (ordinary, geometrically simple, angle rank 6;
  none was found for 12T123 with q <= 10^30); fermat21_jacobi_p43_1_1_19 and fermat21_jacobi_p127_1_1_19 = Jacobi
  sums of order 21 (jacobi21.py; 12T2, angle rank 5).

CONSTRUCTION (main computation: produce.py + produce_lib.gp, PARI/GP through cypari)
  E = F.K, F totally real sextic with Galois group 6Tj, K = Q(i) or Q(sqrt-3). Take p split completely in E; the 12
  primes of E above p are labelled (j, e): j = root of F mod p, e = root of K mod p (residues sorted). A CM type is
  S = (e_1, ..., e_6), e_j in {1, 2}; it is balanced (signature (3,3) with respect to K) iff exactly three e_j = 1.
  a = product of the primes (j, e_j); h = order of [a] in Cl(E); alpha = generator of a^h; u = alpha*conj(alpha)/p^h;
  pi = alpha/v if some unit v has v*conj(v) = u (q = p^h), else pi = alpha^2/u (q = p^(2h)); P = charpoly(pi).
  Checks in produce.py (exact unless stated):
   - the generator of a^h is checked against the ideal; pi*conj(pi) = q;
   - P integral and q-symmetric; all |roots|^2 = q numerically (PARI polroots, tolerance 1e-20);
   - P irreducible; ordinary: p does not divide the coefficient of T^6;
   - geometric simplicity: charpoly(pi^N) squarefree modulo nextprime(10^40) or nextprime(10^41),
     N = lcm{n : phi(n) divides |G|}, |G| = 2|Gal(F~/Q)|;
   - relation N_{E/K}(pi)^12 = q^36, in K;
   - angle rank = rank(M) - 1, M the integer matrix of p-adic valuations v_P(g(pi))/(h or 2h) for
     g in Gal(F~/Q) x Gal(K/Q) and the 12 primes P of E above p. The action of Gal(F~/Q) on the p-adic roots of F is
     computed with nfsplitting and polrootspadic (|Gal(F~/Q)| <= 72), is all of Sym(6) / Alt(6) for 6T16 / 6T15,
     and for 6T12 / 6T14 rank(M) is computed for all 720 relabellings and is the same for all of them.
  bnfinit (class group under GRH) is used only to find h and alpha; every output above is checked unconditionally.

INDEPENDENT VERIFICATION (verify.py: python-flint + mpmath; verify_gap.py: GAP; verify_K.py: FLINT + mpmath;
no PARI and no code shared with produce.py). Input: only P, p, q and |G|.
  V1 q-symmetry: P monic, a_{12-i} = q^(6-i) a_i, q a power of the prime p.
  V2 Weil property (exact): P(T) = T^6 h(T + q/T); Sturm sequences over Q show that h has 6 real roots in [-B, B]
     (B = Cauchy bound) and that h(y)h(-y), as a polynomial in z = y^2, has no root in (4q, B^2]; so all 12 roots
     of P have absolute value sqrt(q).
  V3 irreducibility: FLINT factorisation of P over Z.
  V4 ordinarity: p does not divide the coefficient of T^6.
  V5 geometric simplicity: Q_N = charpoly of x^N on F_l[x]/(P) is squarefree for l = 2^61 - 1 (l = 2^62 - 57 is
     tried if that fails), N = lcm{n : phi(n) divides |G|} (every root of unity in the splitting field has such
     an order), so no quotient of two distinct roots is a root of unity. N ranges from 32760 (12T2) to
     12376209948955935621739200 (12T219).
  V6 Galois group: GAP GaloisType(h) = 6Tj (verify_gap.py), and an exact certificate (verify_K.py): c(y) in Q[y]
     with c(y)^2 + d(y^2 - 4q) = 0 mod h(y), so Q(sqrt(-d)) is contained in E = Q(pi) and E = F(sqrt(-d)), with d as
     in the table (d = 1 for Q(i), 3 for Q(sqrt-3)). Since the Galois closure of F is totally real, the Galois group
     is 6Tj x C2 in the product action; its 12T label comes from combi_rows.out (GAP TransitiveIdentification).
     GaloisType(P) itself (degree 12) was run only for the rows 12T2 and 12T10 and agrees.
  V7 angle rank: LLL (FLINT) on (theta_1/2pi, ..., theta_6/2pi, 1), roots to 420 digits, finds exactly one integer
     relation for each *_deg row (none for the *_ctrl rows). Each relation n is certified: with
     A' = prod_{n_k>0} pi_k^{n_k} prod_{n_k<0} conj(pi_k)^{|n_k|} and m = sum |n_k|, x = A'^2 - q^m is an
     algebraic integer of degree <= |G| whose conjugates have absolute value <= 2q^m, so |x| < (2q^m)^-(|G|-1)
     forces x = 0; |x| was evaluated with 204 to 57,816 digits (mpmath floating point, not interval arithmetic).
     Lower bound (numerical): by the Gram-Schmidt norms of the reduced basis, a second independent relation would
     need Euclidean norm of about 2.3e49 to 7.6e49 or more, depending on the row. Exact lower bound: the
     valuation-matrix rank of the main computation.

FILES
  final_table.txt, final_table.json  the table (31 entries: 15 *_deg, 14 *_ctrl, 2 fermat21) combining the main
      computation (key main_computation / label main_ar), the verifier, GAP and verify_K; VERIFIED=True for all.
      Supports Table 3 and Corollary 8.1 (Sec. 8), and the 12T7 Example (Sec. 4).
  final_table.py      assembles final_table.* from examples_in.json, produce_runs/, verify_A/B.json,
      verify_gap.json, verify_K.json.
  produce.py, produce_lib.gp  main computation (construction and exact checks above). Writes
      produce_runs/produce_<tag>.json / .out.
  produce_runs/       the 15 runs of produce.py from which the table rows were taken (one per row, all primes and
      CM types of that run). Earlier search runs that did not supply a table row are not included.
  collect.py          selects, per row, the example of smallest q passing all checks of the main computation (and
      one control) from produce_runs/*.json and jacobi21.json -> examples_in.json, examples_table.txt.
  examples_in.json    verifier input (31 entries: P, p, q, |G|, F, K, E, S, pi).
  examples_in_A.json, examples_in_B.json  the two parts of examples_in.json used for verify_A.json / verify_B.json.
  examples_in_gap.json  examples_in.json with "try12": true for the 12T2 / 12T10 entries (input of verify_gap.py;
      reconstructed, see RERUN_LOG.txt).
  examples_table.txt  summary printed by collect.py.
  verify.py           independent verifier V1-V5, V7 -> verify_A.json, verify_B.json.
  verify_gap.py       GAP GaloisType (Sage docker) -> verify_gap.json (stdout: verify_gap.out).
  verify_K.py         exact sqrt(-d) certificate and 12T label (V6) -> verify_K.json.
  combi_rows.py       GAP (Sage docker): for each 6Tj, the group 6Tj x C2 on 12 points, its 12T label, all 64 CM
      types (primitivity, rank of the valuation-type matrix = angle rank + 1) and the GAP conjugacy test of the
      degenerate primitive types against the APFV row -> combi_rows.json, combi_rows.out. Shows that for j != 2
      these are the 15 APFV rows and that the balanced primitive types have angle rank 5 and lie in the APFV row
      (for j = 1, 3, 5, 9, 10, 13 there are further degenerate primitive types, unbalanced with respect to K; they
      lie in the same APFV row).
      Reads the APFV data file apfv_data/6-0.000_..._1.00.m (not included; see fetch_apfv_data.py).
  fetch_apfv_data.py  downloads that APFV file from the public APFV repository and checks its SHA-256.
  jacobi21.py, jacobi21.json, jacobi21.out  control: Jacobi sums J(chi^a, chi^b) of order 21 over F_43 and F_127;
      (1,1,19) gives ordinary, geometrically simple sixfolds of angle rank 5 (12T2); PARI valuation rank checked
      against the Stickelberger rank.
  main_check_12T7.py, main_check_12T7.out  independent check of the 12T7 polynomial with sympy and mpmath only
      (q-symmetry, irreducibility, |roots|^2 = q, ordinarity, PSLQ relation theta_1 + theta_2 + theta_3 - theta_4
      - theta_5 - theta_6 = -pi, no second relation with coefficients up to 1e12). Supports the Sec. 4 Example.
  search_F.py, search_F.json, search_F.out  search for totally real sextic F with each Galois group 6Tj (smallest
      discriminant found per j); source of F for all rows except 12T75 and 12T180.
  search_F2.py, extra_F2.json  search for totally real 6T12 / 6T14 sextics (resolvent sextics of quintics); source
      of F for 12T75. F for 12T180 is the LMFDB number field 6.6.170067681.1 (polynomial recorded in the run command;
      disc 170067681 and Galois group 6T15 re-checked by produce.py / verify_gap.py).
  RERUN_LOG.txt       commands, run times and comparison results of the re-run.
  .gitignore          excludes apfv_data/ (third-party data, downloaded by fetch_apfv_data.py).
  Third-party files: none included.

HOW TO RUN (from this folder; versions used for the re-run in RERUN_LOG.txt)
  Python 3.12.9 with cypari 2.5.6 (PARI 2.15.4), python-flint 0.9.0, mpmath 1.3.0 (+ gmpy2), sympy 1.13.1;
  Docker image sagemath/sagemath:latest (SageMath 10.10, GAP 4.15.1). Docker command, below written DOCKER:
    docker run --rm --cpus=1 -v "<path>/examples:/work" -w /work sagemath/sagemath:latest sage -python
  (Git Bash on Windows: prefix MSYS_NO_PATHCONV=1.) Each command overwrites the stored file(s) it reproduces;
  to compare with the stored files, run in a copy of this folder.
   python fetch_apfv_data.py
   DOCKER combi_rows.py > combi_rows.out                                                   (25 s)
   python jacobi21.py > jacobi21.out                                                        (< 1 s)
   python produce.py 6T1_K3 "x^6 - x^5 - 7*x^4 + 2*x^3 + 7*x^2 - 2*x - 1" "x^2+x+1" 1 3 20000
   python produce.py 6T5_K3 "x^6 - x^5 - 6*x^4 + 7*x^3 + 4*x^2 - 5*x + 1" "x^2+x+1" 5 3 20000
   python produce.py 6T6_K3 "x^6 - 8*x^4 + 19*x^2 - 13" "x^2+x+1" 6 3 20000
   python produce.py 6T9_K3 "x^6 - x^5 - 13*x^4 + 10*x^3 + 34*x^2 - 39*x + 9" "x^2+x+1" 9 3 20000
   python produce.py 6T10_K3 "x^6 - 21*x^4 - 21*x^3 + 99*x^2 + 198*x + 99" "x^2+x+1" 10 3 20000
   python produce.py 6T11_K3 "x^6 - 3*x^5 - 2*x^4 + 9*x^3 - x^2 - 4*x + 1" "x^2+x+1" 11 3 20000
   python produce.py 6T13_K3 "x^6 - 2*x^5 - 4*x^4 + 8*x^3 + 2*x^2 - 5*x + 1" "x^2+x+1" 13 3 20000
   python produce.py 6T3_K4_b "x^6 - x^5 - 10*x^4 + x^3 + 12*x^2 - 3*x - 1" "x^2+1" 3 12 200000 Q30
   python produce.py 6T4_K4_b "x^6 - 12*x^4 + 35*x^2 - 25" "x^2+1" 4 12 200000 Q30
   python produce.py 6T7_K4_b "x^6 - 6*x^4 + 8*x^2 - 1" "x^2+1" 7 12 200000 Q30
   python produce.py 6T8_K4_b "x^6 - 11*x^4 + 37*x^2 - 37" "x^2+1" 8 12 200000 Q30
   python produce.py 6T16_K4_b "x^6 - 2*x^5 - 4*x^4 + 6*x^3 + 4*x^2 - 3*x - 1" "x^2+1" 16 12 3000000 Q30
   python produce.py 6T12_K3_a "x^6 - x^5 - 91*x^4 - 123*x^3 + 840*x^2 + 1310*x + 446" "x^2+x+1" 12 12 3000000 Q30
   python produce.py 6T15_K3_lm "x^6 - 3*x^5 - 9*x^4 + 14*x^3 + 27*x^2 - 3*x - 10" "x^2+x+1" 15 12 3000000 Q30
   python produce.py 6T14_K3_b "x^6 - 3*x^5 - 35*x^4 + 30*x^3 + 90*x^2 - 74*x - 13" "x^2+x+1" 14 12 3000000 Q30
      (Q30 = 1000000000000000000000000000000 = 10^30; arguments: tag, F, K, j, number of split primes, largest p,
       largest q; 0.4 s to 9.6 s each; 6T12 / 6T14 read combi_rows.json; PARI stack 4 GB)
   python collect.py                                                                         (< 1 s)
   python verify.py examples_in_A.json verify_A.json                                         (136 s)
   python verify.py examples_in_B.json verify_B.json                                         (25 s)
   DOCKER verify_gap.py examples_in_gap.json verify_gap.json > verify_gap.out                (22 s)
   python verify_K.py examples_in.json verify_gap.json verify_K.json                         (14 s)
   python final_table.py                                                                     (< 1 s)
   python main_check_12T7.py > main_check_12T7.out                                           (1 s)
   python search_F.py > search_F.out                                                         (129 s)
   python search_F2.py                                                                       (4 s)
  Expected output: the stored files are reproduced (JSON by content); only timing values differ
  ("time"/"done" in produce_runs/*.out, "sec" in verify_gap.*, elapsed times in search_F.out).
  final_table.txt: 31 lines "... VERIFIED=True"; *_deg rows main_ar=5 verifier_ar=5 cert=[True].

NOTES
  - Exactness: the Weil polynomials, q-symmetry, Weil property (Sturm), irreducibility, ordinarity, geometric
    simplicity (mod l), the sqrt(-d) certificate and the angle-rank upper bound (certified relation) are exact; GAP
    GaloisType and PARI polgalois give the Galois group of F; the angle-rank lower bound is exact in the main
    computation (valuation-matrix rank) and numerical in the verifier.
  - produce.py uses bnfinit (GRH) only to find generators; GRH could only affect which examples are found
    (e.g. the minimality of q), not the checks.
  - Changes against the original scripts: paths (relative to the script directory; produce.py writes to
    produce_runs/), comments, and in final_table.py the names of the main-computation fields (now
    "main_computation" in the JSON, "main_ar" in the text). No algorithm, seed, parameter or output format was
    changed otherwise.

lmfdb/ (statement in Sec. 8 that the rows 12T23, 12T48, 12T77 and 12T219 already occur in the LMFDB over F_2)
  lmfdb_rows.json            result of a search of the LMFDB isogeny classes of abelian sixfolds over F_2 by the
                             Galois group of the Frobenius field, one entry per row of Table 1, retrieved on
                             5 October 2026: matches for 12T23 (6), 12T48 (12), 12T77 (40) and 12T219 (at least
                             1000), none for the other rows; "_q" records searches over F_43, F_73, F_127, F_139
                             (no matches). The search script is not included; the LMFDB limits automated access.
  lmfdb_examples_in*.json    the Weil polynomials of one LMFDB class for each of the four rows;
  verify_lmfdb*.json, verify_K_lmfdb*.json, verify_gap_lmfdb*.json
                             outputs of verify.py, verify_K.py and verify_gap.py on these four classes (ordinary,
                             geometrically simple, Galois group of the row).
