enumeration/certificate/ -- certificate for Section 5 of the note, with an independent checker
=============================================================================================

PURPOSE
Section 5 of the note (Proposition 5.1, Theorem 5.2, Proposition 5.3 and Table 2) rests on the enumeration of
reduced pairs (G, Phi) with g <= 6.  certificate.json records, for every reduced pair with defect d >= 1 (one
record per G-orbit of CM types, 175 records), the data needed to verify the statements of Section 5 without
the enumeration programs, and check_certificate.py verifies every record from certificate.json alone.

FILES
  certificate.json         the certificate (175 records).
  make_certificate.py      writes certificate.json from ../enum_g{g}.json; uses ../verify_enum.py (pure
                           Python) only to find the Weil characters.
  check_certificate.py     the checker; pure Python 3, standard library only, reads only certificate.json.
  check_certificate_out.txt   its output.
  mutation_test.py         corrupts records in memory and counts how many the checker rejects.
  mutation_test_out.txt    its output.

RECORD FORMAT (certificate.json, key "records")
  g, gid, Nclass           dimension; index of the group in ../enum_g{g}.json; N_W(G)-class (as there)
  gens                     generators of G as permutations of {0, ..., 2g-1} (0-based lists, p[i] = image of
                           i); complex conjugation is rho(i) = i + g mod 2g
  orbits                   the G-orbits (the simple factors A_i)
  Phi                      the CM type (smallest representative of its G-orbit)
  d, LU                    the defect and a basis (Hermite normal form) of Lambda_U in Z^g; coordinate i
                           corresponds to the pair {i, i+g}, as in ../README.txt
  status "kinds a-d":      weil = list of Weil characters, each given by
                             chi  values +1/-1 of a homomorphism G -> {+1,-1} on the generators (it defines
                                  the imaginary quadratic field K; chi(rho) = -1)
                             c    {orbit index: c_i}, the multiplicities c_i of the factors A_i in B (the sign
                                  of c_i chooses the embedding of K into E_i, as in the note)
                             kind a (dim B = 2), b (dim B = 4), c (dim B = 6 with a factor of odd dimension)
                                  or d (B = A_i, a simple factor of dimension 6)
                           T = integer matrix with LU = T * (the Weil characters), row by row
  status "table2":         one of the 48 G-orbits of Table 2 (no Weil characters listed)

WHAT check_certificate.py VERIFIES (for every record)
  0. All entries are integers of the right shape: gens, orbits, Phi (strictly increasing, in 0..2g-1), d, LU, and
     for "kinds a-d" chi (+1/-1), the multiplicities c_i (non-zero) and T.
  1. G = <gens> commutes with rho and contains rho; the stored orbits are the G-orbits; Phi is a CM type.
  2. The pair is reduced: on each orbit Phi is primitive (no block system of G on the orbit, with rho
     permuting the blocks without fixed points, has Phi as a union of blocks), and no G-set isomorphism
     between two orbits maps Phi onto Phi.
  3. d = g - rank{mu_h : h in G}; the rows of LU are orthogonal to every mu_h and span a saturated lattice of
     rank d (gcd of the maximal minors is 1), so they form a basis of Lambda_U.
  4. "kinds a-d": chi is a homomorphism with chi(rho) = -1; each orbit used splits into two blocks for the
     kernel of chi; the sub-product is balanced (sum c_i t_i = 0); its dimension and kind are as stated;
     LU = T * (characters) with T an integer matrix.  Hence Lambda_W = Lambda_U for these Weil characters.
  5. "table2": the Weil characters of all imaginary quadratic fields on all balanced sub-products of powers of
     dimension <= 6 (all multiplicities), which are checked to lie in Lambda_U, do not generate Lambda_U.  For
     these records it also prints the rank of the lattice generated in dimension <= 12 and whether it equals
     Lambda_U (column "IQ" of Table 2).
  Over all records:
  6. No two records with the same group G have CM types in the same G-orbit.
  7. For the "table2" records (family by family) and for the 21 transitive records with g = 6: the classes up
     to conjugacy in W(B_g), computed by brute force over the 46080 elements of W(B_6) ((G_1, Phi_1) and
     (G_2, Phi_2) are conjugate if some w conjugates G_1 onto G_2 and maps Phi_1 into the G_2-orbit of Phi_2),
     coincide with the classes given by the stored labels (gid, Nclass).  The numbers printed for Table 2 and
     for the transitive pairs (15 classes, 21 G-orbits; Proposition 4.1(iii)) are these brute-force counts.
  NOT VERIFIABLE FROM THE FILE: that the records include every reduced pair with d >= 1 and g <= 6.  This rests
  on the enumeration (../README.txt) and its independent recomputation.

HOW TO RUN (from this folder)
  python make_certificate.py                      (writes certificate.json)
  python check_certificate.py > check_certificate_out.txt      (about 3 seconds)
  python mutation_test.py > mutation_test_out.txt              (about 6 seconds)

EXPECTED OUTPUT
  make_certificate.py: "records: 175 [((4, 'kinds a-d'), 6), ((5, 'kinds a-d'), 21), ((6, 'kinds a-d'),
    100), ((6, 'table2'), 48)]"
  check_certificate_out.txt:
    "records: 175 ; with errors: 0"
    Table 2: (5,1) 5 (10), rank 1, equal; (4,1,1) 3 (12), rank 2, not equal; (4,2) 8 (22), rank 0;
    (2,2,2) 1 (4), rank 0
    "transitive pairs with g = 6 (simple sixfolds, d = 1): 15 classes up to W(B_6) (21 G-orbits)"
    "RESULT: all records verified"
  mutation_test_out.txt: every corruption of LU, d, T (a changed entry, or an entry of type float), a kind
    label, a multiplicity (changed, or an extra factor with multiplicity 0), the status or the range of Phi
    (an element x replaced by x - 2g) is rejected; so are the 12 false certificates with a non-integral T
    built from the imaginary-quadratic Weil characters of the Table 2 records whose characters span
    Lambda_U (x) Q (the family (4,1,1)), each only because T is not integral; a second record for the same
    pair is rejected in all 175 cases, and a changed class label in all 14 cases where one can be changed
    (1475 of 1475 must-reject mutations).  A changed sign of chi is rejected in all 127 cases; replacing an
    element of Phi by its conjugate is rejected in 151 of 175 cases, and the 24 accepted ones give a CM type
    in the G-orbit of the original one (the same pair).
    "RESULT: all must-reject mutations detected"
  Steps 0, 6 and 7 and the corresponding mutations were added after an internal review of the code
  (7 October 2026) found that an earlier version accepted a non-integral T and out-of-range CM types.

WORKED EXAMPLE (Remark 5.4(a) of the note): the record g = 6, gid = 270, Phi = [0,1,2,3,4,11] (family
(4,1,1); orbits [0,1,4,5,6,7,10,11], [2,8], [3,9]).  LU = [[1,1,0,-2,-1,-1],[0,0,1,-1,-1,-1]]; the imaginary
quadratic Weil characters (dimension <= 6, and also <= 12) are +-(1,1,0,-2,-1,-1) and +-(-1,-1,2,0,-1,-1),
which generate a sublattice of index 2 in Lambda_U.

SOFTWARE
  Python 3.12.9 (standard library only).  Written and run on 7 October 2026.
