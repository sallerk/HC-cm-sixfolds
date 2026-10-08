enumeration/prop41/ -- computer check of the hand proof of Proposition 4.1 (simple CM sixfolds)
================================================================================================

PURPOSE
Section 4 of the note proves Proposition 4.1 (d <= 1 for a primitive CM type of degree 12, and the
description of the cases d = 1) by hand, in Lemmas 4.2-4.5.  check_prop41.py tests each step of that
proof by direct computation and compares the results with the stored enumeration.  It is pure Python 3
(standard library only) and shares no code with the enumeration programs: it has its own group closure,
orbits, block systems and exact rational ranks.  The proof in the note does not use this program.

FILES
  check_prop41.py        the check; reads ../enum_g6.json and ../../examples/combi_rows.json.
  check_prop41_out.txt   its output.

WHAT IS CHECKED (labels as printed)
  A. All 142 W(B_6)-classes of transitive subgroups G containing rho (../enum_g6.json) and all 64 CM types
     of each (9088 pairs); for each G, that the stored generators commute with rho, generate a transitive
     group and contain rho:
     - Lemma 4.2: for every primitive Phi the 12 linear forms x_i, -x_i are pairwise distinct on the span of
       the Galois conjugates of mu_Phi (so 12 <= 2^s, s the dimension of the span); if two of them agree,
       Phi is imprimitive.  Also the minimum of s and the maximum of d over primitive types.
     - primitive types with d = 1: number and G-orbits, compared with the stored data.
     - Lemma 4.3, sextic case: for every block system with 6 blocks of size 2 permuted freely by rho, every
       CM type with exactly one element in each block is imprimitive.
     - Lemma 4.3(i) (parity): no CM type has signature (3,3) with respect to two imaginary quadratic
       subfields.
  B. Lemma 4.4: G = S_3 x C_2 acting simply transitively on 12 points, all 64 CM types: the formula
     d = [m_1 = 0] + [m_2 = 0] + 2(2 - rank A), where m = sum_s eps_s s in Q[S_3], m_1, m_2 are its
     images under the trivial and the sign character and A = sum_s eps_s std(s) (2 x 2; m_3 in the note); and
     "rank A <= 1 implies Phi imprimitive".
  C. Lemma 4.5 and Table 1: the 16 transitive groups H_0 of degree 6, built from their definitions
     (generators written out in the script), with their orders and their parity (contained in A6 exactly for
     6T4, 6T7, 6T10, 6T12, 6T15) compared with the standard tables; their block systems; G = H_0 x C_2, the
     CM types Phi_P of signature (3,3) with respect to K for the 20 subsets P of size 3; the criterion "Phi_P
     is primitive iff P is neither a block of size 3 nor a transversal of a system of blocks of size 2"
     against a direct primitivity test; the number n of primitive Phi_P, their G-orbits, the second imaginary
     quadratic field K' when H_0 has a system of two blocks of size 3; comparison with
     ../../examples/combi_rows.json (column n_prim_balanced) and with the stored classes in ../enum_g6.json.
     The column "stored 12T" is copied from combi_rows.json (GAP labels); for 12T23 and 12T24 the script
     checks that the stored groups act on the 6 conjugate pairs as 6T7 and 6T8 (by order and parity).

HOW TO RUN (from this folder)
  python check_prop41.py > check_prop41_out.txt
  Runtime: about 9 seconds.

EXPECTED OUTPUT (key lines of check_prop41_out.txt)
  "transitive W(B_6)-classes of subgroups containing rho: 142"
  "primitive types: 8676 ; min s over primitive: 5 ; max d over primitive: 1"
  "functionals +-x_i coincide for a primitive type: 0 ; types where they coincide (all imprimitive): 412"
  "primitive types with d = 1: 332 in 21 G-orbits (stored: 21)"
  "types with one element in each fibre: 320 , of which primitive: 0"
  "parity: types balanced for two imaginary quadratic subfields: 0"
  "CM types: 64 ; d-formula holds for 64 ; rank A <= 1 for 40 (primitive among them: 0) ;
   primitive types: 24, all with d = 0 ; primitive and balanced for K or K': 0"
  Table: n = 12, 0, 12, 12, 18, 12, 12, 12, 18, 18, 12, 20, 18, 20, 20, 20 for 6T1, ..., 6T16.
  "groups with primitive types of signature (3,3): 15 ; G-orbits of primitive types with d = 1: 21"
  "RESULT: all checks passed"

SOFTWARE
  Python 3.12.9 (standard library only).  Written and run on 7 October 2026.
