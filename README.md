# HC-cm-sixfolds

Code and data supporting the note *Split Weil structures on abelian sixfolds and the Hodge conjecture* (Kevin Saller,
2026), DOI [10.5281/zenodo.23199503](https://doi.org/10.5281/zenodo.23199503).

**Results.** Theorems A-D of the note do not use Markman's theorem.
- Theorem A: every abelian variety of CM type with a Weil structure has a split polarization.
- Theorem B: the Hodge group of a simple CM abelian sixfold has codimension at most one in the norm-one torus of its
  CM field, with equality in exactly 15 cases up to conjugacy.
- Theorem C: a criterion for split Weil structures on simple abelian sixfolds of type III.
- Theorem D: an explicit split polarization for abelian sixfolds with multiplication by a biquadratic CM field L with
  multiplicities (2,1) and (1,2); if End^0(A) = L, their Hodge group has dimension 17 and is not semisimple.
- Assuming Markman's theorem on abelian sixfolds of split Weil type (arXiv:2502.03415), Corollary E and further
  results: the Hodge conjecture for all powers of simple CM abelian sixfolds, of certain sixfolds of type III and of
  biquadratic type IV, and of CM abelian varieties of reduced dimension at most 5; the Tate conjecture for all powers
  of some ordinary abelian varieties over finite fields.
- Manuscripts released by OpenAI ([github.com/openai/math](https://github.com/openai/math)), which have not been
  refereed, claim the Hodge conjecture for all CM abelian varieties and for all powers of abelian sixfolds of split
  Weil type. If they are correct, the results above on the Hodge and Tate conjectures follow from theirs, in some
  cases together with the split polarizations above; the note's proofs do not use them (Secs. 1.2 and 1.3 of the
  note).

**Contents.** Each folder has a README.txt with rerun instructions; the reference outputs are stored next to the code.
- `enumeration/`: CM types of abelian varieties of dimension g <= 6: enumeration, defect lattices and the
  Weil-character lattice test; a computer check of the hand proof of Proposition 4.1 (`prop41/`) and a certificate
  for Section 5 with an independent checker (`certificate/`).
- `split_lemma/`: numerical evidence for the split-polarization lemma (Theorem A of the note).
- `type_III/`: abelian sixfolds of type III: discriminant criterion and Hodge-ring checks.
- `type_IV/`: type-IV abelian sixfolds with L = KF (Theorem D of the note).
- `examples/`: explicit Weil polynomials of ordinary, geometrically simple abelian sixfolds of angle rank 5.

**Status:** not peer reviewed.

**Disclosure:** the arguments and computations were produced with the assistance of an AI system (Claude, by
Anthropic). Most computations were re-checked with separately written code, and most arguments were checked in
separate reviews; the code and the reviews were also produced with AI assistance. What was not re-checked or not
reviewed is listed in the note and in the README files.

## Licences

- **Code** (source files such as `.py`, `.sage`, `.gp`, `.sh`): MIT License, see [`LICENSE`](LICENSE).
- **Everything else** (README files, data and computed outputs): Creative Commons Attribution 4.0 International
  (CC BY 4.0), see [`LICENSE-CC-BY-4.0`](LICENSE-CC-BY-4.0).
- The data of Arango-Piñeros, Frengley and Vemulapalli used in `enumeration/` and `examples/` are not included;
  `fetch_apfv_data.py` in each of these folders downloads them.
