# HC-cm-sixfolds

Code and data supporting the note *Split Weil structures and the Hodge conjecture for some abelian sixfolds and their
powers* (Kevin Saller, 2026).

**Results.**
- Every abelian variety of CM type with a Weil structure has a split polarization.
- Assuming Markman's theorem on abelian sixfolds of split Weil type (arXiv:2502.03415): the Hodge conjecture for all
  powers of simple CM abelian sixfolds, of certain sixfolds of type III and of biquadratic type IV, and of CM abelian
  varieties of reduced dimension at most 5; the Tate conjecture for all powers of some ordinary abelian varieties over
  finite fields.

**Contents.** Each folder has a README.txt with rerun instructions; the reference outputs are stored next to the code.
- `enumeration/`: CM types of abelian varieties of dimension g <= 6: enumeration, defect lattices and the
  Weil-character lattice test.
- `split_lemma/`: numerical evidence for the split-polarization lemma (Theorem A of the note).
- `type_III/`: abelian sixfolds of type III: discriminant criterion and Hodge-ring checks.
- `type_IV/`: type-IV abelian sixfolds with L = KF (Theorem D of the note).
- `examples/`: explicit Weil polynomials of ordinary, geometrically simple abelian sixfolds of angle rank 5.

**Status:** draft, not yet peer reviewed.

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
