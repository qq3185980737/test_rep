# Released regression corpus (349 inputs)

This directory holds the released regression corpus used for the bitwise
verification reported in the manuscripts (`docs/paper_cpc_manuscript.md`,
`docs/softwarex_manuscript.md`, `docs/ijqc_manuscript.md`).

## Layout

- `inputs/` — 349 MOPAC input decks (`*.mop`), covering single points,
  geometry optimizations, transition-state searches, COSMO runs, CI, and
  molecular-dynamics/DRC trajectories across the PM7/PM6/AM1/RM1/MNDO
  Hamiltonians and the MOZYME linear-scaling method.
- `out_cpp/` — 349 recorded outputs from the released C++ build
  (`mopac.exe`, Windows x64, single-threaded: `OMP_NUM_THREADS=1;
  MKL_NUM_THREADS=1`).
- `out_ftn/` — 327 `.arc` summary files from the Fortran reference rebuild
  (Intel ifort/ifx, the verification build pair) for the same inputs
  (the remaining 22 inputs are MD/DRC trajectories that do not produce an
  `.arc` summary).

## Bitwise comparison result

`FINAL HEAT OF FORMATION` values were extracted from both sides and compared
string-exactly at printed precision:

- 327 inputs are paired (both sides produced a heat of formation);
- **320 are bitwise identical** (single points and non-EF tasks);
- 7 are EF geometry-optimization tasks (`_mx_opt_*`, `_val_ef`) whose
  endpoints differ at the 1-ULP level (e.g., −8.78954 vs −8.78878,
  −57.79850 vs −57.79870 kcal/mol) — the documented ULP-chaos sensitivity
  of EF trajectories described in §3.2 of the CPC manuscript (build
  sensitivity, not translation error);
- 22 MD/DRC trajectories have no Fortran `.arc` counterpart.

Hence the "320/320 released-regression paired non-EF inputs, zero exceptions"
statement in the manuscripts refers to the 320 paired non-EF inputs
reproduced bitwise (317 of the 349 released decks are 1SCF single points).
The comparison is reproduced by
`cmp_corpus.py` in this directory (run `python cmp_corpus.py`).

## Provenance

The corpus is a snapshot of the working verification corpus
(`D:\mopac\_full_cmp\cur` on the authors' machine, 2026-10-14) with the
recorded C++ and Fortran outputs. Inputs were not modified for this release.
