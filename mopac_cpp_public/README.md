# mopac_cpp — Bitwise-Reproducible C++17 Translation of MOPAC2016

`mopac_cpp` is a complete, 1-based, bitwise-verified translation of the MOPAC2016
semiempirical molecular-orbital package (679 Fortran source files, ~120,000 lines)
into C++17 (MSVC 2022, SSE2, Intel MKL). It accompanies the manuscript

> *Bitwise-Reproducibility Verification and Deterministic Numerical Methods for
> Scientific-Software Migration, with a Full MOPAC2016 (Fortran to C++) Case Study*
> (submitted to Computer Physics Communications)

and provides the verification pipeline, the ULP-attribution protocol, and the two
deterministic methods (FBDR, GAB) described therein.

## What this is

- A **faithful port, not a rewrite**: variable names, array dimensions, initial
  values, and 1-based indexing semantics are preserved from the Fortran source, so
  the C++ is auditable line-by-line against the original.
- **Bitwise fidelity on the nominal path**: for ten small molecules
  (NH₃, HF, HCl, CO, N₂, CO₂, CH₂O, CH₃OH, H₂O, CH₄) the FINAL HEAT of formation is
  bitwise identical to the official MOPAC2016 executable in both single-point (1SCF)
  and geometry-optimization modes, and across CI, DRC, and multi-Hamiltonian
  (PM7/PM6/RM1) trajectories.
- **The reference is the official MOPAC2016 binary** (Intel ifort + MKL, Windows x64),
  built from the original 2016 Fortran sources. The C++ translation is validated
  against it; it is *not* validated against later MOPAC releases (v22/v23), whose
  outputs differ.

## Build

- Visual Studio 2022 (MSVC v143, x64)
- Intel MKL 2025.3 (static) — used for diagonalization (`dsyevd`), density assembly
  (`dgemm`), and related BLAS/LAPACK calls
- Build configuration: `x64/Release` (SSE2, LTCG)

Open `mopac_cpp.sln` and build the `mopac` project, or use MSBuild:

    msbuild mopac_cpp.sln /p:Configuration=Release /p:Platform=x64

## Run

    OMP_NUM_THREADS=1
    MKL_NUM_THREADS=1
    mopac.exe <input>.mop          # output written to <input>.out

Set `OMP_NUM_THREADS=1; MKL_NUM_THREADS=1` for bitwise-comparable single-threaded
runs. The default build follows the bitwise-official numerical path.

## Reproducibility gates (all default OFF)

| Variable | Effect |
|---|---|
| `MOPAC_DSPAR`, `MOPAC_DSPAR_NB` | FBDR deterministic parallel Fock reduction (NB buckets; NB=1 ≡ serial) |
| `MOPAC_GAB` | GAB deterministic degenerate-basis choice in the EF optimizer |
| `EFPROBE`, `HESSPROBE` | non-intrusive diagnostic dumps (no numerical-path change when unset) |
| `MOPAC_ULPPERT` | ±1-ULP perturbation injection (first EF cycle only; research use) |

With all variables unset the binary is bitwise identical to the official path.

## Verification

- `tests/` — representative input decks (`.mop`) used in the paper (10 small
  molecules, CI/DRC cases, and the C₁₂H₂₆ optimization).
- The full comparison scripts (per-layer bitwise comparison, ULP replay,
  ULP-sensitivity scan, FBDR thread-independence matrix, GAB batch) are released
  as supplementary material with the manuscript.

## License

MIT (see `LICENSE`). This is an independent translation project; it is not
endorsed by, and is not affiliated with, the MOPAC authors or the MOPAC project.
MOPAC2016 is © its original authors; the Fortran reference remains governed by its
own terms. The GAMESS ab initio benchmark data referenced in the paper are
reproducible from the released input/output files accompanying the manuscript.

## Citation

(placeholder — to be completed on publication)

    @article{...,
      title  = {Bitwise-Reproducibility Verification and Deterministic Numerical
                Methods for Scientific-Software Migration, with a Full MOPAC2016
                (Fortran to C++) Case Study},
      journal = {Computer Physics Communications},
      year    = {2026}
    }
