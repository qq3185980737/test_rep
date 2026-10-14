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
- **Bitwise fidelity on the nominal path**: the C++ port is bitwise identical to
  the reference for eighteen small molecules in single-point mode — 320/320
  fixed-geometry single-point runs with zero exceptions, covering nine elements
  (H, C, N, O, F, Cl, S, P, Si) and both closed-shell RHF and open-shell UHF —
  and across SCF, configuration-interaction, dynamic-reaction-coordinate, and
  multi-Hamiltonian (PM7/PM6/RM1) trajectories.
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

### Optional accelerated eigensolver: `FULLDIAG`

By default the SCF uses the official pseudo-diagonalization path, which is
bitwise identical to the Fortran reference. To switch to the full LAPACK
(`dsyevd`) eigensolver — faster for large systems (≈2–3.6× in 1SCF) and more
robust to near-degenerate SCF convergence (e.g. C₄₀H₈₂ optimization diverges on
the pseudo path but stays stable on the full path) — either add the `FULLDIAG`
keyword to the input deck (e.g. `PM7 FULLDIAG 1SCF`) or set the environment
variable:

    set MOPAC_FULL_DIAG=1
    mopac.exe <input>.mop

Both paths converge to mathematically equivalent optimized geometries
(NIST-verified: hexadecane ΔHf = -90.18 kcal/mol vs exp. -89.60), while the
default path keeps bitwise fidelity to the Fortran reference. The `FULLDIAG`
keyword and the `MOPAC_FULL_DIAG` variable are equivalent and interchangeable.

## Reproducibility gates (all default OFF)

| Variable | Effect |
|---|---|
| `MOPAC_DSPAR`, `MOPAC_DSPAR_NB` | FBDR deterministic parallel Fock reduction (NB buckets; NB=1 ≡ serial) |
| `MOPAC_GAB` | GAB deterministic degenerate-basis choice in the EF optimizer |
| `MOPAC_ULPPERT` | ±1-ULP perturbation injection (first EF cycle only; research use) |

With all variables unset the binary is bitwise identical to the official path.

## Verification suite

- `tests/` — input decks for the paper's validation molecules (eighteen small
  molecules, the C₁₂H₂₆ ULP-scan family, C₄₀H₈₂/C₈₀H₁₆₂/C₁₂₀H₂₄₂ performance
  set, COSMO regression molecules, and the FBDR/GAB test drivers).
- `tests/unit/` — per-module unit tests (one `test_*.cpp` per translated module).
- `run_test/` — auxiliary analysis scripts (trajectory, energy, and distance
  parsers; ULP-scan and COSMO regressions) together with run inputs and
  reference outputs.
- Layer-level comparisons (Hessian, eigenvectors, projected gradient, step),
  the ULP replay of disputed `dgemm` calls, the FBDR thread-independence matrix
  (1/4/8/16 threads), and the GAB 24-molecule regression were executed during
  the verification campaign; their inputs and recorded outputs are provided
  under `run_test/`, and the procedures are described in the manuscript.

The GAMESS CCSD(T)/CBS benchmark data cited in the paper are summarized in the
manuscript (Fig. 3 and §3.5); the full input/output files are available from the
authors on request.

## License

MIT (see `LICENSE`). This is an independent translation project; it is not
endorsed by, and is not affiliated with, the MOPAC authors or the MOPAC project.
MOPAC2016 is © its original authors; the Fortran reference remains governed by its
own terms.

## Citation

(placeholder — to be completed on publication)

    @article{...,
      title  = {Bitwise-Reproducibility Verification and Deterministic Numerical
                Methods for Scientific-Software Migration, with a Full MOPAC2016
                (Fortran to C++) Case Study},
      journal = {Computer Physics Communications},
      year    = {2026}
    }
