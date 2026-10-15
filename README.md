# MOPAC-CPP — Bitwise-Reproducible C++17 Translation of MOPAC2016

This repository hosts **MOPAC-CPP**, a complete, 1-based, bitwise-verified C++17
translation of the MOPAC2016 semiempirical molecular-orbital package
(679 Fortran source files, ≈120,000 lines), released under the MIT license.

## Repository layout

- `mopac_cpp_public/` — the C++17 source tree, build configuration
  (`MOPAC_Cpp.sln`), documentation (`README.md`, `KEYWORDS.md`,
  `UPSTREAM_NOTICE.md`), and tests.
- `mopac_cpp_public/regression/` — released regression corpus: 349 input decks,
  349 recorded C++ outputs, 327 Fortran `.arc` summaries, and the comparison
  scripts.
- `mopac_cpp_public/docs/` — accompanying manuscripts and figures.

## Highlights

- **Bitwise fidelity on the nominal path**: identical to a local Intel-ifort
  rebuild of the official MOPAC2016 source for 18 small molecules in
  single-point mode (320/320 released-corpus paired non-EF inputs, zero
  exceptions), and across SCF, CI, DRC, and multi-Hamiltonian (PM7/PM6/RM1)
  trajectories.
- **Determinism gates (all default OFF)**: FBDR (deterministic parallel Fock
  construction), GAB (deterministic degenerate-basis choice in the EF
  optimizer), and an optional exact-diagonalization mode (`FULLDIAG`).

## Build

See `mopac_cpp_public/README.md` for build and runtime instructions
(MSVC 2022 x64, Intel MKL 2025.3, SSE2 baseline; single-threaded by default
for bitwise-comparable runs).

## License

MIT — see `LICENSE.txt`. This is an independent translation project, not
affiliated with the MOPAC authors. The Fortran reference (MOPAC2016) is
governed by its own terms (see `mopac_cpp_public/UPSTREAM_NOTICE.md`).

## Manuscripts

Three companion submissions describe this work: a software description
(SoftwareX), a methodology paper (Computer Physics Communications), and a
study of exact versus pseudo diagonalization (Int. J. Quantum Chem.).
