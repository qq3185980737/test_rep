# MOPAC-CPP: a bitwise-reproducible C++17 translation of the MOPAC2016 semiempirical molecular-orbital package

**SoftwareX — Original Software Publication**

## Highlights

- Complete C++17 translation of MOPAC2016 with bitwise-identical default path
- Preserves 1-based semantics, names, and initial values for auditability
- FBDR and GAB gates keep bitwise fidelity under parallelism
- Optional FULLDIAG exact-diagonalization mode for large systems

**Authors:** Yang Liu^a, Xun Wang^a,*

^a Shenyang Jianzhu University, Shenyang, China
*Corresponding author: wangxun71@aliyun.com (X. Wang); liuyang123@stu.sjzu.edu.cn (Y. Liu)

**Keywords:** semiempirical quantum chemistry; PM7; Fortran-to-C++ migration; bitwise reproducibility; scientific software; molecular dynamics

---

## Abstract

MOPAC2016 is one of the most widely used semiempirical molecular-orbital programs, and its Fortran implementation has served as the reference for the PM7 and related NDDO methods for decades. We present MOPAC-CPP, a complete translation of the MOPAC2016 package (679 Fortran source files, ≈120,000 lines) into C++17. The translation is deliberately a faithful port, not a rewrite: names, dimensions, initial values, and 1-based indexing are preserved so every function remains auditable. On the nominal numerical path, the C++ build is bitwise identical to a local Intel-ifort rebuild of the official MOPAC2016 source (the verification build pair), verified on 18 small molecules in single-point mode, on 320/320 released-corpus inputs, and across CI, DRC, and multi-Hamiltonian (PM7/PM6/RM1) trajectories, with geometry-optimization trajectories identical within the verification build pair. The package also ships reproducibility gates that keep this fidelity while enabling parallelism: a fixed-bucket deterministic parallel Fock reduction (FBDR), a deterministic degenerate-basis choice for the optimizer (GAB), and an optional LAPACK exact-diagonalization mode for large systems, with a measured 1.85–2.5× performance improvement over the reference rebuild. MOPAC-CPP is released under the MIT license with its full test suite and verification scripts.

---

## Code metadata

| | |
|---|---|
| Current code version | v1.1.0 |
| Permanent link to code/repository used for this code version | `https://github.com/qq3185980737/MOPAC-CPP` (directory `mopac_cpp_public`) |
| Code Ocean compute capsule | (to be added at acceptance) |
| Legal Code License | MIT License |
| Code versioning system used | git |
| Software code languages, tools, and services used | C++17; MSVC 2022 (x64); Intel MKL 2025.3 (dsyevd/dgemm/BLAS kernels); optional: OpenCL (ICD loader + `opencl.lib`), compile-time `MOPAC_USE_OPENCL` |
| Compilation requirements, operating environments & dependencies | Windows x64; MSBuild/`mopac_cpp.sln`; Intel MKL (static link), SSE2 baseline; optional OpenCL runtime/ICD loader when the `MOPAC_USE_OPENCL` path is enabled |
| If available Link to developer documentation/manual | `https://github.com/qq3185980737/MOPAC-CPP/blob/main/mopac_cpp_public/README.md` |
| Support email for questions | liuyang123@stu.sjzu.edu.cn |

---

## 1. Motivation and significance

Semiempirical quantum-chemical methods (PM7 ¹, PM6 ², AM1 ³) remain a workhorse for systems where density-functional and wave-function methods are too expensive: large organic molecules, condensed-phase sampling, reaction-path screening, and molecular dynamics ⁴. Their reference implementation, MOPAC2016 ⁵⁻⁶, is a mature Fortran codebase (~120,000 lines, 679 source files) whose numerical semantics — including a particular pseudo-diagonalization SCF path and a fixed accumulation order — define the *de facto* values reported in thousands of published studies.

Scientific software in this state faces two well-known problems. First, legacy Fortran code is increasingly difficult to maintain, extend, and integrate into modern C++-based toolchains (AI/ML pipelines, visualization, workflow engines). Second, as reproducibility standards tighten ⁷,⁸, the community needs a defensible answer to the question: *when a program is re-expressed in another language, how do we know the numbers are still the numbers?* Prior work has addressed parallel non-determinism within a single code (e.g., FELTOR ⁹), used differential fault injection to test LLM-modernized code ¹⁰, and evaluated the correctness and reproducibility of LLM-generated parallel C/C++ translations at scale ¹¹; to our knowledge none reports a bitwise-verified full migration of a production quantum-chemistry package at this scale on the nominal path.

MOPAC-CPP addresses both problems. It provides (i) a complete, buildable C++17 port that preserves the original numerical semantics exactly on the default path, so existing inputs and expected outputs carry over unchanged; (ii) a public verification pipeline that certifies bitwise fidelity against the reference rebuild; and (iii) optional deterministic parallelism and an exact-diagonalization mode for users who need speed or robustness without sacrificing auditability. The software is aimed at computational chemists who need a modern, embeddable MOPAC, and at scientific-software engineers who need a worked, fully open example of a verified large-scale migration.

## 2. Software description

### 2.1 Architecture

The translation maps each Fortran module or function to a C++ translation unit of the same name (e.g., `compfg.cpp` ↔ `COMPFG`, `hcore.cpp` ↔ `HCORE`), preserving 1-based indexing through dedicated type aliases and fixed-size array layouts. Fortran common blocks are represented as plain-old-data structs, and the original call graph and data flow are left intact, which makes the code auditable line-by-line and keeps the behavior identical to the reference.

### 2.2 Functionality

The package implements the full MOPAC2016 feature set:

- **Hamiltonians:** PM7 (incl. PM7-TS), PM6 (incl. D3 and D3H4 corrections), PM6-DH+, AM1, PM3, RM1, MNDO, MNDO-d, and the linear-scaling **MOZYME** method for large systems.
- **Tasks:** single-point (1SCF), geometry optimization (EF, BFGS, L-BFGS), transition-state search (TS), IRC/DRC reaction paths, **molecular dynamics** (MOLDYN/DRC with velocity input), normal modes/thermochemistry, CI (MECI), ESP/polarizability, COSMO solvation, and property output (bonds, charges, dipoles, densities).
- **Input/output:** MOPAC keyword decks (full keyword set, documented in `KEYWORDS.md`), XYZ/PDB/Z-matrix geometries, and the classic MOPAC text output format.

### 2.3 Implementation and dependencies

- C++17, MSVC 2022 (x64), SSE2 baseline; build via the provided `mopac_cpp.sln` or MSBuild.
- Intel MKL (2025.3, static) for `dsyevd`/`dgemm`/BLAS kernels — the same kernel family used by the reference binary, which is part of why bitwise fidelity is attainable on the default path.
- Single-threaded by default (`OMP_NUM_THREADS=1; MKL_NUM_THREADS=1`) for bitwise-comparable runs.
- Optional OpenCL offload for the density rank-update kernels (compile-time `MOPAC_USE_OPENCL`, opt-in, experimental): when enabled, the `density_for_GPU` iopc=2/3 and 4/5 paths may run on an OpenCL GPU device; this path is excluded from the bitwise contract and falls back to the CPU path on device failure. Usage and measured results are documented in `docs/OPENCL_ACCEL.md`.

**Input/output.** MOPAC-CPP reads the classic MOPAC keyword deck: one or more keywords on the first line(s) (`PM7 1SCF GEO-OK`), a title line, and a geometry block in Z-matrix, XYZ, or PDB form. The full keyword set is implemented and documented in `KEYWORDS.md` (generated from the same keyword tables as the reference program), so existing MOPAC input decks run unchanged. Output is the classic MOPAC text format (header, `FINAL HEAT OF FORMATION`, population analysis, geometry, property blocks) so downstream parsers continue to work; trajectory and property exports are additionally available in plain XYZ and CSV form for modern toolchains. All files are plain text, which keeps the program scriptable and diff-friendly for regression testing.

### 2.4 Reproducibility gates (all default OFF)

The default binary follows the bitwise-identical-to-reference path exactly. Three optional gates are provided, all isolated behind dedicated keywords or environment variables so that the default path is never silently changed; a fourth, compile-time-only option (the OpenCL offload of §2.3) is likewise opt-in and excluded from the bitwise contract, but is not a keyword/environment-variable gate:

| Gate | Effect |
|---|---|
| `MOPAC_DSPAR` / `MOPAC_DSPAR_NB` | **FBDR** — fixed-bucket deterministic parallel Fock-matrix reduction; thread-count-independent, bitwise identical across 1/4/8/16 threads; `NB=1` equals the serial path |
| `MOPAC_GAB` | **GAB** — deterministic choice of basis inside degenerate eigenspaces of the EF optimizer, removing library-instance-dependent trajectory noise on symmetric molecules |
| `FULLDIAG` (keyword) / `MOPAC_FULL_DIAG` (env) | optional exact LAPACK (`dsyevd`) diagonalization for large systems; the fidelity–performance characterization is reported in the companion study ¹⁴ on exact versus pseudo diagonalization |

In prose: **FBDR** (`MOPAC_DSPAR` / `MOPAC_DSPAR_NB`) makes parallel Fock-matrix construction deterministic and thread-count-independent — the reduction is partitioned into fixed buckets whose partial sums are combined in a fixed order, and the output is byte-identical across 1, 4, 8, and 16 threads. **GAB** (`MOPAC_GAB`) removes library-instance-dependent trajectory noise inside degenerate eigenspaces of the optimizer by concentrating the projected gradient onto a deterministically chosen basis. **FULLDIAG** (`FULLDIAG` keyword / `MOPAC_FULL_DIAG` environment variable) swaps the SCF pseudo-diagonalization eigensolver for the LAPACK divide-and-conquer driver `dsyevd`, targeting large systems; its fidelity–performance trade-off is characterized in the companion study ¹⁴ on exact versus pseudo diagonalization.

### 2.5 Testing

Hundreds of unit tests (`test_*.cpp`) cover individual translated functions; a released regression corpus (349 inputs; 320/320 paired non-EF inputs bitwise identical, 317 of them 1SCF single points) verifies end-to-end behavior; `run_test/` contains executable integration cases, including molecular-dynamics runs.

## 3. Illustrative examples

### 3.1 Single-point PM7 energy (H₂O)

Input deck `h2o.mop` (XYZ):

```
PM7 1SCF
H2O single point

O   0.000000   0.000000   0.000000
H   0.957200   0.000000   0.000000
H  -0.238700   0.927000   0.000000
```

Run: `mopac.exe h2o.mop` → output `h2o.out` contains `FINAL HEAT OF FORMATION = -57.78974 kcal/mol` — bitwise identical to the reference MOPAC2016 rebuild on the verification build pair (measured with `OMP_NUM_THREADS=1; MKL_NUM_THREADS=1`).

### 3.2 Geometry optimization (EF, PM7)

```
PM7 EF GNORM=0.1
hexane optimization

 C ...
```

The EF trajectory and FINAL HEAT match the reference within the verification build pair; the trajectory is 1-ULP build-sensitive in the optimizer layer (documented), while the SCF layer and final converged energy are bitwise identical.

### 3.3 Molecular dynamics (velocity input)

```
PM7 MOLDYN T=300 ...
```

The DRC/MD path (including the MOPAC velocity-input `VELO=` handling) runs end-to-end and reproduces the reference trajectory; a trajectory-export module emits XYZ trajectories for downstream analysis.

### 3.4 Deterministic parallel run (FBDR)

```
set MOPAC_DSPAR=1
set OMP_NUM_THREADS=4
mopac.exe c80h162.mop
```

With FBDR enabled, the Fock-matrix reduction is partitioned into fixed buckets whose partial sums are combined in a deterministic order, so the output file is byte-identical whether the job runs on 1, 4, 8, or 16 threads. The same deck without `MOPAC_DSPAR` reproduces the reference rebuild bit-for-bit; with `MOPAC_DSPAR NB=1` the serial path is recovered exactly. This gives users a single switch between bitwise-parity with the reference and deterministic multi-threaded throughput.

## 4. Impact

**Fidelity.** On the nominal path the C++ build is bitwise identical to a local Intel-ifort rebuild of the official MOPAC2016 source: 18 small molecules (NH₃, HF, HCl, CO, N₂, CO₂, CH₂O, CH₃OH, H₂O, CH₄, SO₂, H₂S, PH₃, SiH₄, CH₃SH, CH₃Cl, O₂, CH₃) in single-point mode, with EF-optimized trajectories identical within the verification build pair; 320/320 released-corpus paired non-EF inputs with zero exceptions (317 of the 349 released decks are 1SCF single points); and CI, DRC, and multi-Hamiltonian (PM7/PM6/RM1) trajectories. The translation therefore contributes zero additional error at every level of the accompanying physical benchmark (NIST experimental heats of formation and GAMESS ¹² CCSD(T)/cc-pV{T,Q}Z CBS atomization energies; PM7's deviations from experiment are intrinsic parameterization errors, not migration errors).

**Performance.** On the same machine (AMD Ryzen 7 7840HS, 8 cores/16 threads, Windows x64) and workload (C₈₀H₁₆₂ PM7 RHF single point), MOPAC-CPP runs 1.85× faster than the reference rebuild with SSE2 + MKL and ≈2.5× with LTCG. With FBDR, deterministic parallelism adds 1.15–1.26× on C₈₀H₁₆₂ and 1.57× on C₁₂₀H₂₄₂ at 8 threads while remaining bitwise identical across thread counts. The optional FULLDIAG exact-diagonalization gate additionally targets large single points; its fidelity–performance trade-off is characterized in the companion study ¹⁴ on exact versus pseudo diagonalization.

**Reproducibility research.** The companion methodology paper ¹³ reports the ULP-attribution protocol and the discovery that EF-optimization trajectories are chaotic at the 1-ULP level even when the SCF layer is bitwise identical; MOPAC-CPP is the open, fully instrumented case study for that work, and its gates (FBDR, GAB) are reference implementations of the two deterministic methods described there.

**Audience value.** MOPAC-CPP lets existing MOPAC users switch to a modern C++ codebase with unchanged inputs and outputs, gives developers a clean embeddable library, and provides educators and researchers a complete, tested example of a verified 120k-line migration.

**Maintainability and community.** The repository ships the full build configuration, hundreds of unit tests, an end-to-end regression corpus (349 inputs; 320/320 paired non-EF inputs bitwise identical), and documentation (`README.md`, `KEYWORDS.md`, build and verification guides), so that a new contributor can rebuild, test, and audit the entire package without external services. The MIT license and the explicit file correspondence (the 679 Fortran source files are consolidated into 374 C++ translation units in the MSBuild project — a further 7 `.cpp` files are standalone test/utility sources — each documenting its source module) lower the barrier for domain experts to review or extend the code, and the reproducibility gates are isolated behind dedicated keywords and environment variables, keeping the default path stable across releases.

## 5. Conclusions

MOPAC-CPP is a complete, bitwise-verified C++17 translation of MOPAC2016, released open-source with tests, documentation, and reproducibility gates. It is a drop-in replacement for the Fortran program on the default path, a platform for extending MOPAC functionality in modern C++, and a worked reference for bitwise-reproducible scientific-software migration.

## Declaration of competing interest

The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Declaration of generative AI and AI-assisted technologies in the writing process

During the preparation of this work the authors used AI-assisted tools to assist with code translation, drafting, and language editing. After using these tools, the authors reviewed and edited the content as needed and take full responsibility for the content of the published article.

## Acknowledgements

This work has been financially supported by the Basic Scientific Research Project for Higher Education Institutions (Project No. LJ212610153014) of Liaoning Provincial Department of Education. The authors thank the MOPAC development community for the open-source release of the reference code, and the GAMESS project for the ab initio reference software used for the physical benchmarks.

## Code availability

The MOPAC-CPP source code, test suite, verification scripts, and documentation are publicly available under the MIT license at `https://github.com/qq3185980737/MOPAC-CPP` (directory `mopac_cpp_public`). MOPAC-CPP is a derivative work of MOPAC2016: the reference Fortran source is the official MOPAC2016 source code, archived in the openmopac/mopac-archive repository (BSD-3-Clause); the current open-source MOPAC release is numerically identical to MOPAC2016 (per the official FAQ) and is described in its JOSS release paper. The provenance and license status of the upstream source are documented in `UPSTREAM_NOTICE.md` in the repository. The C++ translation itself is licensed under the MIT license (see `LICENSE`). All 349 input decks with the 349 recorded C++ outputs and the 327 Fortran `.arc` summaries used in the verification reported here are included in the repository (`regression/`); a permanent Zenodo archive with a DOI will be registered upon acceptance. No commercial dependencies are required beyond the documented MSVC/Intel MKL toolchain (the optional OpenCL path additionally requires an OpenCL runtime/ICD loader).

## CRediT authorship contribution statement

**Yang Liu:** Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Writing – original draft, Visualization. **Xun Wang:** Conceptualization, Supervision, Resources, Funding acquisition, Writing – review & editing.


## References
[1] J. J. P. Stewart, Optimization of parameters for semiempirical methods VI: more modifications to the NDDO approximations and re-optimization of parameters, J. Mol. Model. 19 (2013) 1–32. https://doi.org/10.1007/s00894-012-1667-x
[2] J. J. P. Stewart, Optimization of parameters for semiempirical methods V: modification of NDDO approximations and application to 70 elements, J. Mol. Model. 13 (2007) 1173–1213. https://doi.org/10.1007/s00894-007-0233-4
[3] M. J. S. Dewar, E. G. Zoebisch, E. F. Healy, J. J. P. Stewart, Development and use of quantum mechanical molecular models. 76. AM1: a new general purpose quantum mechanical molecular model, J. Am. Chem. Soc. 107 (1985) 3902–3909. https://doi.org/10.1021/ja00299a024
[4] W. Thiel, Semiempirical quantum-chemical methods, Wiley Interdiscip. Rev.: Comput. Mol. Sci. 4 (2014) 145–157. https://doi.org/10.1002/wcms.1161
[5] J. J. P. Stewart, MOPAC: A semiempirical molecular orbital program, J. Comput.-Aided Mol. Des. 4 (1990) 1–103.
[6] J. J. P. Stewart, MOPAC2016, source code archive: https://github.com/openmopac/mopac-archive (directory 2016_MOPAC2016)
[7] R. D. Peng, Reproducible research in computational science, Science 334 (6060) (2011) 1226–1227. https://doi.org/10.1126/science.1213847
[8] G. K. Sandve, A. Nekrutenko, J. Taylor, T. Hovig, Ten simple rules for reproducible computational research, PLoS Comput. Biol. 9 (10) (2013) e1003285. https://doi.org/10.1371/journal.pcbi.1003285
[9] M. Wiesenberger, L. Einkemmer, M. Held, A. Gutierrez-Milla, X. Sáez, R. Iakymchuk, Reproducibility, accuracy and performance of the FELTOR code and library on parallel computer architectures, Comput. Phys. Commun. 238 (2019) 145–156.
[10] E. Coleman, Y. Shen, M. Sosonkina, P. Xu, Validating LLM-modernized scientific software through differential fault injection, arXiv:2608.14527 (2026).
[11] P. Valero-Lara, W. F. Godoy, J. Gonzalez, A. Huante, H. Gauthier-Chaparro, J. Gonzalez, LLM-driven Fortran-to-C/C++ portability for parallel scientific codes, in: Proceedings of the 2025 IEEE International Conference on eScience (eScience), IEEE, 2025.
[12] M. W. Schmidt et al., General atomic and molecular electronic structure system (GAMESS), J. Comput. Chem. 14 (1993) 1347–1363. https://doi.org/10.1002/jcc.540140112
[13] Y. Liu, X. Wang, Bitwise-reproducibility verification and deterministic numerical methods for scientific-software migration, with a full MOPAC2016 (Fortran to C++) case study, Comput. Phys. Commun. 2026, submitted.
[14] Y. Liu, X. Wang, Pseudo- versus exact diagonalization in semiempirical SCF calculations: a performance–robustness–fidelity trade-off study on the MOPAC2016 C++ translation, Int. J. Quantum Chem., submitted (2026).
