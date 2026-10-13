# MOPAC-CPP: a bitwise-reproducible C++17 translation of the MOPAC2016 semiempirical molecular-orbital package

**SoftwareX — Original Software Publication (draft)**

**Authors:** Yang Liu^a, Xun Wang^b,*

^a ^b Shenyang Jianzhu University, Shenyang, China
*Corresponding author: wangxun71@aliyun.com (X. Wang); liuyang123@stu.sjzu.edu.cn (Y. Liu)

**Keywords:** semiempirical quantum chemistry; PM7; Fortran-to-C++ migration; bitwise reproducibility; scientific software; molecular dynamics

---

## Abstract

MOPAC2016 is one of the most widely used semiempirical molecular-orbital programs, and its Fortran implementation has served as the reference for the PM7 and related NDDO methods for decades. We present MOPAC-CPP, a complete translation of the MOPAC2016 package (679 Fortran source files, ≈120,000 lines) into C++17. The translation is deliberately a faithful port, not a rewrite: names, dimensions, initial values, and 1-based indexing are preserved so every function remains auditable. On the nominal numerical path, the C++ build is bitwise identical to the official MOPAC2016 executable, verified on 18 small molecules in single-point mode, on 320/320 released-regression inputs, and across CI, DRC, and multi-Hamiltonian (PM7/PM6/RM1) trajectories, with geometry-optimization trajectories identical within the verification build pair. The package also ships reproducibility gates that keep this fidelity while enabling parallelism: a fixed-bucket deterministic parallel Fock reduction (FBDR), a deterministic degenerate-basis choice for the optimizer (GAB), and an optional LAPACK exact-diagonalization mode for large systems, and a first-line performance improvement of 1.85–2.5× over the reference binary. MOPAC-CPP is released under the MIT license with its full test suite and verification scripts.

---

## Code metadata

| | |
|---|---|
| Current code version | v1.1.0 |
| Permanent link to code/repository used for this code version | https://github.com/qq3185980737/test_rep (directory `mopac_cpp_public`) |
| Code Ocean compute capsule | (to be added at acceptance) |
| Legal Code License | MIT License |
| Code versioning system used | git |
| Software code languages, tools, and services used | C++17; MSVC 2022 (x64); Intel MKL 2025.3 (dsyevd/dgemm/BLAS kernels) |
| Compilation requirements, operating environments & dependencies | Windows x64; MSBuild/`mopac_cpp.sln`; Intel MKL (static link), SSE2 baseline |
| If available Link to developer documentation/manual | https://github.com/qq3185980737/test_rep/blob/main/README.md |
| Support email for questions | liuyang123@stu.sjzu.edu.cn |

---

## 1. Motivation and significance

Semiempirical quantum-chemical methods (PM7 ¹, PM6 ², AM1 ³) remain a workhorse for systems where density-functional and wave-function methods are too expensive: large organic molecules, condensed-phase sampling, reaction-path screening, and molecular dynamics ⁴. Their reference implementation, MOPAC2016 ⁵,⁶, is a mature Fortran codebase (~120,000 lines, 679 source files) whose numerical semantics — including a particular pseudo-diagonalization SCF path and a fixed accumulation order — define the *de facto* values reported in thousands of published studies.

Scientific software in this state faces two well-known problems. First, legacy Fortran code is increasingly difficult to maintain, extend, and integrate into modern C++-based toolchains (AI/ML pipelines, visualization, workflow engines). Second, as reproducibility standards tighten ⁷,⁸, the community needs a defensible answer to the question: *when a program is re-expressed in another language, how do we know the numbers are still the numbers?* Prior work has addressed parallel non-determinism within a single code (e.g., FELTOR ⁹) or used differential fault injection to test LLM-modernized code ¹⁰, but no published work reports a bitwise-verified full migration of a production quantum-chemistry package at this scale.

MOPAC-CPP addresses both problems. It provides (i) a complete, buildable C++17 port that preserves the original numerical semantics exactly on the default path, so existing inputs and expected outputs carry over unchanged; (ii) a public verification pipeline that certifies bitwise fidelity against the official binary; and (iii) optional deterministic parallelism and an exact-diagonalization mode for users who need speed or robustness without sacrificing auditability. The software is aimed at computational chemists who need a modern, embeddable MOPAC, and at scientific-software engineers who need a worked, fully open example of a verified large-scale migration.

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

**Input/output.** MOPAC-CPP reads the classic MOPAC keyword deck: one or more keywords on the first line(s) (`PM7 1SCF GEO-OK`), a title line, and a geometry block in Z-matrix, XYZ, or PDB form. The full keyword set is implemented and documented in `KEYWORDS.md` (generated from the same keyword tables as the reference program), so existing MOPAC input decks run unchanged. Output is the classic MOPAC text format (header, `FINAL HEAT OF FORMATION`, population analysis, geometry, property blocks) so downstream parsers continue to work; trajectory and property exports are additionally available in plain XYZ and CSV form for modern toolchains. All files are plain text, which keeps the program scriptable and diff-friendly for regression testing.

### 2.4 Reproducibility gates (all default OFF)

The default binary follows the bitwise-official path exactly. Three optional gates are provided:

| Gate | Effect |
|---|---|
| `MOPAC_DSPAR` / `MOPAC_DSPAR_NB` | **FBDR** — fixed-bucket deterministic parallel Fock-matrix reduction; thread-count-independent, bitwise identical across 1/4/8/16 threads; `NB=1` equals the serial path |
| `MOPAC_GAB` | **GAB** — deterministic choice of basis inside degenerate eigenspaces of the EF optimizer, removing library-instance-dependent trajectory noise on symmetric molecules |
| `FULLDIAG` (keyword) / `MOPAC_FULL_DIAG` (env) | optional exact LAPACK (`dsyevd`) diagonalization: ≈2.1–3.6× faster 1SCF on large alkanes and robust to near-degenerate SCF convergence; converges to mathematically equivalent optimized geometries (NIST-verified) while the default pseudo-diagonalization path keeps bitwise fidelity |

### 2.5 Testing

Hundreds of unit tests (`test_*.cpp`) cover individual translated functions; a released-regression corpus (320 inputs) verifies end-to-end single-point behavior; `run_test/` contains executable integration cases, including molecular-dynamics runs.

## 3. Illustrative examples

### 3.1 Single-point PM7 energy (H₂O)

Input deck `h2o.mop`:

```
PM7 1SCF GEO-OK
H2O single point

 O
 H 1 r
 H 1 r 2 a
 r=0.957 a=104.5
```

Run: `mopac.exe h2o.mop` → output `h2o.out` contains `FINAL HEAT OF FORMATION = -57.799 kcal/mol` — bitwise identical to the official MOPAC2016 binary on the verification build pair.

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

With FBDR enabled, the Fock-matrix reduction is partitioned into fixed buckets whose partial sums are combined in a deterministic order, so the output file is byte-identical whether the job runs on 1, 4, 8, or 16 threads. The same deck without `MOPAC_DSPAR` reproduces the official binary bit-for-bit; with `MOPAC_DSPAR NB=1` the serial path is recovered exactly. This gives users a single switch between bitwise-parity with the reference and deterministic multi-threaded throughput.

## 4. Impact

**Fidelity.** On the nominal path the C++ build is bitwise identical to the official MOPAC2016 executable: 18 small molecules (NH₃, HF, HCl, CO, N₂, CO₂, CH₂O, CH₃OH, H₂O, CH₄, SO₂, H₂S, PH₃, SiH₄, CH₃SH, CH₃Cl, O₂, CH₃) in single-point mode, with EF-optimized trajectories identical within the verification build pair; 320/320 released-regression single-point inputs with zero exceptions; and CI, DRC, and multi-Hamiltonian (PM7/PM6/RM1) trajectories. The translation therefore contributes zero additional error at every level of the accompanying physical benchmark (NIST experimental heats of formation and GAMESS ¹¹ CCSD(T)/cc-pV{T,Q}Z CBS atomization energies; PM7's deviations from experiment are intrinsic parameterization errors, not migration errors).

**Performance.** On the same machine and workload (C₈₀H₁₆₂ PM7 RHF single point), MOPAC-CPP runs 1.85× faster than the official binary with SSE2 + MKL and ≈2.5× with LTCG. With FBDR, deterministic parallelism adds 1.15–1.26× on C₈₀H₁₆₂ and 1.57× on C₁₂₀H₂₄₂ at 8 threads while remaining bitwise identical across thread counts. On large near-degenerate systems the optional FULLDIAG path is 2.1–3.6× faster in 1SCF and fixes SCF-divergence cases (e.g., a C₄₀H₈₂ optimization that diverges on the pseudo-diagonalization path but converges on the full path).

**Reproducibility research.** The companion methodology paper (submitted) reports the ULP-attribution protocol and the discovery that EF-optimization trajectories are chaotic at the 1-ULP level even when the SCF layer is bitwise identical; MOPAC-CPP is the open, fully instrumented case study for that work, and its gates (FBDR, GAB) are reference implementations of the two deterministic methods described there.

**Audience value.** MOPAC-CPP lets existing MOPAC users switch to a modern C++ codebase with unchanged inputs and outputs, gives developers a clean embeddable library, and provides educators and researchers a complete, tested example of a verified 120k-line migration.

**Maintainability and community.** The repository ships the full build configuration, hundreds of unit tests, an end-to-end regression corpus (320 released inputs), and documentation (`README.md`, `KEYWORDS.md`, build and verification guides), so that a new contributor can rebuild, test, and audit the entire package without external services. The MIT license and the 1:1 Fortran-to-C++ file correspondence (each translated unit documents its source module) lower the barrier for domain experts to review or extend the code, and the reproducibility gates are isolated behind dedicated keywords and environment variables, keeping the default path stable across releases.

## 5. Conclusions

MOPAC-CPP is a complete, bitwise-verified C++17 translation of MOPAC2016, released open-source with tests, documentation, and reproducibility gates. It is a drop-in replacement for the Fortran program on the default path, a platform for extending MOPAC functionality in modern C++, and a worked reference for bitwise-reproducible scientific-software migration.

## Declaration of competing interest

The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Acknowledgements

This work has been financially supported by the Basic Scientific Research Project for Higher Education Institutions (Project No. LJ212610153014) of Liaoning Provincial Department of Education. The authors thank the MOPAC development community for the open-source release of the reference code, and the GAMESS project for the ab initio reference software used for the physical benchmarks.

## Code availability

The MOPAC-CPP source code, test suite, verification scripts, and documentation are publicly available under the MIT license at https://github.com/qq3185980737/test_rep (directory `mopac_cpp_public`). All input decks and recorded outputs used in the verification reported here are included in the repository, and the repository is archived on Zenodo with a DOI (https://doi.org/10.5281/zenodo.6511958 is the DOI of the reference MOPAC2016 archive; the MOPAC-CPP archive DOI will be added at acceptance). No commercial dependencies are required beyond the documented MSVC/Intel MKL toolchain.

## CRediT authorship contribution statement

**Yang Liu:** Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Writing – original draft, Visualization. **Xun Wang:** Conceptualization, Supervision, Resources, Funding acquisition, Writing – review & editing.

## References
[1] J. J. P. Stewart, Optimization of parameters for semiempirical methods VI: more modifications to the NDDO approximations and re-optimization of parameters, J. Mol. Model. 19 (2013) 1–32. https://doi.org/10.1007/s00894-012-1667-x
[2] J. J. P. Stewart, Optimization of parameters for semiempirical methods V: modification of NDDO approximations and application to 70 elements, J. Mol. Model. 13 (2007) 1173–1213. https://doi.org/10.1007/s00894-007-0233-4
[3] M. J. S. Dewar, E. G. Zoebisch, E. F. Healy, J. J. P. Stewart, Development and use of quantum mechanical molecular models. 76. AM1: a new general purpose quantum mechanical molecular model, J. Am. Chem. Soc. 107 (1985) 3902–3909. https://doi.org/10.1021/ja00299a024
[4] W. Thiel, Semiempirical quantum-chemical methods, Wiley Interdiscip. Rev.: Comput. Mol. Sci. 4 (2014) 145–157. https://doi.org/10.1002/wcms.1161
[5] J. J. P. Stewart, MOPAC: A semiempirical molecular orbital program, J. Comput.-Aided Mol. Des. 4 (1990) 1–103.
[6] MOPAC2016 DOI: https://doi.org/10.5281/zenodo.6511958
[7] R. D. Peng, Reproducible research in computational science, Science 334 (6060) (2011) 1226–1227. https://doi.org/10.1126/science.1213847
[8] G. K. Sandve, A. Nekrutenko, J. Taylor, T. Hovig, Ten simple rules for reproducible computational research, PLoS Comput. Biol. 9 (10) (2013) e1003285. https://doi.org/10.1371/journal.pcbi.1003285
[9] M. Wiesenberger, L. Einkemmer, M. Held, A. Gutierrez-Milla, X. Sáez, R. Iakymchuk, Reproducibility, accuracy and performance of the FELTOR code and library on parallel computer architectures, Comput. Phys. Commun. 238 (2019) 145–156.
[10] Yuan et al., Differential fault-injection validation of LLM-modernized scientific software (GAMESS), 2026. https://hub.baai.ac.cn/paper/27140009-4896-40d9-863e-ba267771d133
[11] M. W. Schmidt et al., General atomic and molecular electronic structure system (GAMESS), J. Comput. Chem. 14 (1993) 1347–1363. https://doi.org/10.1002/jcc.540140112
