# Pseudo- versus exact diagonalization in semiempirical SCF calculations: a performance–robustness–fidelity trade-off study on the MOPAC2016 C++ translation

**Yang Liu^a, Xun Wang^a,\***
^a *Shenyang Jianzhu University, Shenyang, China*
*\*Corresponding author: wangxun71@aliyun.com (X. Wang); liuyang123@stu.sjzu.edu.cn (Y. Liu)*

**Keywords:** semiempirical methods; PM7; SCF convergence; diagonalization; LAPACK; reproducible computing

---

## Abstract

Semiempirical molecular-orbital methods such as PM7 solve the self-consistent-field (SCF) problem by repeated diagonalization of the Fock matrix. The reference MOPAC2016 implementation uses a traditional pairwise-rotation (Jacobi) eigensolver as its default SCF path, whereas modern LAPACK divide-and-conquer solvers (`dsyevd`) are substantially faster on dense matrices. We report a controlled study of this trade-off in MOPAC-CPP, a bitwise-faithful C++17 translation of MOPAC2016 in which the official pseudo-diagonalization path is preserved bit-for-bit and an exact-diagonalization mode (FULLDIAG) is added as an opt-in gate. For non-degenerate systems the two paths differ by only ≈2×10⁻⁴ kcal/mol in the single-point energy — at the level of floating-point rounding — and converge to mathematically equivalent optimized geometries (hexadecane ΔHf = −90.18 kcal/mol vs. NIST experimental −89.60). For a near-degenerate 122-atom alkane (C₄₀H₈₂), the pseudo-diagonalization path exhibits SCF/optimization divergence while the exact path converges stably. In a C₈₀H₁₆₂ PULAY-accelerated single point the exact path is 2.1–3.6× faster than the pseudo path. We frame the choice as a three-way trade-off between bitwise fidelity, convergence robustness, and performance, and give concrete guidance on when each path should be used.

## 1. Introduction

Semiempirical quantum-chemical methods [1] remain indispensable for large systems and for dynamics and screening applications where first-principles methods are prohibitively expensive [2]. Their computational core is the self-consistent-field (SCF) procedure: given a Fock matrix **F**, the density matrix is built from the occupied eigenvectors of **F**, and the process is iterated to self-consistency. The cost and numerical behavior of the SCF iteration therefore depend critically on the eigensolver.

MOPAC2016 [3,4], the standard open reference implementation of the PM7 [5], PM6 [6], and AM1 [7] family of methods, uses a *pseudo-diagonalization* scheme based on pairwise Jacobi rotations as its default SCF eigensolver. This choice is historically motivated: for the small matrices typical of semiempirical calculations (a few hundred orbitals), a purpose-built Jacobi sweeper avoids the overhead of a general LAPACK call and preserves a fixed, reproducible accumulation order — a property that has made MOPAC outputs stable and comparable across decades of published work.

Modern LAPACK divide-and-conquer eigensolvers (`dsyevd`) are, however, significantly faster on dense symmetric matrices and are numerically robust, including in the presence of (near-)degenerate eigenvalues, where the outcome of any eigensolver becomes basis-dependent. For a modernized, reproducible re-implementation of MOPAC, the natural question is whether the traditional pseudo-diagonalization path is still the right default, or whether an exact-diagonalization path can deliver better performance and robustness without compromising the numerical contract that users rely on.

We answer this question in the context of MOPAC-CPP, a complete C++17 translation of MOPAC2016 in which the official numerical path — including the pseudo-diagonalization eigensolver and its accumulation order — is preserved bit-for-bit against the reference binary [8]. We add an opt-in exact-diagonalization gate (FULLDIAG, activated by an input keyword or an environment variable) and characterize the resulting trade-off along three axes: fidelity (how much the energy changes), robustness (convergence on near-degenerate systems), and performance (wall-clock time). The default path is unchanged, so all bitwise-fidelity guarantees of the translation are retained unless the gate is explicitly enabled.

## 2. Methods

### 2.1 Semiempirical SCF and the Fock-matrix eigenproblem

In the NDDO-based semiempirical methods, the RHF SCF cycle solves, at each iteration,

**F** C = **S** C **ε**,  (1)

where **F** is the Fock matrix, **S** the overlap matrix, **C** the coefficient matrix, and **ε** the diagonal orbital energies. The density matrix is **P** = C_occ C_occᵀ, where C_occ contains the occupied columns of **C**, and the cycle is iterated until the density and energy are self-consistent (subject to the PULAY/DIIS accelerator when enabled [9]).

For closed-shell calculations on organic molecules with n ≈ 50–500 orbitals, the diagonalization of **F** is a significant fraction of the per-iteration cost, and its numerical result — in particular the *sign* and *composition* of eigenvectors inside near-degenerate subspaces — feeds directly into the density matrix and hence into the next Fock matrix.

### 2.2 Pseudo-diagonalization in MOPAC (pairwise Jacobi rotations)

The MOPAC default SCF eigensolver performs classical Jacobi sweeps: it annihilates off-diagonal elements by applying a sequence of plane rotations,

**F** ← **R**(p,q)ᵀ **F** **R**(p,q),

one pair (p,q) at a time, until the off-diagonal norm is below a threshold. This is a *pseudo-diagonalization* in the sense that the sweeps are interleaved with the SCF iteration: the eigensolver is invoked on the current Fock matrix and only the matrix elements that change are re-rotated in later iterations, an optimization that exploits the small incremental change of **F** between SCF cycles. The fixed rotation schedule gives the method a deterministic accumulation order, which is why the official MOPAC2016 binary is exactly reproducible run-to-run on a given build.

The weakness of the pseudo path is exposed when eigenvalues are (near-)degenerate: the rotation-based sweeper may struggle to separate a near-degenerate pair cleanly, and the resulting eigenvector composition can fluctuate, leading to slow or non-convergent SCF cycles in pathological cases.

### 2.3 Exact diagonalization (LAPACK dsyevd)

The exact path uses the LAPACK divide-and-conquer symmetric eigensolver `dsyevd` (via Intel MKL, the same kernel library family used by the reference build). `dsyevd` computes the full eigensystem of the dense symmetric matrix **F** with documented error bounds, handles (near-)degenerate spectra robustly, and is heavily optimized (BLAS-3 kernels, cache-blocked divide-and-conquer). For n in the hundreds it is typically faster than the Jacobi sweeper despite the higher asymptotic flop count, because the sweeper's O(n³) rotation schedule is not cache-optimized.

### 2.4 The FULLDIAG gate in MOPAC-CPP

MOPAC-CPP is a 1-based, line-auditable C++17 translation of the 679-file MOPAC2016 source [8]. Its default SCF path is bitwise identical to the official binary (verified on 320/320 released-regression inputs and on 18 small molecules in single-point mode, with geometry-optimization trajectories identical within the verification build pair). The exact-diagonalization mode is exposed as an opt-in gate:

- input keyword: `PM7 FULLDIAG 1SCF` (or any deck containing `FULLDIAG`), or
- environment variable: `set MOPAC_FULL_DIAG=1`

Both are equivalent. When the gate is off — the default — the binary follows the official pseudo-diagonalization path exactly. When it is on, the Fock-matrix eigensolver is replaced by `dsyevd` at each SCF iteration; everything else (integral evaluation, Fock assembly, density construction, DIIS/PULAY accelerator, EF optimizer) is unchanged.

## 3. Results

### 3.1 Fidelity of the default path

With the FULLDIAG gate off, MOPAC-CPP reproduces the official MOPAC2016 binary bit-for-bit on the nominal path: 320/320 released-regression single-point inputs with zero exceptions; FINAL HEAT of formation bitwise identical on 18 small molecules in single-point mode, and on the same molecules in EF-optimized mode within the verification build pair [8]. All fidelity statements below refer to the *default* path; enabling FULLDIAG is an explicit, documented departure from bitwise parity.

### 3.2 Energy difference on non-degenerate systems

On a C₈₀H₁₆₂ PM7 RHF single point (PULAY-accelerated, 482 basis functions; C: 4 AOs × 80, H: 1 AO × 162), the single-point energy computed with FULLDIAG differs from the pseudo-diagonalization result by ≈2×10⁻⁴ kcal/mol — at the level of floating-point rounding and far below the PM7 method error (≈±5–10 kcal/mol against experiment). The optimized geometries and FINAL HEAT are mathematically equivalent: for hexadecane (C₁₆H₃₄), FULLDIAG gives ΔHf = −90.18 kcal/mol vs. the NIST experimental value −89.60 kcal/mol (PM7 method error −0.58, within the method's expected accuracy band).

### 3.3 Near-degenerate robustness: C₄₀H₈₂

The 122-atom alkane C₄₀H₈₂ provides a stress case. On the pseudo-diagonalization path the geometry optimization diverges: the SCF cycle fails to converge stably and the EF trajectory breaks down. On the exact path the same input converges smoothly to a stable optimized geometry. This case demonstrates that the exact eigensolver's robust treatment of (near-)degenerate spectra translates directly into SCF/optimization robustness — the regime where the traditional pseudo path is weakest.

### 3.4 Performance

On C₈₀H₁₆₂ PULAY single point (the workload of §3.2), the exact path is **2.1×** faster (SSE2 build) and **3.6×** faster (LTCG build) than the pseudo path, reflecting `dsyevd`'s cache-blocked divide-and-conquer efficiency versus the Jacobi sweeper's O(n³) rotation schedule. The gap grows with matrix size, as the pseudo path's per-iteration rotation work scales poorly relative to BLAS-3 kernels.

### 3.5 Summary of the trade-off

| System | Energy diff (exact vs pseudo) | Convergence | Speed (exact vs pseudo) |
|---|---|---|---|
| C₈₀H₁₆₂ (non-degenerate) | ≈2×10⁻⁴ kcal/mol | equivalent | 2.1× (SSE2) / 3.6× (LTCG) |
| C₁₆H₃₄ (optimization) | equivalent geometry, ΔHf −90.18 vs NIST −89.60 | equivalent | — |
| C₄₀H₈₂ (near-degenerate) | — | pseudo diverges; exact converges | — |

## 4. Discussion

### 4.1 When to sacrifice bitwise fidelity

The results define two clean regimes.

- **Default path (pseudo-diagonalization)** — required when bitwise reproducibility against the official MOPAC2016 binary is the contract: regression testing against archived results, verification of a translation, or any workflow where published reference values must be matched exactly.
- **Exact path (FULLDIAG)** — preferred when (i) the system is large and the diagonalization dominates the SCF cost (2–4× speedup), or (ii) the spectrum is (near-)degenerate and the pseudo path diverges. The energy cost is at most ≈2×10⁻⁴ kcal/mol on non-degenerate systems, which is negligible against PM7's intrinsic error; the geometry/energy end-points are mathematically equivalent.

The gate design — opt-in, keyword or environment variable, default unchanged — is deliberate: it lets a reproducible-software project offer a performance/robustness upgrade without ever silently changing the numbers a user gets by default.

### 4.2 Implications for reproducible computing

The study is a worked example of a more general principle: in a faithfully translated scientific code, a numerical-algorithm substitution (here, the eigensolver) can be characterized *exactly* — bitwise parity of the default path, and a quantified, negligible deviation of the alternative. This turns an otherwise fuzzy "does the new solver change my results?" question into a precise, published statement. The same protocol transfers to any other kernel substitution (linear solvers, integral screening, reduction order) in legacy scientific software.

### 4.3 Limitations

- The near-degenerate robustness case (C₄₀H₈₂) is reported for a specific 122-atom alkane; a broader survey of pathological spectra is left to future work.
- Performance numbers are single-platform (Windows x64, MSVC 2022, Intel MKL 2025.3); relative gains are expected to carry to other platforms, but absolute times do not.
- FULLDIAG changes the numerical path; users requiring bit-for-bit parity with archived MOPAC2016 outputs must keep the default gate off.

## 5. Conclusions

We characterized the pseudo- versus exact-diagonalization trade-off in semiempirical SCF calculations using MOPAC-CPP, a bitwise-faithful C++17 translation of MOPAC2016. The exact LAPACK path costs ≈2×10⁻⁴ kcal/mol on non-degenerate systems, delivers mathematically equivalent optimized geometries (NIST-verified), fixes near-degenerate SCF divergence on a 122-atom alkane where the pseudo path fails, and is 2.1–3.6× faster on a large PULAY-accelerated single point. The opt-in gate design preserves the default path's bitwise fidelity, providing a template for performance upgrades in reproducible scientific software.

## Acknowledgements

This work has been financially supported by the Basic Scientific Research Project for Higher Education Institutions (Project No. LJ212610153014) of Liaoning Provincial Department of Education.

## Data and code availability

The MOPAC-CPP source code, the FULLDIAG implementation, the test suite, and all input decks used in this study are available in the public repository (MIT license) at https://github.com/qq3185980737/test_rep (mopac_cpp_public).

## References

[1] W. Thiel, Semiempirical quantum-chemical methods, Wiley Interdiscip. Rev.: Comput. Mol. Sci. 4 (2) (2014) 145–157. https://doi.org/10.1002/wcms.1161
[2] J. J. P. Stewart, Application of the PM6 method to modeling the solid state, J. Mol. Model. 14 (2008) 499–535.
[3] J. J. P. Stewart, MOPAC: A semiempirical molecular orbital program, J. Comput.-Aided Mol. Des. 4 (1990) 1–103.
[4] MOPAC2016 DOI: https://doi.org/10.5281/zenodo.6511958
[5] J. J. P. Stewart, Optimization of parameters for semiempirical methods VI: more modifications to the NDDO approximations and re-optimization of parameters, J. Mol. Model. 19 (2013) 1–32. https://doi.org/10.1007/s00894-012-1667-x
[6] J. J. P. Stewart, Optimization of parameters for semiempirical methods V: modification of NDDO approximations and application to 70 elements, J. Mol. Model. 13 (12) (2007) 1173–1213. https://doi.org/10.1007/s00894-007-0233-4
[7] M. J. S. Dewar, E. G. Zoebisch, E. F. Healy, J. J. P. Stewart, Development and use of quantum mechanical molecular models. 76. AM1: a new general purpose quantum mechanical molecular model, J. Am. Chem. Soc. 107 (13) (1985) 3902–3909. https://doi.org/10.1021/ja00299a024
[8] Y. Liu, X. Wang, Bitwise-Reproducibility Verification and Deterministic Numerical Methods for Scientific-Software Migration, with a Full MOPAC2016 (Fortran to C++) Case Study, Comput. Phys. Commun. (submitted).
[9] P. Pulay, Convergence acceleration of iterative sequences. The case of SCF iteration, Chem. Phys. Lett. 73 (2) (1980) 393–398. https://doi.org/10.1016/0009-2614(80)80396-4
