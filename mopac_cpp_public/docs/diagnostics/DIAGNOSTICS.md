# Probe-development diagnostics archive (2026-10-14)

This directory preserves diagnostic outputs produced while developing the
probe-instrumented comparison between the C++ translation and the Fortran
reference rebuild. They document one intermediate comparison artifact and are
retained for auditability; they do **not** describe the final verified state of
the code.

## hcore comparison artifact (07:43–07:48 local time)

Files: `_hcore_out.txt`, `_geo_lines.txt`, `_coord_lines.txt`,
`_heat_lines.txt`, `_iter_heat_cmp.txt`.

Observed discrepancy at the time of the probe run:

- CORE-CORE REPULSION: Fortran rebuild `14373.76725 eV` vs C++ `14759.50729 eV`
  (a difference of 385.74 eV).
- Geometry block: the Fortran output shows the **input** Cartesian coordinates
  (regular chain: 0.000, 0.889, 1.778, 2.667, 3.556, 4.445 Å), while the C++
  output shows **EF-optimized** coordinates (−2.084 … 3.930 Å).

Diagnosis: the two runs were compared at **different trajectory states**
(input geometry in the Fortran output versus an optimized geometry in the C++
output), i.e., the probe captured output blocks that do not correspond to the
same point of the computation. The 385.74 eV core-core difference is a
state/comparison artifact of the probe setup, not a translation error in the
hcore energy terms.

Resolution history: this artifact belongs to the earlier probe-development
phase. The MOZYME hcore/Fock-assembly alignment fixes (commit chain preceding
`0b58f9b`; see git history) closed the underlying probe-placement issues, and
the current verification state is described in the manuscripts
(`paper_cpc_manuscript.md` §3, `softwarex_manuscript.md` §2.5, and
`ijqc_manuscript.md` §3).

## Re-use

The files are retained verbatim (unmodified) for traceability. When comparing
C++ and Fortran outputs, always compare the same run mode (e.g., both 1SCF, or
both fully converged EF) and the same output block, as established by the
verification pipeline in §2 of the CPC manuscript.
