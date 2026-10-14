# Optional OpenCL acceleration (opt-in, off by default)

`opencl_accel.{h,cpp}` provides an **optional, opt-in** OpenCL GPU offload for
the BLAS-3-shaped hotspots of MOPAC-CPP.  It is **not compiled by default**,
it is **not part of the bitwise (L0) verification contract**, and enabling it
changes the floating-point reduction order, so results on the GPU path are
deterministic per device but **not** guaranteed bit-identical to the reference
build.

## Why this exists

On hardware with a strong FP64 GPU (discrete NVIDIA/AMD/Intel Pro parts) and
for large-molecule workloads, offloading the density rank-k updates and the
symmetric rank update (`density_for_GPU` iopc=2/3 and iopc=4/5) can be faster
than the CPU path.  The module is shipped so that such machines can try it
without touching the default bitwise path.  On iGPUs such as the AMD 780M
measured here, FP64 throughput is too low to beat the CPU: if enabled on
such a part, `init()` succeeds and the program does use the GPU path (which
is slower than the CPU), so the module is aimed at discrete-GPU machines
with strong FP64.  On machines without a suitable device the module reports
`available()==false` and the CPU path is used unchanged.

## Files

| File | Role |
|---|---|
| `opencl_accel.h` | Public API (`OpenCLAccel` class): init / gemm / syrk / rank2 / selftest |
| `opencl_accel.cpp` | Implementation. Compiled with `MOPAC_USE_OPENCL`; without it a stub build is produced (module inert, `available()==false`) |
| `density_for_GPU.cpp` | Opt-in hooks: iopc=2/3 rank-2 update and iopc=4/5 DSYRK try the GPU first, fall back to CPU on device failure |
| `tests/test_opencl_accel.cpp` | Stand-alone self-test harness |

## Enabling

1. Install an OpenCL SDK (headers + ICD loader) and make sure a GPU driver
   with `cl_khr_fp64` is present.
2. Define `MOPAC_USE_OPENCL` when compiling **all** translation units that
   include the module (the accelerator object and `density_for_GPU.cpp`).
3. Link against the OpenCL ICD loader:
   - Windows/MSVC: `opencl.lib`
   - Linux/macOS: `-lOpenCL`
4. Build and run the self-test:
   ```
   # Windows (MSVC, with OCL_ROOT set to the SDK path)
   cl /O2 /DMOPAC_USE_OPENCL /I %OCL_ROOT%\include tests\test_opencl_accel.cpp opencl_accel.cpp /link %OCL_ROOT%\lib\x86_64\opencl.lib
   tests\test_opencl_accel.exe
   ```
   Expected output: `SELFTEST PASS` (gemm/syrk/rank2 relative error ~1e-16).
5. Rebuild the full program with `MOPAC_USE_OPENCL`; the SCF/density path
   will use the GPU automatically when a device is found.

## Behaviour

- Initialisation is lazy and one-time (`ocl_accel()` in `density_for_GPU.cpp`).
  If no usable FP64 GPU device exists, the module stays inert and the exact
  CPU reference path runs unchanged.
- GPU kernels compute only the upper triangle that the packed density needs
  (`i <= j`), matching the CPU serial semantics; the lower triangle of the
  scratch matrix is never read.
- On a GPU call failure (allocation / kernel error) the code falls back to the
  CPU path for that call, so a device hiccup cannot corrupt the density.

## Bitwise / determinism caveat

The kernels accumulate in a different order than the CPU loops, so a GPU-run
FINAL HEAT may differ from the reference build at the last few ULP (typical
relative deviation ~1e-16 per operation, growing with SCF iterations).  Use
this module **only** for performance experiments on fast-FP64 hardware, never
on the bitwise-verification corpus or in a reproducibility claim.

## Measured on this machine (AMD Radeon 780M, gfx1103, 6 CU, 2700 MHz)

Naive OpenCL DGEMM kernel vs CPU (MSVC /O2, single thread):

| N | CPU 1-thread | GPU OpenCL | GPU vs CPU1 |
|---|---|---|---|
| 256 | 15.5 GFLOPS | 8.8 GFLOPS | 0.56x |
| 512 | 10.8 GFLOPS | 9.1 GFLOPS | 0.84x |
| 1024 | 11.9 GFLOPS | 9.2 GFLOPS | 0.77x |

Conclusion on this iGPU: GPU FP64 is not faster than the CPU; the module is
retained for discrete-GPU machines where it can help.  Tiled/local-memory
kernels are a natural next step if a suitable machine becomes available.
