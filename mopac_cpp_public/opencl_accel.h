// opencl_accel.h - Optional OpenCL acceleration for MOPAC-CPP.
//
// PURPOSE
//   Offload BLAS-3-shaped hotspots (density rank-k updates / DGEMM) to an
//   OpenCL GPU when one is present and faster than the CPU path.  This is an
//   OPT-IN module: it is NOT compiled by default, and when enabled the GPU
//   path intentionally uses its own reduction order, so results may differ
//   from the bitwise reference at the 1-ULP level.  It must never be used on
//   the bitwise-verification path.
//
// ENABLE
//   - Define MOPAC_USE_OPENCL when compiling.
//   - Link against an OpenCL ICD loader (opencl.lib on Windows, -lOpenCL on
//     Linux/macOS) and include the OpenCL headers (CL/cl.h).
//   - At runtime, OpenCLAccel::init() finds a GPU platform/device and builds
//     the kernels; if no usable device exists the module reports
//     available()==false and callers silently fall back to the CPU path.
//
// SUITABLE HARDWARE
//   Discrete GPUs with strong FP64 (NVIDIA Tesla/Quadro/RTX with
//   cl_khr_fp64, AMD CDNA/RDNA Pro cards, Intel Arc Pro) or large-molecule
//   workloads where even consumer cards amortise the transfer cost.  iGPU
//   parts (e.g. AMD 780M) usually lose to the CPU for FP64 and are detected
//   but not preferred unless no other device exists.
//
// BITWISE NOTE
//   The kernels below compute the same mathematical quantities as the serial
//   CPU code but accumulate in a different order; on a given GPU the result is
//   deterministic (same kernel + same input => same bits), but it is NOT
//   guaranteed bit-identical to the reference build.  Treat this module as a
//   performance-only, opt-in experiment (like the FULLDIAG gate), never as
//   part of the L0 bitwise contract.
#pragma once

#include <string>
#include <vector>

class OpenCLAccel {
public:
    OpenCLAccel();
    ~OpenCLAccel();

    // Find and initialise a GPU OpenCL device.  If platform_substr is
    // non-empty, only platforms whose name contains it are considered.
    // Returns true when a device was found AND all kernels compiled.
    bool init(const std::string& platform_substr = "");

    bool available() const { return ok_; }
    const std::string& device_name() const { return dev_name_; }
    int compute_units() const { return cu_; }
    int max_clock_mhz() const { return clock_; }

    // DGEMM (row-major):  C[i*n+j] += sum_k A[i*n+k] * B[k*n+j].
    // C may already contain a value (beta=1 semantics); the kernel reads the
    // existing C element and adds the accumulated dot product.
    bool gemm(int n, const double* A, const double* B, double* C);

    // DSYRK upper (row-major A):  C(i,j) = alpha * sum_{k=0..kk-1} A[i*kk+k]*A[j*kk+k]
    // for 0<=i<=j<n.  Only the upper triangle is written (C[j*n+i] left
    // untouched); caller packs the triangle if needed.
    bool syrk_upper(int n, int kk, double alpha, const double* A, double* C);

    // Rank-2 update used by density_for_GPU iopc=2/3:
    //   C(i,j) = a * sum_k2 c2[i*k2+k]*c2[j*k2+k]  +  b * sum_k1 c1[i*k1+k]*c1[j*k1+k]
    // for 0<=i<=j<n (upper triangle; C[j*n+i] untouched).  c2/c1 are
    // n x k2 / n x k1 matrices in row-major layout.
    bool rank2_update(int n, int k2, const double* c2, double a,
                      int k1, const double* c1, double b, double* C);

    // Simple built-in self-test: run gemm/syrk/rank2 on small random inputs
    // and compare against a serial CPU reference.  Returns true on success.
    bool selftest();

private:
    struct Impl;
    Impl* impl_;

    bool ok_;
    std::string dev_name_;
    int cu_, clock_;
};
