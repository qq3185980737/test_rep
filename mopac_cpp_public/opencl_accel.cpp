// opencl_accel.cpp - Optional OpenCL acceleration for MOPAC-CPP (see header).
#include "opencl_accel.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#endif

// ---- OpenCL is only compiled when the module is enabled ----
#ifdef MOPAC_USE_OPENCL

#include <CL/cl.h>

namespace {

const char* err_str(cl_int e) {
    switch (e) {
        case CL_SUCCESS: return "CL_SUCCESS";
        case CL_DEVICE_NOT_FOUND: return "CL_DEVICE_NOT_FOUND";
        case CL_DEVICE_NOT_AVAILABLE: return "CL_DEVICE_NOT_AVAILABLE";
        case CL_COMPILER_NOT_AVAILABLE: return "CL_COMPILER_NOT_AVAILABLE";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
        case CL_OUT_OF_RESOURCES: return "CL_OUT_OF_RESOURCES";
        case CL_OUT_OF_HOST_MEMORY: return "CL_OUT_OF_HOST_MEMORY";
        case CL_BUILD_PROGRAM_FAILURE: return "CL_BUILD_PROGRAM_FAILURE";
        case CL_INVALID_KERNEL_NAME: return "CL_INVALID_KERNEL_NAME";
        default: return "OPENCL_ERR";
    }
}

// ---------------------------------------------------------------- kernels
const char* KERNEL_SRC = R"CLC(
// C[i*n+j] += sum_k A[i*n+k]*B[k*n+j]   (row-major)
__kernel void k_gemm(__global const double* A, __global const double* B,
                     __global double* C, int n) {
    int i = get_global_id(0), j = get_global_id(1);
    if (i >= n || j >= n) return;
    double s = 0.0;
    for (int k = 0; k < n; ++k) s += A[i*n+k] * B[k*n+j];
    C[i*n+j] += s;
}

// upper-triangle DSYRK:  C(i,j) = alpha * sum_k A[i*kk+k]*A[j*kk+k]
__kernel void k_syrk_u(__global const double* A, __global double* C,
                       int n, int kk, double alpha) {
    int i = get_global_id(0), j = get_global_id(1);
    if (i >= n || j > i) return;
    double s = 0.0;
    for (int k = 0; k < kk; ++k) s += A[i*kk+k] * A[j*kk+k];
    C[i*n+j] = alpha * s;
}

// rank-2 upper update:  C(i,j) = a*sum_k2 c2[i*k2+k]*c2[j*k2+k] + b*sum_k1 c1[i*k1+k]*c1[j*k1+k]
__kernel void k_rank2(__global const double* c2, int k2, double a,
                      __global const double* c1, int k1, double b,
                      __global double* C, int n) {
    int i = get_global_id(0), j = get_global_id(1);
    if (i >= n || j > i) return;
    double s = 0.0;
    for (int k = 0; k < k2; ++k) s += c2[i*k2+k] * c2[j*k2+k];
    s *= a;
    double s1 = 0.0;
    for (int k = 0; k < k1; ++k) s1 += c1[i*k1+k] * c1[j*k1+k];
    C[i*n+j] = s + b * s1;
}
)CLC";

struct OclImpl {
    cl_platform_id platform = nullptr;
    cl_device_id device = nullptr;
    cl_context ctx = nullptr;
    cl_command_queue q = nullptr;
    cl_program prog = nullptr;
    cl_kernel k_gemm = nullptr, k_syrk = nullptr, k_rank2 = nullptr;
    ~OclImpl() {
        if (k_gemm) clReleaseKernel(k_gemm);
        if (k_syrk) clReleaseKernel(k_syrk);
        if (k_rank2) clReleaseKernel(k_rank2);
        if (prog) clReleaseProgram(prog);
        if (q) clReleaseCommandQueue(q);
        if (ctx) clReleaseContext(ctx);
    }
};

} // namespace

struct OpenCLAccel::Impl {
    OclImpl o;
    bool built = false;
};

OpenCLAccel::OpenCLAccel() : impl_(new Impl), ok_(false), cu_(0), clock_(0) {}
OpenCLAccel::~OpenCLAccel() { delete impl_; }

bool OpenCLAccel::init(const std::string& platform_substr) {
    cl_int err;
    cl_uint nplat = 0;
    if (clGetPlatformIDs(0, nullptr, &nplat) != CL_SUCCESS || nplat == 0)
        return false;

    std::vector<cl_platform_id> plats(nplat);
    clGetPlatformIDs(nplat, plats.data(), nullptr);

    // Prefer a GPU device; fall back to any device if nothing else.
    cl_device_id chosen = nullptr;
    for (cl_uint p = 0; p < nplat; ++p) {
        if (!platform_substr.empty()) {
            char pname[256] = {0};
            clGetPlatformInfo(plats[p], CL_PLATFORM_NAME, sizeof(pname), pname, nullptr);
            if (!strstr(pname, platform_substr.c_str())) continue;
        }
        cl_uint ndev = 0;
        if (clGetDeviceIDs(plats[p], CL_DEVICE_TYPE_GPU, 0, nullptr, &ndev) == CL_SUCCESS && ndev > 0) {
            std::vector<cl_device_id> devs(ndev);
            clGetDeviceIDs(plats[p], CL_DEVICE_TYPE_GPU, ndev, devs.data(), nullptr);
            for (cl_uint d = 0; d < ndev; ++d) {
                char exts[2048] = {0};
                clGetDeviceInfo(devs[d], CL_DEVICE_EXTENSIONS, sizeof(exts), exts, nullptr);
                if (strstr(exts, "cl_khr_fp64")) { chosen = devs[d]; impl_->o.platform = plats[p]; break; }
            }
            if (chosen) break;
        }
    }
    if (!chosen) return false;               // no fp64 GPU device
    impl_->o.device = chosen;

    cl_uint cu = 0, clk = 0;
    clGetDeviceInfo(chosen, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cu), &cu, nullptr);
    clGetDeviceInfo(chosen, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(clk), &clk, nullptr);
    cu_ = (int)cu; clock_ = (int)clk;
    char dname[256] = {0};
    clGetDeviceInfo(chosen, CL_DEVICE_NAME, sizeof(dname), dname, nullptr);
    dev_name_ = dname;

    impl_->o.ctx = clCreateContext(nullptr, 1, &chosen, nullptr, nullptr, &err);
    if (err != CL_SUCCESS) { std::fprintf(stderr, "[opencl_accel] context: %s\n", err_str(err)); return false; }
    impl_->o.q = clCreateCommandQueue(impl_->o.ctx, chosen, 0, &err);
    if (err != CL_SUCCESS) { std::fprintf(stderr, "[opencl_accel] queue: %s\n", err_str(err)); return false; }

    impl_->o.prog = clCreateProgramWithSource(impl_->o.ctx, 1, &KERNEL_SRC, nullptr, &err);
    if (err != CL_SUCCESS) { std::fprintf(stderr, "[opencl_accel] source: %s\n", err_str(err)); return false; }
    err = clBuildProgram(impl_->o.prog, 1, &chosen, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS) {
        char log[8192] = {0};
        clGetProgramBuildInfo(impl_->o.prog, chosen, CL_PROGRAM_BUILD_LOG, sizeof(log), log, nullptr);
        std::fprintf(stderr, "[opencl_accel] build: %s\n%s\n", err_str(err), log);
        return false;
    }
    impl_->o.k_gemm = clCreateKernel(impl_->o.prog, "k_gemm", &err);
    impl_->o.k_syrk = clCreateKernel(impl_->o.prog, "k_syrk_u", &err);
    impl_->o.k_rank2 = clCreateKernel(impl_->o.prog, "k_rank2", &err);
    if (err != CL_SUCCESS || !impl_->o.k_gemm || !impl_->o.k_syrk || !impl_->o.k_rank2) {
        std::fprintf(stderr, "[opencl_accel] kernel create: %s\n", err_str(err));
        return false;
    }
    impl_->built = true;
    ok_ = true;
    return true;
}

bool OpenCLAccel::gemm(int n, const double* A, const double* B, double* C) {
    if (!ok_) return false;
    cl_int err;
    size_t nb = (size_t)n * n;
    cl_mem dA = clCreateBuffer(impl_->o.ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, nb * sizeof(double), (void*)A, &err);
    cl_mem dB = clCreateBuffer(impl_->o.ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, nb * sizeof(double), (void*)B, &err);
    cl_mem dC = clCreateBuffer(impl_->o.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nb * sizeof(double), (void*)C, &err);
    if (err != CL_SUCCESS) return false;
    clSetKernelArg(impl_->o.k_gemm, 0, sizeof(cl_mem), &dA);
    clSetKernelArg(impl_->o.k_gemm, 1, sizeof(cl_mem), &dB);
    clSetKernelArg(impl_->o.k_gemm, 2, sizeof(cl_mem), &dC);
    clSetKernelArg(impl_->o.k_gemm, 3, sizeof(int), &n);
    size_t gws[2] = {(size_t)n, (size_t)n};
    err = clEnqueueNDRangeKernel(impl_->o.q, impl_->o.k_gemm, 2, nullptr, gws, nullptr, 0, nullptr, nullptr);
    if (err == CL_SUCCESS) err = clEnqueueReadBuffer(impl_->o.q, dC, CL_TRUE, 0, nb * sizeof(double), C, 0, nullptr, nullptr);
    clReleaseMemObject(dA); clReleaseMemObject(dB); clReleaseMemObject(dC);
    return err == CL_SUCCESS;
}

bool OpenCLAccel::syrk_upper(int n, int kk, double alpha, const double* A, double* C) {
    if (!ok_) return false;
    cl_int err;
    size_t nb = (size_t)n * n, ab = (size_t)n * kk;
    cl_mem dA = clCreateBuffer(impl_->o.ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, ab * sizeof(double), (void*)A, &err);
    cl_mem dC = clCreateBuffer(impl_->o.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nb * sizeof(double), (void*)C, &err);
    if (err != CL_SUCCESS) return false;
    clSetKernelArg(impl_->o.k_syrk, 0, sizeof(cl_mem), &dA);
    clSetKernelArg(impl_->o.k_syrk, 1, sizeof(cl_mem), &dC);
    clSetKernelArg(impl_->o.k_syrk, 2, sizeof(int), &n);
    clSetKernelArg(impl_->o.k_syrk, 3, sizeof(int), &kk);
    clSetKernelArg(impl_->o.k_syrk, 4, sizeof(double), &alpha);
    size_t gws[2] = {(size_t)n, (size_t)n};
    err = clEnqueueNDRangeKernel(impl_->o.q, impl_->o.k_syrk, 2, nullptr, gws, nullptr, 0, nullptr, nullptr);
    if (err == CL_SUCCESS) err = clEnqueueReadBuffer(impl_->o.q, dC, CL_TRUE, 0, nb * sizeof(double), C, 0, nullptr, nullptr);
    clReleaseMemObject(dA); clReleaseMemObject(dC);
    return err == CL_SUCCESS;
}

bool OpenCLAccel::rank2_update(int n, int k2, const double* c2, double a,
                               int k1, const double* c1, double b, double* C) {
    if (!ok_) return false;
    cl_int err;
    size_t nb = (size_t)n * n;
    size_t ab2 = (size_t)n * (size_t)(k2 > 0 ? k2 : 1), ab1 = (size_t)n * (size_t)(k1 > 0 ? k1 : 1);
    cl_mem dC2 = clCreateBuffer(impl_->o.ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, ab2 * sizeof(double), (void*)c2, &err);
    cl_mem dC1 = clCreateBuffer(impl_->o.ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, ab1 * sizeof(double), (void*)c1, &err);
    cl_mem dC = clCreateBuffer(impl_->o.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, nb * sizeof(double), (void*)C, &err);
    if (err != CL_SUCCESS) return false;
    clSetKernelArg(impl_->o.k_rank2, 0, sizeof(cl_mem), &dC2);
    clSetKernelArg(impl_->o.k_rank2, 1, sizeof(int), &k2);
    clSetKernelArg(impl_->o.k_rank2, 2, sizeof(double), &a);
    clSetKernelArg(impl_->o.k_rank2, 3, sizeof(cl_mem), &dC1);
    clSetKernelArg(impl_->o.k_rank2, 4, sizeof(int), &k1);
    clSetKernelArg(impl_->o.k_rank2, 5, sizeof(double), &b);
    clSetKernelArg(impl_->o.k_rank2, 6, sizeof(cl_mem), &dC);
    clSetKernelArg(impl_->o.k_rank2, 7, sizeof(int), &n);
    size_t gws[2] = {(size_t)n, (size_t)n};
    err = clEnqueueNDRangeKernel(impl_->o.q, impl_->o.k_rank2, 2, nullptr, gws, nullptr, 0, nullptr, nullptr);
    if (err == CL_SUCCESS) err = clEnqueueReadBuffer(impl_->o.q, dC, CL_TRUE, 0, nb * sizeof(double), C, 0, nullptr, nullptr);
    clReleaseMemObject(dC2); clReleaseMemObject(dC1); clReleaseMemObject(dC);
    return err == CL_SUCCESS;
}

// ------------------------------------------------------------ self test
namespace {
void cpu_gemm(int n, const double* A, const double* B, double* C) {
    for (int i = 0; i < n; ++i) {
        const double* Ai = A + (size_t)i * n;
        double* Ci = C + (size_t)i * n;
        for (int j = 0; j < n; ++j) {
            double s = 0.0;
            const double* Bj = B + (size_t)j;
            for (int k = 0; k < n; ++k) s += Ai[k] * Bj[(size_t)k * n];
            Ci[j] += s;
        }
    }
}
void cpu_syrk(int n, int kk, double alpha, const double* A, double* C) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) {
            double s = 0.0;
            for (int k = 0; k < kk; ++k) s += A[(size_t)i * kk + k] * A[(size_t)j * kk + k];
            C[(size_t)i * n + j] = alpha * s;
        }
}
void cpu_rank2(int n, int k2, const double* c2, double a, int k1, const double* c1, double b, double* C) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) {
            double s = 0.0;
            for (int k = 0; k < k2; ++k) s += c2[(size_t)i * k2 + k] * c2[(size_t)j * k2 + k];
            double s1 = 0.0;
            for (int k = 0; k < k1; ++k) s1 += c1[(size_t)i * k1 + k] * c1[(size_t)j * k1 + k];
            C[(size_t)i * n + j] = a * s + b * s1;
        }
}
double maxrel(const double* a, const double* b, size_t n, double& refmax) {
    double m = 0.0; refmax = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double d = std::fabs(a[i] - b[i]);
        if (d > m) m = d;
        if (std::fabs(b[i]) > refmax) refmax = std::fabs(b[i]);
    }
    return m;
}
} // namespace

bool OpenCLAccel::selftest() {
    if (!ok_) return false;
    const int n = 96, kk = 40;
    std::vector<double> A((size_t)n * kk), B((size_t)n * n), C1((size_t)n * n), C2((size_t)n * n);
    for (auto& x : A) x = (double)(rand() % 2000) / 13.0 - 70.0;
    for (auto& x : B) x = (double)(rand() % 2000) / 17.0 - 60.0;
    for (auto& x : C1) x = (double)(rand() % 2000) / 19.0 - 50.0;
    for (auto& x : C2) x = (double)(rand() % 2000) / 23.0 - 40.0;

    std::vector<double> ref((size_t)n * n, 0.0), got((size_t)n * n, 0.0);
    cpu_gemm(n, A.data(), B.data(), ref.data());
    gemm(n, A.data(), B.data(), got.data());
    double rm = 0.0, md = maxrel(got.data(), ref.data(), (size_t)n * n, rm);
    bool ok = md / (rm + 1e-300) < 1e-12;
    std::printf("[opencl_accel] gemm  rel-err %.2e %s\n", md / (rm + 1e-300), ok ? "OK" : "FAIL");

    std::fill(ref.begin(), ref.end(), 0.0); std::fill(got.begin(), got.end(), 0.0);
    cpu_syrk(n, kk, 2.0, A.data(), ref.data());
    syrk_upper(n, kk, 2.0, A.data(), got.data());
    rm = 0.0; md = maxrel(got.data(), ref.data(), (size_t)n * n, rm);
    bool ok2 = md / (rm + 1e-300) < 1e-12;
    std::printf("[opencl_accel] syrk  rel-err %.2e %s\n", md / (rm + 1e-300), ok2 ? "OK" : "FAIL");

    std::fill(ref.begin(), ref.end(), 0.0); std::fill(got.begin(), got.end(), 0.0);
    cpu_rank2(n, kk, A.data(), 2.0, 24, B.data(), 0.5, ref.data());
    rank2_update(n, kk, A.data(), 2.0, 24, B.data(), 0.5, got.data());
    rm = 0.0; md = maxrel(got.data(), ref.data(), (size_t)n * n, rm);
    bool ok3 = md / (rm + 1e-300) < 1e-12;
    std::printf("[opencl_accel] rank2 rel-err %.2e %s\n", md / (rm + 1e-300), ok3 ? "OK" : "FAIL");

    return ok && ok2 && ok3;
}

#else // !MOPAC_USE_OPENCL -------------------------------------------------

// Stub build: module present but inert.  Callers always see available()==false
// and use the CPU path, so the bitwise reference behaviour is untouched.

OpenCLAccel::OpenCLAccel() : impl_(nullptr), ok_(false), cu_(0), clock_(0) {}
OpenCLAccel::~OpenCLAccel() {}
bool OpenCLAccel::init(const std::string&) { return false; }
const std::string& OpenCLAccel::device_name() const { static const std::string s("(OpenCL disabled)"); return s; }
bool OpenCLAccel::gemm(int, const double*, const double*, double*) { return false; }
bool OpenCLAccel::syrk_upper(int, int, double, const double*, double*) { return false; }
bool OpenCLAccel::rank2_update(int, int, const double*, double, int, const double*, double, double*) { return false; }
bool OpenCLAccel::selftest() { return false; }

#endif // MOPAC_USE_OPENCL
