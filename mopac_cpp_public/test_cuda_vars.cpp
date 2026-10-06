// test_cuda_vars.cpp
#include <cstdio>
#include "mod_vars_cuda.h"
int main() {
    mod_vars_cuda::lgpu=true; mod_vars_cuda::ngpus=1;
    std::printf("threads=%d ngpus=%d gpu=%d PASS\n",mod_vars_cuda::nthreads_gpu,mod_vars_cuda::ngpus,mod_vars_cuda::gpu_id);
    return 0;
}
