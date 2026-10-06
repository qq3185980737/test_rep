// mod_vars_cuda.h — C++ translation.
#pragma once
namespace mod_vars_cuda {
    constexpr int nthreads_gpu = 256;
    constexpr int nblocks_gpu = 256;
    extern bool lgpu;
    constexpr int prec = 8;
    extern int ngpus, gpu_id;
}
