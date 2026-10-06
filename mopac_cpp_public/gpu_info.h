// gpu_info.h — C++ binding for the external GPU info / device-selection C API.
// Mirrors mod_gpu_info.F90 (getGPUInfo / setDevice).
#pragma once
#include <cstddef>
#include <cstdint>

extern "C" {
// Query GPU availability. Arrays are sized up to 6 devices.
void getGPUInfo(bool* hasGpu, bool* hasDouble, int* nDevices,
                char* name, int* name_size, std::size_t* totalMem,
                int* clockRate, int* major, int* minor);

// Select the active CUDA device.
void setGPU(int idevice, bool* stat);
}
