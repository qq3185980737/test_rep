// mod_gpu_info.h — C++ translation.
#pragma once
#include <cstddef>
extern "C" {
void getGPUInfo(int* hasGpu, int* hasDouble, int* nDevices, char* name,
                int* name_size, size_t* totalMem, int* clockRate, int* major,
                int* minor);
void setDevice(int idevice, int* stat);
}
