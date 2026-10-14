// test_opencl_accel.cpp - self-test harness for the optional OpenCL module.
// Build (enabled):  cl /O2 /DMOPAC_USE_OPENCL /I <SDK>\include test_opencl_accel.cpp ..\opencl_accel.cpp /link <SDK>\lib\x86_64\opencl.lib
// Build (stub):    cl /O2 test_opencl_accel.cpp ..\opencl_accel.cpp
#include "../opencl_accel.h"

#include <cstdio>
#include <string>

int main(int argc, char** argv) {
    OpenCLAccel acc;
    std::string plat = argc > 1 ? argv[1] : "";
    if (!acc.init(plat)) {
        std::printf("OpenCL unavailable (platform filter '%s')\n", plat.c_str());
        return 2;
    }
    std::printf("device : %s\n", acc.device_name().c_str());
    std::printf("compute units : %d, clock : %d MHz\n", acc.compute_units(), acc.max_clock_mhz());
    bool ok = acc.selftest();
    std::printf(ok ? "SELFTEST PASS\n" : "SELFTEST FAIL\n");
    return ok ? 0 : 1;
}
