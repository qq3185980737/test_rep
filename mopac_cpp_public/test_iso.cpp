// test_iso.cpp
#include <cstdio>
#include "mod_iso_c_binding.h"
int main() {
    static_assert(iso_c_binding::c_double==8 && sizeof(iso_c_binding::c_ptr)==sizeof(void*), "");
    std::printf("iso_c_binding PASS\n");
    return 0;
}
