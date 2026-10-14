// test_vast.cpp
#include <cstdio>
#include "vastkind.h"
int main() {
    static_assert(sizeof(vast_kind)==8, "");
    std::printf("vastkind PASS\n"); return 0;
}
