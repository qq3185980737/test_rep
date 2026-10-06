// to_screen.cpp — screen/log messages (F90: write(0,...)); keep them off the
// .out file, which is stdout.
#include "to_screen.h"
#include "chanel_C.h"
#include <cstdio>
void to_screen(const std::string& text) {
    if (chanel_C::iw0 > -1) { std::fprintf(stderr, "%s\n", text.c_str()); }
}
