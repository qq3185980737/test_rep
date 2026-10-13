// _test_stubs.cpp -- link-time stubs for interface tests only.
// Provides pinout() so tidy.cpp links without the full driver deps.
#include <cstdio>
void pinout(int) {
    std::fprintf(stderr, "PINOUT CALLED\n");
    std::fflush(stderr);
}
