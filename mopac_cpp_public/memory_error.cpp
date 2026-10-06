// memory_error.cpp — C++ translation of the runtime support routine.
#include "memory_error.h"

#include <cstdio>
#include <cstdlib>

void memory_error(const char* where) {
    std::printf(" MOPAC  Memory allocation failure in %s\n", where ? where : "?");
    std::fflush(stdout);
    std::exit(1);
}
