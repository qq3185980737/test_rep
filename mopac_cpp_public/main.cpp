// main.cpp — C++ translation of "MOPAC for Windows.F90" (program MOPAC_win).
// The original uses Intel QuickWin (IFQWIN/DFLIB) for a windowed menu GUI.
// That platform glue is dropped; this provides a console entry point that
// performs the same core initialization and runs MOPAC.

#include <cstdio>
#include <cstdlib>
#include <string>

#include "chanel_C.h"
#include "molkst_C.h"

// Forward declarations of the computational entry points.
void GetDateStamp(std::string& line, std::string& verson);
void run_mopac();

using namespace chanel_C;
using namespace molkst_C;

int main(int argc, char** argv) {
    gui = false;
    iw0 = 0;
    ijulian = 123;
    site_no = 88;
    run = 1;                    // argv[1] = data-set filename
    program_name = "Standalone MOPAC ";

    if (argc <= 1) {
        std::printf("\n");
        std::printf("  MOPAC: supply a data-set filename as the first argument.\n");
        return EXIT_SUCCESS;
    }

    GetDateStamp(line, verson);
    run_mopac();
    return EXIT_SUCCESS;
}
