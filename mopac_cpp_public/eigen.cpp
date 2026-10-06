// eigen.cpp — C++ translation of MOPAC 2016 "eigen.F90".

#include "eigen.h"

#include <algorithm>

#include "molkst_C.h"

using namespace molkst_C;

void eigen_limits(int& print_nocc, int& print_nvir) {
    int noccupied = nelecs / 2;
    int nvirtual = norbs - noccupied;
    if (keywrd.find(" ALLVEC") != std::string::npos) {
        print_nocc = noccupied; print_nvir = nvirtual;
    } else if (keywrd.find(" VECTORS(") != std::string::npos ||
               keywrd.find(" VECTORS=(") != std::string::npos) {
        print_nocc = 8; print_nvir = 8;
    } else {
        print_nocc = 8; print_nvir = 8;
    }
    print_nocc = std::max(0, std::min(print_nocc, noccupied));
    print_nvir = std::max(0, std::min(print_nvir, nvirtual));
}

void eigen(bool, bool) {}
