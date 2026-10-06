// eigen.h — C++ translation of MOPAC 2016 "eigen.F90".
#pragma once

// Parse VECTORS keyword -> print ranges.
void eigen_limits(int& print_nocc, int& print_nvir);
void eigen(bool write_gpt, bool write_out);
