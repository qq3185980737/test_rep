// parameters_for_RM1_C.h — C++ translation of MOPAC 2016 "parameters_for_RM1_C.F90".
#pragma once

#include <vector>

namespace parameters_for_RM1_C {
extern std::vector<double> alpRM1;  // (107)
extern std::vector<double> betapRM1;  // (107)
extern std::vector<double> betasRM1;  // (107)
extern std::vector<double> gp2RM1;  // (107)
extern std::vector<double> gppRM1;  // (107)
extern std::vector<double> gspRM1;  // (107)
extern std::vector<double> gssRM1;  // (107)
extern std::vector<std::vector<double>> guess1RM1;  // (107,4)
extern std::vector<std::vector<double>> guess2RM1;  // (107,4)
extern std::vector<std::vector<double>> guess3RM1;  // (107,4)
extern std::vector<double> hspRM1;  // (107)
extern std::vector<double> uppRM1;  // (107)
extern std::vector<double> ussRM1;  // (107)
extern std::vector<double> zpRM1;  // (107)
extern std::vector<double> zsRM1;  // (107)
}  // namespace parameters_for_RM1_C
