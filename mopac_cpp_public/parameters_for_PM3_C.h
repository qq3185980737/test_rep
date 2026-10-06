// parameters_for_PM3_C.h — C++ translation of MOPAC 2016 "parameters_for_PM3_C.F90".
#pragma once

#include <vector>

namespace parameters_for_PM3_C {
extern std::vector<double> alppm3;  // (107)
extern std::vector<double> betadp;  // (107)
extern std::vector<double> betapp;  // (107)
extern std::vector<double> betasp;  // (107)
extern std::vector<double> f0sdpm3;  // (107)
extern std::vector<double> g2sdpm3;  // (107)
extern std::vector<double> gp2pm3;  // (107)
extern std::vector<double> gpppm3;  // (107)
extern std::vector<double> gsppm3;  // (107)
extern std::vector<double> gsspm3;  // (107)
extern std::vector<std::vector<double>> guesp1;  // (107,4)
extern std::vector<std::vector<double>> guesp2;  // (107,4)
extern std::vector<std::vector<double>> guesp3;  // (107,4)
extern std::vector<double> hsppm3;  // (107)
extern std::vector<double> polvolpm3;  // (107)
extern std::vector<double> upppm3;  // (107)
extern std::vector<double> usspm3;  // (107)
extern std::vector<double> zdnpm3;  // (107)
extern std::vector<double> zpnpm3;  // (107)
extern std::vector<double> zppm3;  // (107)
extern std::vector<double> zsnpm3;  // (107)
extern std::vector<double> zspm3;  // (107)
}  // namespace parameters_for_PM3_C
