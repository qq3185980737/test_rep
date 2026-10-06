// parameters_for_mndo_C.h — C++ translation of MOPAC 2016 "parameters_for_mndo_C.F90".
#pragma once

#include <vector>

namespace parameters_for_mndo_C {
extern std::vector<double> alpm;  // (107)
extern std::vector<double> betadm;  // (107)
extern std::vector<double> betapm;  // (107)
extern std::vector<double> betasm;  // (107)
extern std::vector<double> f0sdm;  // (107)
extern std::vector<double> g2sdm;  // (107)
extern std::vector<double> gp2m;  // (107)
extern std::vector<double> gppm;  // (107)
extern std::vector<double> gspm;  // (107)
extern std::vector<double> gssm;  // (107)
extern std::vector<std::vector<double>> guesm1;  // (107,4)
extern std::vector<std::vector<double>> guesm2;  // (107,4)
extern std::vector<std::vector<double>> guesm3;  // (107,4)
extern std::vector<double> hspm;  // (107)
extern std::vector<double> pocm;  // (107)
extern std::vector<double> polvolm;  // (107)
extern std::vector<double> polvom;  // (107)
extern std::vector<double> uddm;  // (107)
extern std::vector<double> uppm;  // (107)
extern std::vector<double> ussm;  // (107)
extern std::vector<double> zdm;  // (107)
extern std::vector<double> zdnm;  // (107)
extern std::vector<double> zpm;  // (107)
extern std::vector<double> zpnm;  // (107)
extern std::vector<double> zsm;  // (107)
extern std::vector<double> zsnm;  // (107)
}  // namespace parameters_for_mndo_C
