// parameters_for_PM3_Sparkles_C.h — C++ translation of MOPAC 2016 "parameters_for_PM3_Sparkles_C.F90".
#pragma once

#include <vector>

namespace parameters_for_PM3_Sparkles_C {
extern std::vector<double> alpPM3sp;  // (107)
extern std::vector<double> gssPM3sp;  // (107)
extern std::vector<std::vector<double>> guesPM3sp1;  // (107,4)
extern std::vector<std::vector<double>> guesPM3sp2;  // (107,4)
extern std::vector<std::vector<double>> guesPM3sp3;  // (107,4)
}  // namespace parameters_for_PM3_Sparkles_C
