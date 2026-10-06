// parameters_for_AM1_Sparkles_C.h — C++ translation of MOPAC 2016 "parameters_for_AM1_Sparkles_C.F90".
#pragma once

#include <vector>

namespace parameters_for_AM1_Sparkles_C {
extern std::vector<double> alpam1sp;  // (107)
extern std::vector<double> gssam1sp;  // (107)
extern std::vector<std::vector<double>> guesam1sp1;  // (107,4)
extern std::vector<std::vector<double>> guesam1sp2;  // (107,4)
extern std::vector<std::vector<double>> guesam1sp3;  // (107,4)
}  // namespace parameters_for_AM1_Sparkles_C
