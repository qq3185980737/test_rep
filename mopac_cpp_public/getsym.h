// getsym.h — C++ translation of MOPAC 2016 "getsym.F90".
#pragma once
#include <vector>
void getsym(std::vector<int>& locpar, std::vector<int>& idepfn,
            std::vector<int>& locdep, std::vector<double>& depmul);
