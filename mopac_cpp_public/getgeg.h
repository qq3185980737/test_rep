// getgeg.h — C++ translation of MOPAC 2016 "getgeg.F90".
#pragma once
#include <vector>
void getgeg(int iread, std::vector<int>& labels, std::vector<std::vector<double>>& geo,
            std::vector<std::vector<int>>& lopt, std::vector<int>& na,
            std::vector<int>& nb, std::vector<int>& nc);
