// getgeo.h — C++ translation of MOPAC 2016 "getgeo.F90".
#pragma once
#include <vector>
void getgeo(int iread, std::vector<int>& labels,
            std::vector<std::vector<double>>& geo,
            std::vector<std::vector<double>>& xyz,
            std::vector<std::vector<int>>& lopt,
            std::vector<int>& na, std::vector<int>& nb, std::vector<int>& nc,
            bool int_);
