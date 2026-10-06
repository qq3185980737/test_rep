// formxy.h — C++ translation of MOPAC 2016 "formxy.F90".
#pragma once
#include <vector>
void formxy(const std::vector<double>& w, int& kr,
            std::vector<double>& wca, std::vector<double>& wcb,
            const std::vector<double>& ca, const std::vector<double>& cb,
            int na, int nb);
