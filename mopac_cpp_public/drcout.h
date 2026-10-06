// drcout.h — C++ translation of MOPAC 2016 "drcout.F90".
#pragma once
#include <string>
#include <vector>

void drcout(const std::vector<std::vector<double>>& xyz3,
            const std::vector<std::vector<double>>& geo3,
            const std::vector<std::vector<double>>& vel3,
            int nvar, double time,
            const std::vector<double>& escf3, const std::vector<double>& ekin3,
            const std::vector<double>& etot3, const std::vector<double>& xtot3,
            int iloop, const std::vector<double>& charge, double fract,
            const std::string& text1, const std::string& text2,
            int ii, int& jloop);
