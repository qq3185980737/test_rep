// genun.h — C++ translation of MOPAC 2016 "genun.F90".
#pragma once
#include <vector>
void genun(std::vector<std::vector<double>>& u, int& n);
bool collid(double rw, const std::vector<double>& cw,
            const std::vector<std::vector<double>>& cnbr,
            const std::vector<double>& rnbr, int nnbr, int ishape);
