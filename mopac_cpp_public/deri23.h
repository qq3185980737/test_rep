// deri23.h — C++ translation of MOPAC 2016 "deri23.F90".
#pragma once
#include <vector>

// Unpack CI-active MO derivatives into cmo; eigenvalue relaxation into emo.
// f/fd/fci are column-major supervectors over the full basis (1-based
// padding); ivar selects the current geometric-variable column so the caller
// passes the same container for every column (as F90 f(1,ivar)).
void deri23(const std::vector<double>& f,
            const std::vector<double>& fd,
            const std::vector<double>& e,
            const std::vector<double>& fci,
            std::vector<std::vector<double>>& cmo,
            std::vector<double>& emo, int minear, int ninear, int ivar);
