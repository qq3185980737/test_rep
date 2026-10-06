// digit.h — C++ translation of MOPAC 2016 "digit.F90".
#pragma once
#include <string>

// Convert a clean numeric substring to double, 1-based istart.
double digit(const std::string& s, int istart);
