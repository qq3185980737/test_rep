// datin.h — C++ translation of MOPAC 2016 "datin.F90".
#pragma once
#include <string>
#include <vector>

// Parse the EXTERNAL=/PARAMS= file list out of the keyword string.
std::vector<std::string> split_param_files(const std::string& keywrd);

// Element symbol by atomic number (1..107), matching MOPAC elemnt table.
std::string elemnt_sym(int z);

// Full parameter-file read-in (file I/O; body reads files and calls update()).
void datin(int iw);
