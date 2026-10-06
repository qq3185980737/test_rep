// convert_storage.h — C++ translation of MOPAC 2016 "convert_storage.F90".
#pragma once
#include <vector>

// MOZYME packed -> lower-triangle packed matrix.
void convert_mat_packed_to_triangle(const std::vector<double>& matrix_packed,
                                    std::vector<double>& matrix_triangle);

// MOZYME packed LMOs -> conventional square eigenvector matrix.
void convert_lmo_packed_to_square(std::vector<std::vector<double>>& c_square);
