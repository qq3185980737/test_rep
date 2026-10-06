// geo_diff.h — C++ translation of readmo.F90 "geo_diff" (MOPAC 2016).
#pragma once

// geo_diff calculates the total distance between geometries "geo" and "geoa":
//   sum = addition of all differences in position (sum of motions of atoms).
//   rms = sum of squares of differences in positions.
// If prt is true, lists the atoms that move a lot.
void geo_diff(double& sum, double& rms, bool prt);
