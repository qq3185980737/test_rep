// axis.h — C++ translation of MOPAC 2016 "axis.F90".
#pragma once

// axis: compute the three principal moments of inertia and molecular weight.
// Rotational constants (cm^-1) are returned in a, b, c; eigenvectors evec(1..3,1..3).
void axis(double& a, double& b, double& c, double evec[4][4]);
