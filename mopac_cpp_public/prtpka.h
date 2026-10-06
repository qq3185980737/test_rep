// prtpka.h — C++ translation of MOPAC 2016 "prtpka.F90".
#pragma once

// Calculates and prints the pKa values for an organic molecule with an
// ionizable hydrogen attached to an oxygen atom.
void prtpka(int* ipKa_sorted, double* pKa_sorted, int* ipKa_unsorted,
            double* pKa_unsorted, int& no);

// pKa parameter set + NDDO parameter-table overrides (H/C/N/O/F/Cl/Br/I).
void Parameters_for_PKA(double& c1, double& c2, double& c3);
