// diagg.h — C++ translation of diagg.F90 / diagg1.F90 / diagg2.F90 /
// density_for_MOZYME.F90 / epseta.F90 (MOPAC 2016).
#pragma once
#include <vector>

// diagg: most important operation of MOZYME - the matrix elements involving
// occupied and virtual LMO's are annihilated by a single pass Euler rotation.
//   fao:   FOCK matrix over atomic orbitals (packed, mpack)
//   partp: partial density matrix (mpack)
void diagg(const std::vector<double>& fao, int nocc, int nvir, int idiagg,
           std::vector<double>& partp, int indi);

// diagg1: build the significant occupied-virtual matrix elements (stored in
// fmo/ifmo, count nij).  eigv is the eigs(nocc+1:) slice (full eigs passed,
// offset nocc).
void diagg1(const std::vector<double>& fao, int nocc, int nvir,
            std::vector<double>& eigv, std::vector<double>& ws,
            std::vector<bool>& latoms, std::vector<std::vector<int>>& ifmo,
            std::vector<double>& fmo, int fmo_dim, int& nij, int idiagg,
            std::vector<double>& avir, std::vector<double>& aocc,
            std::vector<double>& aov);

// diagg2: annihilate the significant elements by a two-by-two Jacobi rotation.
void diagg2(int nocc, int nvir, std::vector<double>& eigv,
            std::vector<int>& iused, std::vector<bool>& latoms, int nij,
            int idiagg, std::vector<double>& storei, std::vector<double>& storej);

// density_for_MOZYME: build (partial/full) density matrix from LMO coefficients.
// mode: 0 = full density (zero P first), -1 = remove old-LMO density, else add.
void density_for_MOZYME(std::vector<double>& p, int mode, int nclose_loc,
                        const std::vector<double>& partp);

// epseta: machine constants; eps = smallest number with 1+eps != 1,
// eta = smallest representable magnitude (guarded to 1e-39 / 1e-17).
void epseta(double& eps, double& eta);