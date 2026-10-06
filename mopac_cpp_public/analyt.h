// analyt.h — C++ translation of MOPAC 2016 "analyt.F90".
#pragma once
#include <vector>

// analyt: analytical derivatives between two atoms i=2, j=1.
//   psum/palpha/pbeta: packed 2D density / Pi arrays (size mpack+1).
//   coord[dim][atom] (atom 1..2). nat[1..2] = atomic numbers.
//   eng[1..3] output derivative w.r.t. each Cartesian component.
void analyt(const std::vector<double>& psum, const std::vector<double>& palpha,
            const std::vector<double>& pbeta,
            const std::vector<std::vector<double>>& coord,
            const int nat[3], int jja, int jjd, int iia, int iid,
            double eng[4]);

// Internal (also called from Fortran as separate routines).
void delmol(const std::vector<std::vector<double>>& coord, int i, int j,
            int ni, int nj, int ia, int id, int ja, int jd, int ix,
            double rij, double tomb, int& isp);
void delri(double dg[23], int ni, int nj, double rr, double del1);
void rotat(const std::vector<std::vector<double>>& coord, int i, int j, int ix,
           int idx, double rij, double del1);
void ders(int m, int n, double rr, double del1, double del2, double del3,
          int is, int iol);
