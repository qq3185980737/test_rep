// pmep.h — C++ translation of "pmep.F90" (12 program units).
#pragma once
#include <vector>

// pmep: driving routine for PMEP — Parametric Molecular Electrostatic
// Potential and MEP charges (AM1 only).
void pmep();

// pmepco: MEP at an arbitrary position w(3) (core subroutine).
void pmepco(const std::vector<double>& pp, const std::vector<double>& ria,
            const std::vector<double>& w, double& ui,
            const std::vector<std::vector<double>>& co, int nonzo, int id);

// surfa: Connolly molecular surface sample points (four vdW layers).
void surfa(const std::vector<std::vector<double>>& co,
           std::vector<std::vector<double>>& potpt, int& nmep);

// collis: collision check of probe with neighboring atoms.
bool collis(const std::vector<double>& cw, double rw,
            const std::vector<std::vector<double>>& cnbr,
            const std::vector<double>& rnbr, int nnbr, int ishape);

// drepp2: two-electron repulsion integrals for PMEP (core[4][3] mirrors
// rotate_C::ccore, 1-based).
void drepp2(int ni, double rij, std::vector<double>& ri, double core[4][3]);

// drotat: rotate integrals from local to molecular coordinates.
void drotat(int ni, const std::vector<double>& xi, const std::vector<double>& xj,
            std::vector<double>& e1b, double& enuc, double rij);

// genvec: generate unit vectors over a sphere (Connolly surface).
void genvec(std::vector<std::vector<double>>& u, int& n);

// grids: Williams surface (grid) sampling.
void grids(const std::vector<std::vector<double>>& co,
           std::vector<std::vector<double>>& potpt, int& nmep);

// mepchg: MEP charges by fitting the quantum potential to Coulomb.
void mepchg(const std::vector<std::vector<double>>& co,
            const std::vector<std::vector<double>>& potpt, int nmep);

// mepmap: MEP minima and contour data in a defined plane.
void mepmap(std::vector<std::vector<double>>& c1);

// meprot: rotate the molecule to choose the MEP plane.
void meprot(const std::vector<std::vector<double>>& c, std::vector<double>& r0,
            std::vector<std::vector<double>>& t, double& xm, double& ym,
            int& icase, double step, std::vector<std::vector<double>>& c1,
            double& z0, int& iback);

// packp: re-write the density matrix (block-triangular per atom).
void packp(const std::vector<double>& p, std::vector<double>& pp, int& mn);
