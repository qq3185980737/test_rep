// H_bond_correction_bits.h — C++ translation of MOPAC 2016.
#pragma once
#include <vector>

double truncation(double R, double limit, double spread);

// Distance (Angstrom) between atoms a,b; periodic-image minimum if id!=0.
double distance_hb(int a, int b);
// Angle (radians) a-b-c.
double angle_hb(int a, int b, int c);
// Torsion (radians) i-j-k-l.
double torsion_hb(int i, int j, int k, int l);
// Sum of covalent radii for atoms x,y.
double bonding(int x, int y, const double* covrad);

// connected: true if atoms a,b within sqrt(criterion) Angstrom.
bool connected_hb(int atom_i, int atom_j, double criterion);
// find all H atoms bonded to N/O (or S).
void find_XH_bonds(std::vector<int>& acc, int& nacc,
                   std::vector<int>& h_b, int& nhb);
// Find all H-bond triples (donor-O/N, H, acceptor-O/N).
void find_H__Y_bonds(const std::vector<int>& acc_a, int nacc_a,
                     const std::vector<int>& acc_b, int nacc_b,
                     const std::vector<int>& bonding_a_h, int nb_a_h,
                     std::vector<int>& hblist1, std::vector<int>& hblist2,
                     std::vector<int>& hblist3, int max_h_bonds, int& nrpairs);
// Top-level: allocate and find all H-bonds.
void all_h_bonds(std::vector<int>& hblist1, std::vector<int>& hblist2,
                 std::vector<int>& hblist3, int max_h_bonds, int& nrpairs);

// Print H-bond diagnostics (populates H_txt / H_energy / P_Hbonds).
void prt_hbonds(int D, int H, int A, double energy);

