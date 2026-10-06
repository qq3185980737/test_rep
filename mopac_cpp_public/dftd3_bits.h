// dftd3_bits.h — C++ translation of MOPAC 2016 "dftd3_bits.F90" + "gdisp.F90".
#pragma once
#include <string>
#include <vector>

using c6ab_t = std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>;

// Analytical derivative dC6/dr (pure blend of C6ab reference values).
void anagrdc6(int max_elem, int maxc, int n,
              const std::vector<double>& cn,
              const std::vector<std::vector<double>>& dcn2,
              const std::vector<std::vector<std::vector<double>>>& dcn3,
              const std::vector<int>& nat,
              const std::vector<int>& mxc,
              int iat, int jat, int kat,
              const c6ab_t& c6ab,
              std::vector<double>& anag);

// limit(iat,jat,iadr,jadr) — coordination-number range (pure).
void d3limit(int iat, int jat, int& iadr, int& jadr);

// lin(i1,i2) — packed symmetric pair index  (pure).
int lin(int i1, int i2);

// ESYM(i) — two-letter element symbol.
std::string esym(int i);

// Hydrogen-bond two-body energy (pure).
double eabh(int n, int A, int B, int H,
            const std::vector<std::vector<double>>& xyz,
            double shortcut, double cab);

// hbpar(elem): element -> hydrogen-bond type (1..6: N,O,F,P,S,Cl).
int hbpar(int elem);

// Interpolated C6(iat,jat) from reference table, exponential weighting.
void getc6(int maxc, int max_elem, const c6ab_t& c6ab, const std::vector<int>& mxc,
           int iat, int jat, double nci, double ncj, double& c6);

// Grimme D3 dispersion energy (two-body C6/C8 with damping).
void edisp(int max_elem, int maxc, int n, const std::vector<std::vector<double>>& xyz,
           const std::vector<int>& nat, const c6ab_t& c6ab, const std::vector<int>& mxc,
           const std::vector<double>& r2r4, const std::vector<std::vector<double>>& r0ab,
           const std::vector<double>& rcov, double rs6, double rs8, double alp6, double alp8,
           double& e6, double& e8);

// Hydrogen-bond correction (Korth-type, A-H-B).
void hbsimple(int n, const std::vector<int>& at, std::vector<std::vector<double>>& xyz,
              double hbscale, double& energy, bool l_grad, std::vector<std::vector<double>>& g);

// setr0ab — fill r(i,j) cut-off radii table from packed data (Bohr->angstrom /autoang).
void setr0ab(int max_elem, double autoang, std::vector<std::vector<double>>& r);

// Coordination numbers (inverse damping counting function).
void ncoord(int natoms, const std::vector<double>& rcov, const std::vector<int>& nat,
            const std::vector<std::vector<double>>& xyz, std::vector<double>& cn);

// dC6/dCN interpolation (used by gdisp).
void get_dC6_dCNij(int maxc, int max_elem, const c6ab_t& c6ab, int mxci, int mxcj,
                   double cni, double cnj, int izi, int izj,
                   double& c6check, double& dc6i, double& dc6j);

// Grimme D3 dispersion gradient (gdisp.F90).
void gdisp(const std::vector<std::vector<double>>& xyz,
           const std::vector<std::vector<double>>& r0ab,
           double rs6, double alp6, const c6ab_t& c6ab, double s6,
           const std::vector<int>& mxc, const std::vector<double>& rcov,
           std::vector<std::vector<double>>& dxyz_temp);
