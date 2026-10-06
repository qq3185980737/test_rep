// common_arrays_C.h — C++ translation of the Fortran module "Common_arrays_C"
// (MOPAC 2016).
//
// Only the members referenced by already-translated sources are declared so
// far; the header will be extended as further .F90 files are ported.
//
// Indexing convention: these Fortran allocatable arrays keep their Fortran
// (1-based) indexing in C++ — element 0 is padding, and Fortran accesses
// such as eigs(k) / f(k) are written verbatim as eigs[k] / f[k].
#pragma once

#include <string>
#include <vector>

namespace common_arrays_C {

// f    — Fock matrix (alpha Fock if UHF), packed lower triangle.
// eigs — M.O. eigenvalues (alpha M.O. eigenvalues if UHF).
extern std::vector<double> f;
extern std::vector<double> eigs; extern std::vector<double> eigb;

// Globals introduced by add_hydrogen_atoms.F90.
extern std::vector<int> nat;          // atomic numbers of real atoms (1-based)
extern std::vector<int> nbonds;       // number of bonds on each atom
extern std::vector<int> labels;
extern std::vector<int> na;
extern std::vector<int> nb, nc;
extern std::vector<int> nfirst;
extern std::vector<double> uspd;   // one-electron integrals / ionization potentials (hcore.F90)
extern std::vector<int> nlast;
extern std::vector<int> breaks;       // breaks(400)
extern std::vector<std::vector<int>> ibonds;   // (15, maxatoms)
extern std::vector<std::vector<int>> lopt;
extern std::vector<std::vector<int>> loc;      // (2, 3*maxatoms)
extern std::vector<double> atmass;
extern std::vector<double> xparam;
extern std::vector<std::vector<double>> coord;     // (3, maxatoms)
extern std::vector<std::vector<double>> geo;       // (3, maxatoms)
extern std::vector<std::vector<double>> coorda;     // (3, maxatoms)
extern std::vector<std::vector<double>> geoa;        // (3, maxatoms)
extern std::vector<std::vector<double>> fcint;       // (4, maxatoms) intfc finite-difference forces
extern std::vector<std::vector<double>> break_coords; // (3, 400)
extern std::vector<char> l_atom; extern std::vector<double> hesinv; extern std::vector<double> profil;
extern int time_start[8];  // date_and_time VALUES (year..ms)
extern std::vector<double> ch;
extern std::vector<std::string> txtatm;     // *26 per atom
extern std::vector<std::string> txtatm1;    // *26 per atom
extern std::vector<std::string> all_comments; // *81 per comment line
extern std::vector<std::vector<double>> tvec; // (3,3) periodic lattice vectors
extern std::string chains;                  // chains(100)*1
extern std::vector<double> Vab;
extern std::vector<int> cell_ijk;
extern std::vector<int> acceptor_a, acceptor_b;
extern std::vector<int> bonding_a_h, bonding_b_h;
extern std::vector<std::string> H_txt;      // hydrogen-bond diagnostic text (prt_hbonds)
extern std::vector<double> H_energy;        // hydrogen-bond energies (prt_hbonds)

// Globals introduced by add_more_interactions.F90: packed interaction
// vectors (double precision, 1-based) — f already declared above.
extern std::vector<double> p;
extern std::vector<double> pdiag;
extern std::vector<std::vector<double>> c;
extern std::vector<double> pa;
extern std::vector<double> pb;
extern std::vector<std::vector<double>> cb;
extern std::vector<double> bondab;
extern std::vector<double> h;
extern std::vector<double> w;
extern std::vector<std::vector<double>> po;
extern std::vector<double> wk;
extern std::vector<std::vector<int>> pibonds;
extern std::vector<double> q;
extern std::vector<int> ifact;   // triangular offsets (mullik.F90)
extern std::vector<int> i1fact;  // second triangular offsets (setup_mopac_arrays)
extern std::vector<double> dxyz;
extern std::vector<int> hblist;  // (max_h_bonds,10) column-major flat
extern std::vector<double> grad;
extern std::vector<double> gnext1, gmin1, aicorr;
extern std::vector<std::string> simbol;       // symbolic geometry labels (getgeg.F90)
extern std::vector<double> errfn, fb;         // error function / beta Fock (MOZYME)
extern std::vector<int> na_store;             // stored NA connectivity (getgeo IRC/DRC)
extern std::vector<double> xparef;            // reference xparam (setup_mopac_arrays)
extern std::vector<double> ptot2;             // total density copy (setup_mopac_arrays)
extern std::vector<int> nw;                   // atom-pair bookkeeping (setup_mopac_arrays)
extern std::vector<double> fmatrx;            // force matrix, packed lower triangle (force.F90)
extern std::vector<double> T_range, HOF_tot, H_tot, Cp_tot, S_tot;  // thermo.F90
extern std::vector<std::string> pibonds_txt; // user-supplied pi bonds (geout.F90)
}  // namespace common_arrays_C