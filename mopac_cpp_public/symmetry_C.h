// symmetry_C.h — C++ mapping of Fortran module "symmetry_C".
#pragma once
#include <string>
#include <vector>
namespace symmetry_C {
// jndex(k) principal q.n. of irreducible representation of vibration k.
extern std::vector<int> jndex;
// namo(k) names of irreducible representations (e.g. " A  ").
extern std::vector<std::string> namo;
// cub(3,3): cubic symmetry matrix.
extern double cub[4][4];
extern int ielem[21]; // symmetry operations present (1..20)
// elem(1..3,1..3,iplace): point-group operation matrices.
extern double elem[4][4][21];
// jelem(ioper, iatom): image atom of iatom under operation ioper.
extern std::vector<std::vector<int>> jelem;

extern int nclass;
extern int nirred;   // number of irreducible representations (symoir.F90)
extern std::vector<double> group;   // characters of I.R.s, packed by class
extern std::vector<std::string> jx; // I.R. names (symoir.F90)
extern std::vector<int> jy;
// Point-group tables (symmetry_C.F90 DATA statements), 1-based.
extern const int ntbs;                    // 38 character tables
extern std::vector<int> nallop;           // group descriptors (348 entries)
extern std::vector<int> ntab;             // character-table element counts
extern std::vector<int> nallg;            // character-table data (764 entries)
extern std::vector<std::string> allrep;   // group + I.R. names (406 entries)
extern std::vector<int> locdep;
extern std::vector<int> locpar;
extern std::vector<int> idepfn;
extern std::vector<double> depmul;   // symmetry multipliers (geout.F90)
extern int nsym;
extern int igroup;   // point-group index (force.F90)
extern std::vector<std::vector<int>> ipo;
extern double r[10][121];
extern int nent;
extern std::string name; extern std::string state_spin; extern std::string state_Irred_Rep; extern int state_QN;
}