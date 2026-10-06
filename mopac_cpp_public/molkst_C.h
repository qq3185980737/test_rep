// molkst_C.h — C++ translation of the Fortran module "molkst_C" (MOPAC 2016).
//
// Only the members referenced by already-translated sources are declared so
// far; the header will be extended as further .F90 files are ported.
#pragma once

#include <string>
#include <vector>

namespace molkst_C {
extern std::string program_name;

// numat — number of real atoms in the system.
// norbs — number of atomic orbitals in the system.
extern int numat;
extern int numat_old;   // number of atoms in the input geometry (geo_ref)
extern int norbs;

// Globals introduced by add_hydrogen_atoms.F90 and its subprograms.
extern int natoms;
extern int lm61;
extern double enuclr;             // number of atoms (incl. pseudo) in input
extern std::string keywrd;     // keyword line
extern double mol_weight;     // molecular weight (AMU)
extern std::vector<int> mers;  // (3) unit-cell replication counts
extern int nopen;
extern double fract;
extern int nclose;
extern int nelecs;
extern int nalpha;
extern int nbeta;
extern bool mozyme;
extern bool method_pm7;
extern bool method_pm7_ts;      // PM7-TS (transition-state) method flag
extern bool N_3_present;        // any N atom with 3 ligands (PM6 sp2 correction)
extern bool Si_O_H_present;     // any Si-O-H substructure (PM7 correction)
extern bool method_pm6;
extern bool method_am1;
extern bool method_mndod; extern bool is_PARAM, sparkle, method_rm1, method_mndo, method_pm3;
extern bool in_house_only;
extern int maxtxt;
extern int nl_atoms;
extern bool isok;
extern int maxatoms;           // allocated upper bound for atoms
extern bool moperr;
extern bool limscf;    // SCF early-exit flag (set by efstr; used by SCF)            // fatal-error flag
extern bool units;             // true if distance units are defined
extern bool Angstroms;         // true if units are Angstroms (default)
extern std::vector<std::string> refkey;
extern std::string koment;
extern std::string title;
extern std::string line;       // scratch line buffer
extern int nbreaks;
extern int ncomments;          // number of comment lines
extern int id;                 // periodic-boundary flag (0 = isolated)
extern int l11, l21, l31;     // supercell search bounds

// Globals introduced by add_more_interactions.F90.
extern int mpack;              // packed-size of interaction arrays
extern bool lgpu;
extern int numcal;             // MOPAC run / job index
extern double clower, cupper, cutofp;
extern int l1u, l2u, l3u, l123;
extern int nvar; extern double pressure; extern double cosine; extern double escf;
extern double tleft, time0; extern int nscf, tdump, iflepo, last;
extern double rjkab1; extern int msdel;
extern int ndep;
extern double efield[4];
extern double E_disp, E_hb, E_hh; extern bool method_PM6_DH2X; extern int N_Hbonds; extern bool method_PM6; extern bool method_pm6_d3_not_h4;
extern int n2elec, ispd, itemp_2;
extern int itemp_1;   // alias: current optimization / MD cycle index (Fortran jloop)
extern int step_num;  // current step number (write_trajectory)
extern std::string allkey;  // full keyword block (wrtkey input)
extern bool rhf, uhf;       // RHF / UHF spin flags (wrtchk)
extern std::string geo_dat_name, geo_ref_name;  // GEO_DAT / GEO_REF file names
extern double step;
extern double Rab;
extern int P_Hbonds;
extern bool method_pm6_d3h4x, method_pm6_d3h4, method_pm6_d3;
extern bool method_pm6_dh_plus, method_pm6_dh2, method_pm6_dh2x;
extern bool method_pm7_hh, method_pm7_minus;
extern double ux, uy, uz;
extern int na1;
extern std::string jobnam;
extern bool gui;
extern std::string errtxt;
// writmo results-summary globals.
extern bool use_ref_geo;
extern double stress, density, elect, atheat;
extern std::string verson;
extern int ijulian, iscf;
extern double hpress, nsp2_corr, Si_O_H_corr, sum_dihed;
extern double emin;              // lowest energy found (compfg.F90)
extern int nalpha_open, nbeta_open;
extern bool prt_gradients, prt_coords, prt_cart, prt_pops, prt_charges;
extern bool prt_topo;            // print atom-topography table (lewis.F90)
extern double gnorm;       // gradient norm (force.F90)
extern double zpe;         // zero-point energy (force.F90)
extern std::string formula;
extern double sz, ss2;
extern int ltxt;                    // length of longest atom text label (getgeg/geoutg)
extern double temp_1, temp_2;       // cp298 / s298 from thermo (Fortran aliases)
extern int ilim;                    // number of temperatures in thermo table
extern int num_threads;             // MKL thread count (diag_for_GPU)
extern double param_constant;   // run_mopac: overall scaling of parameters
extern double trunc_1, trunc_2; // run_mopac: point-charge truncation limits
extern bool lxfac;              // run_mopac: XFAC keyword
extern bool lgpu;               // run_mopac: GPU flag
extern bool in_house_only;      // run_mopac: in-house-only run flag
extern bool MM_corrections;     // run_mopac: MM corrections flag
extern int no_pKa;              // run_mopac: number of pKa values
extern int run;                     // command-line argument index (getdat)
extern double arc_hof_1, arc_hof_2; // archive heats of formation (getdat)
extern double tleft;                // time remaining (grid)
extern int ijulian;                 // days since 1-Jan-2009 (negative until password() runs)
extern int site_no;                 // license site number (password.F90)
extern bool academic;               // academic vs commercial license
extern std::string verson;          // version string, e.g. "2016    " (password.F90 writes [6])
extern int num_bits;                // 32 or 64 (setup_mopac_arrays hint)

}  // namespace molkst_C
