// MOZYME_C.h — C++ mapping of Fortran module "MOZYME_C" (subset used).
#ifndef MOZYME_C_H
#define MOZYME_C_H
#include <string>
#include <vector>
namespace MOZYME_C {
// afn(j): residue name (3-char) for residue j (1..20).
extern std::vector<std::string> afn;
// atomname(j,k): atom name inside residue j (flattened [j][k]); k in 1..n_add(j).
extern std::vector<std::vector<std::string>> atomname;
// n_add(j): number of named atoms in residue j.
extern std::vector<int> n_add;

// Globals introduced by add_more_interactions.F90.
extern bool lijbo;
extern std::vector<std::vector<int>> nijbo;
extern bool direct;   // direct SCF
extern bool semidr;   // semi-direct SCF
extern bool rapid;    // rapid (partitioned) SCF path
// thresh: MOZYME sparsity threshold (cutoff = thresh*10 in adjvec).
extern double thresh;
extern std::vector<double> partp;  // partitioned p vector
extern std::vector<double> parth; // partitioned h vector
extern std::vector<double> partf; // partitioned f vector

// Globals introduced by atomrs.F90 (residue identification).
extern std::vector<int> nbackb;   // nbackb(4): backbone atoms
extern int jatom;                  // previous backbone N'
extern int iatom;                  // current residue start atom
extern int loop;                   // extra-residue loop index
extern int mxeno;                  // number of extra residues
extern std::vector<std::string> txeno;  // extra-residue names
extern std::vector<std::vector<int>> nxeno; // (4, mxeno) element counts
extern std::vector<std::string> allres;  // residue names, per residue
extern int maxres; extern int nres; extern std::vector<std::vector<int>> bbone; extern std::vector<std::vector<double>> angles; extern std::vector<std::string> allr; extern int k;                 // max residues
extern std::vector<int> ions;      // formal charge per atom
extern std::vector<int> res_start; // first atom of each residue
extern std::vector<int> start_res; // first atom index of each fragment (geochk/ligand)
extern std::vector<int> at_res;    // residue index of each atom (modchg.F90)
extern std::vector<std::string> tyres; // residue type names (len 3)
extern std::vector<std::string> tyr;   // one-letter residue type codes (23)
extern int size_mres;              // number of residue templates
extern bool odd_h;                 // odd-hydrogen warning flag
extern std::vector<int> iorbs;
extern std::vector<double> ws;
extern int Lewis_tot;
extern std::vector<std::vector<int>> Lewis_elem;
extern std::vector<int> iz;
extern std::vector<int> ib;
extern int icharges;     // number of occupied AOs on atom i
extern double cutofs;
extern double pmax;
extern bool use_three_point_extrap;
extern std::vector<int> iij, numij, ijall, iijj;
extern int ij_dim;
extern int morb;
extern std::vector<double> tvec;
extern int ipad2, ipad4;
// shift: energy shift between occupied and virtual sets (diagg2.F90).
extern double shift;
extern int nvirtual, noccupied;
extern double energy_diff;
extern double ovmax;
extern double sumt;
extern int cvir_dim;
extern int icocc_dim;
extern int icvir_dim;
extern int ipad4;
extern int nelred;
extern int noccupied;

extern std::vector<int> isort, nce, ncf, ncocc, ncvir, nnce, nncf, icocc, icvir;
extern std::vector<double> cocc, cvir;
extern int cocc_dim, cvir_dim, icocc_dim, icvir_dim;
extern std::vector<std::vector<double>> geo_1, geo_2;        // (3, numat) 1-based rows
extern std::vector<double> dxyz_1, dxyz_2;                   // (3*numat) 1-based
extern std::vector<double> xparam_1, xparam_2;               // (nvar) 1-based
extern std::vector<double> partf_1, partp_1, parth_1;        // system 1 partition arrays
extern std::vector<double> partf_2, partp_2, parth_2;        // system 2 partition arrays
extern std::vector<double> f_1, p_1, f_2, p_2;               // Fock / density of each system
extern std::vector<int> nbonds_1, nbonds_2;                  // (numat) 1-based
extern std::vector<std::vector<int>> ibonds_1, ibonds_2;     // (15, numat) 1-based rows
extern std::vector<int> icocc_1, icvir_1, ncocc_1, ncvir_1;  // system 1 LMO bookkeeping
extern std::vector<int> nncf_1, nnce_1, ncf_1, nce_1;
extern std::vector<int> iij_1, iijj_1, ijall_1, numij_1, iorbs_1;
extern std::vector<std::vector<int>> nijbo_1;                // (numat,?) 1-based
extern std::vector<double> cocc_1, cvir_1;                   // system 1 LMO coeffs
extern std::vector<int> icocc_2, icvir_2, ncocc_2, ncvir_2;
extern std::vector<int> nncf_2, nnce_2, ncf_2, nce_2;
extern std::vector<int> iij_2, iijj_2, ijall_2, numij_2, iorbs_2;
extern std::vector<std::vector<int>> nijbo_2;
extern std::vector<double> cocc_2, cvir_2;
extern double tiny, sumt, ovmax; extern int ijc;
extern double sumb; extern std::vector<double> fmo;
extern std::vector<std::vector<int>> ifmo;  // ifmo(2,*) sparse-MO index pairs
extern std::vector<int> jopt;
extern int numred;
extern int norred, nelred;                     // reduced atom/orb/electron counts
extern int mode;                               // MOZYME mode flag (set_up_MOZYME_arrays)
extern double refnuc;                        // reference nuclear energy (hcore_for_MOZYME)
extern std::vector<int> kopt;
extern std::vector<double> part_dxyz, p1, p2, p3;
extern std::vector<int> idiag, nfmo;
extern int fmo_dim;
}

// ---- dual-geometry storage (big_swap.F90: store/extract two systems) ----



#endif // MOZYME_C_H
