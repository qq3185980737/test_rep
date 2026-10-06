// parameters_C.h — C++ mapping of Fortran module "parameters_C" (MOPAC 2016).
#pragma once
#include <vector>
namespace parameters_C {
// ams(i) — atomic mass of element i (used for hydrogen: ams(1)).
extern std::vector<double> ams;
// xfac(i,j), alpb(i,j) — solvent (ALPB) parameters, 1-based up to 100.
// Filled by alpb_and_xfac_<model> routines; zeroed on entry.
constexpr int N_XFAC = 100;
extern double xfac[N_XFAC + 1][N_XFAC + 1];
extern double alpb[N_XFAC + 1][N_XFAC + 1];

// Analytical-derivative (analyt.F90) parameters. 1-based, size 107.
constexpr int N_PARAM = 107;
extern double alp[N_PARAM + 1];
extern int natorb[N_PARAM + 1];
extern double tore[N_PARAM + 1];
extern double guess1[N_PARAM + 1][5];
extern double guess2[N_PARAM + 1][5];
extern double guess3[N_PARAM + 1][5];
extern double betas[N_PARAM + 1];
extern double betap[N_PARAM + 1];
extern double betad[N_PARAM + 1];
extern double am[N_PARAM + 1];
extern double ad[N_PARAM + 1];
extern double aq[N_PARAM + 1];
extern double dd[N_PARAM + 1];
extern double qq[N_PARAM + 1];

// calpar.F90 derived/raw NDDO parameters (1-based, size 107).
extern int ios[N_PARAM + 1];   // number of s electrons
extern int iop[N_PARAM + 1];   // number of p electrons
extern int iod[N_PARAM + 1];   // number of d electrons
extern double gpp[N_PARAM + 1];
extern double gp2[N_PARAM + 1];
extern double hsp[N_PARAM + 1];
extern double gss[N_PARAM + 1]; extern double gss6sp[N_PARAM + 1]; extern double gssam1sp[N_PARAM + 1]; extern double gssPM3sp[N_PARAM + 1];
extern double gsp[N_PARAM + 1];
extern double zs[N_PARAM + 1];
extern double zp[N_PARAM + 1];
extern double zd[N_PARAM + 1];
// npq(i,1..3) — principal quantum numbers of the s, p, d shells of element i.
extern int npq[N_PARAM + 1][4];
extern double uss[N_PARAM + 1];
extern double upp[N_PARAM + 1];
extern double udd[N_PARAM + 1];
extern double polvol[N_PARAM + 1];
extern double pocord[N_PARAM + 1];
extern double eisol[N_PARAM + 1];
extern double eheat[N_PARAM + 1];   // heat of formation of atom
extern double eheat_sparkles[N_PARAM + 1];  // sparkle heat of formation
extern double f0sd[N_PARAM + 1];
extern double g2sd[N_PARAM + 1];
extern double f0sd_store[N_PARAM + 1];
extern double g2sd_store[N_PARAM + 1];
extern bool main_group[N_PARAM + 1];
extern double zsn[N_PARAM + 1];
extern double zpn[N_PARAM + 1];
extern double zdn[N_PARAM + 1];
extern bool dorbs[N_PARAM + 1];
// ccrep special-case constants (par1..par4)
extern double par1, par2, par3, par4;
extern int ndelec[N_PARAM + 1];

// PM6-DH2 hydrogen-bond acceptor parameters (H_bond_correction_EC_plus_ER).
// 1..6 = N, O(generic), O(acid), O(peptide), O(water), S. Set on first call.
extern double dh2_a_parameters[7];

extern double ddp[7][108];
// REF (PARAM) metadata, parameters_C.F90:301-337.
extern const char* partyp[37 + 1];
extern double defmin[37 + 1];
extern double defmax[37 + 1];
extern double f0dd[N_PARAM + 1];
extern double f2dd[N_PARAM + 1];
extern double f4dd[N_PARAM + 1];
extern double f0pd[N_PARAM + 1];
extern double f2pd[N_PARAM + 1];
extern double g1pd[N_PARAM + 1];
extern double g3pd[N_PARAM + 1];
}
