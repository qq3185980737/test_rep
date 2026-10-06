// linear_cosmo.h — C++ translation of MOPAC 2016 "linear_cosmo.F90".
// module cosmo_mini data + module linear_cosmo public interface.
// Batch A implements: ini_linear_cosmo, bpnew_vec, bz_vec, get_bvec,
// mult_triangle_vec, some_norm. Remaining subroutines declared here are
// implemented in later batches.
#pragma once
#include <vector>
#include "afmm_mod.h"

namespace linear_cosmo {

// ---- module cosmo_mini data (1-based F90 logic maps directly) ----
extern bool new_surface, new_iteration;
extern std::vector<double> a_block, r_vec, p_vec, q_vec, z_vec, a_diag,
                           m_vec, a_part;
extern std::vector<double> tm;      // flat (3*3*numat), F90 tm(3,3,numat)
extern std::vector<int> iblock_pos;

// ---- module linear_cosmo data ----
extern double c_proc;               // F90: double precision, public :: c_proc
extern std::vector<double> rsc;     // flat (4*maxrs)
extern std::vector<int> nipsrs, nset, npoints, iatom_pos, ijbo_diag;
extern int max_block_size, atom_handle, surface_handle;

// F90 integer, parameter
constexpr int na1max = 6000;        // O(n*log(n)) threshold for A*q
constexpr int na2max = 6000;        // O(n) threshold for A*q
constexpr int nb1max = 8000;        // O(n*log(n)) threshold for B*q
constexpr int nb2max = 8000;        // O(n) threshold for B*q
constexpr bool compute_a_part = true;   // precompute A matrix
constexpr bool use_a_blocks = true;     // block-diagonal preconditioner

// ---- subroutines (batch A: implemented) ----
void ini_linear_cosmo();
void bpnew_vec(std::vector<double>& v);                 // interaction of electron densities with induced charges
void bz_vec(std::vector<double>& v);                    // interaction of atom cores with induced charges
void get_bvec(const double* x1, const double* x2, int nao, int ni, double* w);
void mult_triangle_vec(const double* t, const double* x, int n, double* y);
double some_norm(const double* v, int n);

// ---- subroutines (batch B: coscavz + coscanz surface construction) ----
void coscavz();                       // F90 coscavz(coord, nat): driver of coscanz + tessellation
void coscavz(std::vector<double>& c, int nat);  // F90 signature: flat coord(3*nat)
void coscanz();                       // F90 coscanz(srad, coord, nat, cosurf, iatsp, nar_loc, nsetf, phinet, arat)
void simulate_aq_vec(const std::vector<std::vector<double>>& coord,
                     const std::vector<double>& srad, int numat1,
                     const std::vector<std::vector<double>>& cosurf, int nps1,
                     const std::vector<int>& nar_csm, const std::vector<int>& nsetf,
                     const std::vector<int>& nset, std::vector<double>& rsc,
                     std::vector<int>& nipsrs, const std::vector<int>& iatsp,
                     const std::vector<double>& tm, int ioldcv, int maxrs,
                     int lenabc, bool compute, int& ndim);  // F90 simulate_aq_vec
void simulate_aq_dir_int(const int* ind1, int n1, const int* ind2, int n2,
                         const double* c, int nc, bool same, bool compute, int& n);

// ---- batch C: A-part machinery ----
void aq_vec(const std::vector<std::vector<double>>& coord,
            const std::vector<double>& srad, int numat1,
            const std::vector<std::vector<double>>& cosurf, int nps1,
            const std::vector<int>& nar_csm, const std::vector<int>& nsetf,
            const std::vector<int>& nset, std::vector<double>& rsc,
            std::vector<int>& nipsrs, const std::vector<int>& iatsp,
            const std::vector<double>& tm, int ioldcv, int maxrs, int lenabc,
            const std::vector<double>& a_diag, const std::vector<double>& q,
            std::vector<double>& v);                                // F90 aq_vec @1159
void aq_dir_int(const int* ind1, int n1, const int* ind2, int n2,
                const double* c, int nc, const double* q, double* r,
                bool same);                                         // F90 aq_dir_int @1342

// ---- batch D: cgm_solve solver family ----
void precond_aq_dir_int(const int* ind1, int n1, const int* ind2, int n2,
                        const double* c, int nc, bool same, bool compute, int& npos);
void aq_mult_int(const int* ind, int n, const afmm_mod::CArr& psi, int p,
                 double x0, double y0, double z0, const afmm_mod::RArr& y_norm,
                 afmm_mod::RArr& pmn, const double* c, int n3,
                 const double* q, double* r);
void aq_far_int(const int* ind, int n, const afmm_mod::CArr& psi, int p,
                double x0, double y0, double z0, const afmm_mod::RArr& y_norm,
                afmm_mod::RArr& pmn, const double* c, int nc,
                const double* q, double* r);
void amat_diag(const std::vector<std::vector<double>>& coord,
               const std::vector<double>& srad, int numat,
               const std::vector<std::vector<double>>& cosurf, int nps,
               const std::vector<int>& nar_csm, const std::vector<int>& nsetf,
               const std::vector<int>& nset, std::vector<double>& rsc,
               std::vector<int>& nipsrs, const std::vector<int>& iatsp,
               const std::vector<double>& tm, int ioldcv, int maxrs,
               int lenabc, std::vector<double>& m);
void precondition(const std::vector<std::vector<double>>& cosurf, int nps,
                  const std::vector<int>& iatsp, int numat,
                  const std::vector<double>& a_diag, std::vector<double>& m);
void precondition_solve(const std::vector<double>& m, std::vector<double>& z,
                        const std::vector<double>& r, int n);
void cgm_solve(std::vector<double>& x, int n, const std::vector<double>& b,
               std::vector<double>& r, std::vector<double>& p,
               std::vector<double>& q, std::vector<double>& z,
               std::vector<double>& m, bool new_x, bool new_points);

// ---- batch E: BZ/BP expansion family ----
void bz_far_int(const int* ind, int n, const afmm_mod::CArr& psi, int p,
                double x0, double y0, double z0, const afmm_mod::RArr& y_norm,
                afmm_mod::RArr& pmn, const double* c, int nc,
                const double* q, double* r);        // F90 bz_far_int @1828
void bz_mult_int(const int* ind, int n, const afmm_mod::CArr& psi, int p,
                 double x0, double y0, double z0, const afmm_mod::RArr& y_norm,
                 afmm_mod::RArr& pmn, const double* c, int nc,
                 const double* q, double* r);       // F90 bz_mult_int @1885
void bp_dir_int(const int* ind1, int n1, const int* ind2, int n2,
                const double* c, int nc, const double* q, double* r,
                bool same);                         // F90 bp_dir_int @1951

// ---- batch F: nuclear/Fock corrections (addnucz/addfckz) ----
void addnucz(std::vector<std::vector<double>>& phinet,
             std::vector<std::vector<double>>& qscnet,
             std::vector<std::vector<double>>& qdenet);   // F90 addnucz @2544
void addfckz();                                           // F90 addfckz @2573
// ---- Fock far-field family (stubs in batch F, to be translated) ----
void fock_dir_int(const int* ind1, int n1, const int* ind2, int n2,
                  const double* c, int nc, const double* q, double* r,
                  bool same);                             // F90 fock_dir_int @2747
void fock_mult_int(const int* ind, int n, const afmm_mod::CArr& psi, int p,
                   double x0, double y0, double z0, const afmm_mod::RArr& y_norm,
                   afmm_mod::RArr& pmn, const double* c, int nc,
                   const double* q, double* r);           // F90 fock_mult_int @2850
void fock_far_int(const int* ind, int n, const afmm_mod::CArr& psi, int p,
                  double x0, double y0, double z0, const afmm_mod::RArr& y_norm,
                  afmm_mod::RArr& pmn, const double* c, int nc,
                  const double* q, double* r);            // F90 fock_far_int @2923
void sphere_f_multipoles(double x0, double y0, double z0, const int* ind, int n,
                         afmm_mod::RArr& pmn, const afmm_mod::RArr& y_norm,
                         const double* c, int nc, const double* q,
                         afmm_mod::CArr& psi, int p);     // F90 sphere_f_multipoles @2997
void am1dft_solve(std::vector<double>& qsct, std::vector<double>& phit, int nps);

}  // namespace linear_cosmo
