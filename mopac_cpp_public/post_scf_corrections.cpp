// post_scf_corrections.cpp — C++ translation of post_scf_corrections.F90.
//
// Adds dispersion, hydrogen-bonding, extra H-H repulsion and other
// corrections to improve intermolecular interaction energies/geometries.
// Branch structure follows the Fortran exactly (each method flag is an
// independent if/else-if arm; dh_plus/dh2/dh2x are separate arms).
#include "post_scf_corrections.h"
#include "dftd3.h"
#include "H_bonds4.h"
#include "H_bond_correction_PM6_DH_Dispersion.h"
#include "H_bond_correction_PM6_DH_type.h"
#include "disp_DnX.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include <string>
#include <vector>

using namespace molkst_C;
using namespace common_arrays_C;

void post_scf_corrections(double& correction, bool l_grad) {
    bool prt = ((keywrd.find(" 0SCF ") != std::string::npos ||
                 keywrd.find(" PRT ") != std::string::npos) &&
                keywrd.find(" DISP") != std::string::npos);
    correction = 0.0;
    E_hb = 0.0;
    E_hh = 0.0;
    E_disp = 0.0;
    P_Hbonds = 0;
    int n = numat;
    // Column-major (3, numat+1) workspace consumed by H_bonds4 /
    // energy_corr_hh_rep (their F90 signature is dxyz(3,numat)).
    std::vector<double> g_flat((n + 1) * 3, 0.0);
    // dftd3 consumes a (4, numat+1) vector-of-vectors with 1-based rows.
    std::vector<std::vector<double>> dxyz_vv(4, std::vector<double>(n + 1, 0.0));
    // Global dxyz is atom-major 1-based: component c (1..3) of atom a sits at
    // (a-1)*3 + c.  c=0..2 here is mapped to (a-1)*3 + c + 1 below.
    if (static_cast<int>(dxyz.size()) < n * 3 + 1) dxyz.resize(n * 3 + 1);

    if (method_pm6_d3h4x) {
        correction += dftd3(l_grad, dxyz_vv);
        correction += H_bonds4(l_grad, g_flat.data());
        correction += energy_corr_hh_rep(l_grad, g_flat.data());
        correction += disp_DnX(l_grad);
    } else if (method_pm6_d3h4) {
        correction += dftd3(l_grad, dxyz_vv);
        correction += H_bonds4(l_grad, g_flat.data());
        correction += energy_corr_hh_rep(l_grad, g_flat.data());
    } else if (method_pm6_d3) {
        correction += dftd3(l_grad, dxyz_vv);
    } else if (method_pm6_dh_plus) {
        correction += PM6_DH_Dispersion(l_grad, dxyz);
        correction += PM6_DH_H_bond_corrections(l_grad, prt);
    } else if (method_pm6_dh2) {
        correction += PM6_DH_Dispersion(l_grad, dxyz);
        correction += PM6_DH_H_bond_corrections(l_grad, prt);
    } else if (method_pm6_dh2x) {
        correction += PM6_DH_Dispersion(l_grad, dxyz);
        correction += PM6_DH_H_bond_corrections(l_grad, prt);
        correction += disp_DnX(l_grad);
    } else if (method_pm7_hh) {
        correction += energy_corr_hh_rep(l_grad, g_flat.data());
        correction += PM6_DH_Dispersion(l_grad, dxyz);
        correction += PM6_DH_H_bond_corrections(l_grad, prt);
    } else if (method_pm7_minus) {
        return;
    } else if (method_pm7) {
        correction += PM6_DH_Dispersion(l_grad, dxyz);
        correction += PM6_DH_H_bond_corrections(l_grad, prt);
    }
    // Merge the column-major workspaces back into the global atom-major dxyz.
    if (l_grad) {
        for (int a = 1; a <= n; ++a)
            for (int c = 0; c < 3; ++c) {
                dxyz[(a - 1) * 3 + c + 1] += g_flat[c * (n + 1) + a];
                dxyz[(a - 1) * 3 + c + 1] += dxyz_vv[c + 1][a];
            }
    }
    if (keywrd.find(" SILENT") == std::string::npos) {
        if (prt && P_Hbonds > 0) print_post_scf_corrections();
    }
}
