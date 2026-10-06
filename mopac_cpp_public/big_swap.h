// big_swap.h — C++ translation of MOPAC 2016 "big_swap.F90".
//
// Protein transition-state location module.  Contains the full driver chain:
// Locate_TS_for_Proteins (multi-step TS search with two stationary points),
// lbfgs_TS (dual-geometry L-BFGS optimization), compfg_TS (dual-system
// energy/gradient), big_swap (store/extract per-system MOZYME arrays),
// build_active_site / select_opt (active-site topology), Refine_TS_for_Proteins
// (iterative refinement), and the keyword editor l_control.
#pragma once
#include <string>
#include <vector>

// big_swap(mode, system): mode=0 store / mode=1 extract data for system 1 or 2.
void big_swap(int mode, int system);

// Protein transition-state multi-step driver (run_mopac dispatch "LOCATE-TS").
void Locate_TS_for_Proteins();

// Dual-geometry L-BFGS optimization.  big_xparam holds 2*nvar coordinates
// (system 1 then system 2), 1-based.
void lbfgs_TS(std::vector<double>& big_xparam, int big_nvar, double& escf_tot,
              bool extra_print);

// Evaluate heat of formation (and gradients) of the two systems.
void compfg_TS(const std::vector<double>& big_xparam, bool int_flag, double& escf1,
               double& escf2, bool fulscf, std::vector<double>& big_grad, bool lgrad);

// Calibration reader (unit 33).  Arrays are 1-based.
void get_pars(std::vector<double>& stresses, std::vector<double>& gradients,
              std::vector<double>& relscf, std::vector<double>& cutoff, int& nloop);

// Identify the atoms whose bonding topology differs between the two systems.
// active_site is 1-based (size 200+1), ninsite holds Set1/Set2 sizes.
void build_active_site(std::vector<int>& active_site, std::vector<int>& ninsite);

// Sort the first shell atoms of active_site into ascending order (no side
// effect on the optimization variables; the F90 "if (.false.)" block is dead).
void select_opt(int shell, std::vector<int>& active_site);

// Iterative refinement of the transition state (energy min of inactive atoms,
// gradient min of active site, up to 5 cycles).
void Refine_TS_for_Proteins();

// l_control(txt, nt, mode): add (mode=1) or remove (mode=-1) keyword words
// from keywrd.  If a word already exists when adding, the old occurrence is
// replaced.  nt is the declared length of txt.
void l_control(const std::string& txt, int nt, int mode);

// Two-argument convenience form (nt inferred from the string length).
void l_control(const std::string& txt, int mode);
