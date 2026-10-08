// md_driver.cpp — C++ port of the user's Fortran velocity-Verlet MD driver.
//
// Faithful translation of "subroutine run_md_velocity_verlet" and
// "subroutine gen_md_vel" from C:\Users\15632\Desktop\md.txt, operating on the
// existing C++ globals. Two deliberate fixes vs the Fortran original:
//   1. The box size is a real keyword-controlled variable (default 10.0 A),
//      and the printed box always matches the box actually used (the Fortran
//      hard-coded the print as "[0.0, 40.0]" while using 10 A).
//   2. T_target / dt / box can be set by TEMP=, DT=, BOX= keywords instead of
//      only by editing the global defaults.
// Momentum / angular-momentum removal stay disabled, exactly as in the Fortran
// original (they are commented out there).
#include "md_driver.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "common_arrays_C.h"
#include "compfg.h"
#include "cosmo_C.h"
#include "funcon_C.h"
#include "maps_C.h"
#include "md_C.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "reada.h"

using namespace common_arrays_C;   // nat, loc, xparam
using namespace cosmo_C;           // solv_energy
using namespace funcon_C;          // fpc_9 (1 eV in kcal/mol)
using namespace maps_C;            // react (initial-velocity vector from VELO)
using namespace md_C;              // md_vel, md_force, md_vel_half, dt_md, T_target, md_step_max, md_stop
using namespace molkst_C;          // numat, nvar, escf, elect, enuclr
using namespace parameters_C;      // ams

namespace {
// Atomic-unit conversion constants (identical to the Fortran reference).
constexpr double BOHR_TO_ANG   = 0.529177210903;
constexpr double ANG_TO_BOHR   = 1.0 / BOHR_TO_ANG;
constexpr double HARTREE_TO_EV = 27.21073973;
constexpr double EV_TO_HARTREE = 1.0 / HARTREE_TO_EV;
constexpr double AU_TIME_FS    = 0.024188843265;   // 1 atomic unit of time in fs
constexpr double KB_AU         = 3.166811563e-6;   // Hartree / K (md_C::KB)
constexpr double U_TO_ME       = 1822.888486;      // 1 atomic mass unit in electron masses
}  // namespace

// ---------------------------------------------------------------------------
// gen_md_vel: build the initial MD velocities.
//  - VELO keyword present + velocity block read: use the user's velocities
//    (Angstroms/fs in maps_C::react), then zero the net momentum and rescale
//    to T_target — both DRC-referenced (drc.F90 L290-329, L378-381). This keeps
//    the user's relative directions while setting the temperature.
//  - otherwise: Maxwell-Boltzmann at T_target (Box-Muller), rescaled to T_target.
// ---------------------------------------------------------------------------
static void gen_md_vel(const std::string& keywrd, int iw) {
  md_vel.assign(4, std::vector<double>(numat + 1, 0.0));
  const bool use_velo  = (keywrd.find(" VELO") != std::string::npos);
  const bool have_velo = use_velo &&
                         (static_cast<int>(react.size()) >= 3 * numat + 1);

  if (have_velo) {
    // 1 Å/fs = AU_TIME_FS / BOHR_TO_ANG  Bohr/au_time.
    for (int ia = 1; ia <= numat; ++ia)
      for (int c = 1; c <= 3; ++c)
        md_vel[c][ia] = react[(ia - 1) * 3 + c] * AU_TIME_FS / BOHR_TO_ANG;
    // Zero the net momentum: subtract the mass-weighted mean velocity.
    double mx = 0.0, my = 0.0, mz = 0.0, tot = 0.0;
    for (int ia = 1; ia <= numat; ++ia) {
      double m = ams[nat[ia]] * U_TO_ME;
      tot += m; mx += m * md_vel[1][ia]; my += m * md_vel[2][ia]; mz += m * md_vel[3][ia];
    }
    if (tot > 0.0) { mx /= tot; my /= tot; mz /= tot; }
    for (int ia = 1; ia <= numat; ++ia) {
      md_vel[1][ia] -= mx; md_vel[2][ia] -= my; md_vel[3][ia] -= mz;
    }
  } else {
    if (use_velo)
      std::fprintf(stdout,
                   " WARNING: VELO requested but no velocity block read; using random.\n");
    // Maxwell-Boltzmann velocities at T_target (Box-Muller, seed 12345).
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    for (int ia = 1; ia <= numat; ++ia) {
      double mass = ams[nat[ia]] * U_TO_ME;          // electron masses
      double sd   = std::sqrt(KB_AU * T_target / mass);
      for (int c = 1; c <= 3; ++c) {
        double u1 = uni(rng), u2 = uni(rng);
        u1 = std::max(1.0e-12, u1);
        u1 = std::min(1.0 - 1.0e-12, u1);
        double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * 3.14159265358979323846 * u2);
        md_vel[c][ia] = z * sd;
      }
    }
  }

  // Rescale (both paths) so the instantaneous temperature equals T_target.
  int f_vib = 3 * numat;
  double Ek = 0.0;
  for (int ia = 1; ia <= numat; ++ia) {
    double mass = ams[nat[ia]] * U_TO_ME;
    Ek += 0.5 * mass * (md_vel[1][ia] * md_vel[1][ia] +
                        md_vel[2][ia] * md_vel[2][ia] +
                        md_vel[3][ia] * md_vel[3][ia]);
  }
  double T_calc = 2.0 * Ek / (f_vib * KB_AU);

  std::fprintf(stdout, "==== 温度信息 ====\n");
  std::fprintf(stdout, "目标温度 T_target = %10.4f\n", T_target);
  std::fprintf(stdout, "计算温度 T_calc   = %10.4f\n", T_calc);
  if (T_calc > 0.0) {
    double lambda = std::sqrt(T_target / T_calc);
    for (int ia = 1; ia <= numat; ++ia)
      for (int c = 1; c <= 3; ++c) md_vel[c][ia] *= lambda;
    std::fprintf(stdout, "速度缩放因子 = %10.4f\n", lambda);
    Ek = 0.0;
    for (int ia = 1; ia <= numat; ++ia) {
      double mass = ams[nat[ia]] * U_TO_ME;
      Ek += 0.5 * mass * (md_vel[1][ia] * md_vel[1][ia] +
                          md_vel[2][ia] * md_vel[2][ia] +
                          md_vel[3][ia] * md_vel[3][ia]);
    }
    T_calc = 2.0 * Ek / (f_vib * KB_AU);
    std::fprintf(stdout, "最终温度 T_final = %10.4f\n", T_calc);
  } else {
    std::fprintf(stdout, "警告: T_calc = 0，跳过速度缩放\n");
  }
  if (use_velo)
    std::fprintf(stdout,
                 " (初始速度来自 VELO 速度块，已清零净动量并按目标温度缩放)\n");
}

// ---------------------------------------------------------------------------
// run_md_velocity_verlet: velocity-Verlet MD with elastic-wall box.
// ---------------------------------------------------------------------------
void run_md_velocity_verlet(const std::string& keywrd, int iw) {
  // ---- parse control keywords ----
  // T_target (K), default DEF_TTARGET
  T_target = DEF_TTARGET;
  {
    std::size_t p = keywrd.find("TEMP=");
    if (p == std::string::npos) p = keywrd.find("TTARGET=");
    if (p != std::string::npos) T_target = reada(keywrd, static_cast<int>(p) + 1);
  }
  // dt_md (fs), default DEF_DTMD
  dt_md = DEF_DTMD;
  {
    std::size_t p = keywrd.find("DT=");
    if (p != std::string::npos) dt_md = reada(keywrd, static_cast<int>(p) + 1);
  }
  // number of MD steps, default DEF_MDSTEP
  md_step_max = DEF_MDSTEP;
  {
    std::size_t p = keywrd.find("MDSTEPS=");
    if (p == std::string::npos) p = keywrd.find("CYCLES=");
    if (p != std::string::npos) md_step_max = static_cast<int>(reada(keywrd, static_cast<int>(p) + 1) + 0.5);
  }
  // box [0, box_max] in Angstrom, default 10.0 (matches the Fortran constant)
  double box_min_a = 0.0;
  double box_max_a = 10.0;
  {
    std::size_t p = keywrd.find("BOX=");
    if (p != std::string::npos) box_max_a = reada(keywrd, static_cast<int>(p) + 1);
  }
  const double box_min_bohr = box_min_a * ANG_TO_BOHR;
  const double box_max_bohr = box_max_a * ANG_TO_BOHR;

  std::vector<double> x_new = xparam;   // 1-based, size nvar+1
  std::vector<double> grad(nvar + 1, 0.0);
  const double dt = dt_md / AU_TIME_FS;            // fs -> au_time

  // ---- helper to extract the (x,y,z) of one atom from x_new (in Angstrom) ----
  auto atom_pos = [&](std::vector<double>& v, int ia) {
    double x = 0.0, y = 0.0, z = 0.0;
    for (int i = 1; i <= nvar; ++i) {
      if (loc[1][i] != ia) continue;
      if (loc[2][i] == 1) x = v[i];
      else if (loc[2][i] == 2) y = v[i];
      else if (loc[2][i] == 3) z = v[i];
    }
    return std::vector<double>{x, y, z};
  };

  // ---- print initial coordinates (in Bohr) ----
  std::fprintf(stdout, " 初始坐标\n");
  std::fprintf(stdout, " Atom       X(au)            Y(au)            Z(au)\n");
  for (int ia = 1; ia <= numat; ++ia) {
    std::vector<double> p = atom_pos(x_new, ia);
    std::fprintf(stdout, "%4d %16.6E %16.6E %16.6E\n", ia,
                 p[0] * ANG_TO_BOHR, p[1] * ANG_TO_BOHR, p[2] * ANG_TO_BOHR);
  }
  std::fprintf(stdout, "===============================================\n\n");

  gen_md_vel(keywrd, iw);
  std::fprintf(stdout, ">>>>> END INIT md_vel <<<<<\n");
  std::fprintf(stdout, "\n==================== MD SIMULATION START ====================\n");
  std::fprintf(stdout, " Total MD steps      = %6d\n", md_step_max);
  std::fprintf(stdout, " Timestep DT         = %8.4f au_times\n", dt);
  std::fprintf(stdout, " Target temperature  = %8.2f K\n", T_target);
  std::fprintf(stdout, " Simulation Box      = [%.1f, %.1f] (Elastic Wall)\n",
               box_min_a, box_max_a);
  std::fprintf(stdout, "============================================================\n");

  // velocity Verlet main loop
  for (int istep = 1; istep <= md_step_max; ++istep) {
    std::vector<std::vector<int>> hit_flag(4, std::vector<int>(numat + 1, 0));

    // ---- forces at time t ----
    compfg(x_new, true, escf, true, grad, true);
    if (md_stop) break;

    md_force.assign(4, std::vector<double>(numat + 1, 0.0));
    for (int i = 1; i <= nvar; ++i)
      md_force[loc[2][i]][loc[1][i]] =
          -grad[i] * BOHR_TO_ANG / (HARTREE_TO_EV * fpc_9);

    // ---- positions in Bohr ----
    std::vector<std::vector<double>> r_bohr(4, std::vector<double>(numat + 1, 0.0));
    for (int ia = 1; ia <= numat; ++ia) {
      std::vector<double> p = atom_pos(x_new, ia);
      r_bohr[1][ia] = p[0] * ANG_TO_BOHR;
      r_bohr[2][ia] = p[1] * ANG_TO_BOHR;
      r_bohr[3][ia] = p[2] * ANG_TO_BOHR;
    }

    // ---- 1. half-step velocity ----
    md_vel_half.assign(4, std::vector<double>(numat + 1, 0.0));
    for (int ia = 1; ia <= numat; ++ia) {
      double mass = ams[nat[ia]] * U_TO_ME;
      for (int c = 1; c <= 3; ++c)
        md_vel_half[c][ia] = md_vel[c][ia] + md_force[c][ia] / mass * dt * 0.5;
    }

    // ---- 2. position update in Bohr + elastic-wall reflection ----
    for (int ia = 1; ia <= numat; ++ia) {
      double xb = r_bohr[1][ia] + md_vel_half[1][ia] * dt;
      double yb = r_bohr[2][ia] + md_vel_half[2][ia] * dt;
      double zb = r_bohr[3][ia] + md_vel_half[3][ia] * dt;
      if (xb > box_max_bohr) { xb = box_max_bohr - (xb - box_max_bohr); hit_flag[1][ia] = 1; }
      else if (xb < box_min_bohr) { xb = box_min_bohr + (box_min_bohr - xb); hit_flag[1][ia] = 1; }
      if (yb > box_max_bohr) { yb = box_max_bohr - (yb - box_max_bohr); hit_flag[2][ia] = 1; }
      else if (yb < box_min_bohr) { yb = box_min_bohr + (box_min_bohr - yb); hit_flag[2][ia] = 1; }
      if (zb > box_max_bohr) { zb = box_max_bohr - (zb - box_max_bohr); hit_flag[3][ia] = 1; }
      else if (zb < box_min_bohr) { zb = box_min_bohr + (box_min_bohr - zb); hit_flag[3][ia] = 1; }
      r_bohr[1][ia] = xb; r_bohr[2][ia] = yb; r_bohr[3][ia] = zb;
    }
    // write Bohr coords back to x_new (Angstrom) for the next compfg
    for (int ia = 1; ia <= numat; ++ia) {
      for (int i = 1; i <= nvar; ++i) {
        if (loc[1][i] != ia) continue;
        if (loc[2][i] == 1) x_new[i] = r_bohr[1][ia] * BOHR_TO_ANG;
        else if (loc[2][i] == 2) x_new[i] = r_bohr[2][ia] * BOHR_TO_ANG;
        else if (loc[2][i] == 3) x_new[i] = r_bohr[3][ia] * BOHR_TO_ANG;
      }
    }

    // ---- forces at time t+dt ----
    compfg(x_new, true, escf, true, grad, true);
    if (md_stop) break;

    std::fprintf(stdout, " 新坐标重新计算受力\n");
    std::fprintf(stdout, " FUNCTION (Heat of formation, kcal/mol) = %16.8E\n", escf);
    std::fprintf(stdout, " Atom       X(A)            Y(A)            Z(A)\n");
    for (int ia = 1; ia <= numat; ++ia) {
      std::vector<double> p = atom_pos(x_new, ia);
      std::fprintf(stdout, "%4d %20.8E %20.8E %20.8E\n", ia, p[0], p[1], p[2]);
    }

    std::vector<std::vector<double>> f_atom(4, std::vector<double>(numat + 1, 0.0));
    for (int i = 1; i <= nvar; ++i)
      f_atom[loc[2][i]][loc[1][i]] =
          -grad[i] * BOHR_TO_ANG / (HARTREE_TO_EV * fpc_9);

    // ---- 3. full velocity update ----
    for (int ia = 1; ia <= numat; ++ia) {
      double mass = ams[nat[ia]] * U_TO_ME;
      for (int c = 1; c <= 3; ++c)
        md_vel[c][ia] = md_vel_half[c][ia] + f_atom[c][ia] / mass * dt * 0.5;
    }
    // elastic wall: reflect the (already-updated) velocity component
    for (int ia = 1; ia <= numat; ++ia)
      for (int c = 1; c <= 3; ++c)
        if (hit_flag[c][ia]) md_vel[c][ia] = -md_vel[c][ia];

    // ---- instantaneous temperature ----
    int f_vib = 3 * numat;
    double Ek = 0.0;
    for (int ia = 1; ia <= numat; ++ia) {
      double mass = ams[nat[ia]] * U_TO_ME;
      Ek += 0.5 * mass * (md_vel[1][ia] * md_vel[1][ia] +
                          md_vel[2][ia] * md_vel[2][ia] +
                          md_vel[3][ia] * md_vel[3][ia]);
    }
    double T_inst = 2.0 * Ek / (f_vib * KB_AU);

    std::fprintf(stdout, "\n ---- MD Step %6d ----\n", istep);
    std::fprintf(stdout, " Potential Energy (MD mechanical energy, Hartree) = %12.4f\n",
                 (elect + enuclr + solv_energy) * EV_TO_HARTREE);
    std::fprintf(stdout, " Total Energy (Potential + Kinetic)原子单位能量 = %12.4f\n",
                 (elect + enuclr + solv_energy) * EV_TO_HARTREE + Ek);
    std::fprintf(stdout, " Instant Temperature = %10.2f K\n", T_inst);
    std::fprintf(stdout, " Total Kinetic Energy = %12.6f\n", Ek);

    xparam = x_new;
  }

  std::fprintf(stdout, "\n==================== MD SIMULATION FINISHED ====================\n");
}
