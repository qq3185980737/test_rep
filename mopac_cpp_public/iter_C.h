// iter_C.h — C++ mapping of Fortran module "iter_C" (subset used).
#pragma once
#include <vector>
namespace iter_C {
// pold(:): previous-iteration density/packed vector (1-based padding).
extern std::vector<double> pold;
extern std::vector<double> pold2, pold3;
extern std::vector<double> pbold, pbold2, pbold3;
// Camp-King converger scratch (norbs x norbs, norbs*norbs)
extern std::vector<double> vec_ai, fock_ai, p_ai, h_ai, vecl_ai;
extern std::vector<double> vec_bi, fock_bi, p_bi, h_bi, vecl_bi;
extern double eold_alpha, eold_beta;
extern std::vector<double> theta;
extern bool md_test, is_PARAM;
extern bool lgpu;}
