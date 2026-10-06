// parameters_C.cpp — storage for parameters_C.
#include "parameters_C.h"
namespace parameters_C {
// tore(i) = atomic number Z(i) (effective nuclear charge). MOPAC fills these
// from the parameter data sets; for main-group and transition elements Z=Z.
struct ToreInit {
    // MOPAC fills tore(i) = ios(i) + iop(i) + iod(i)  (run_mopac.F90 main:
    // "tore = ios + iop + iod"), the s/p/d shell occupancies from the
    // parameters_C data tables. Then run_mopac.F90:349 sets tore(i)=3 for
    // lanthanides 57..71 with zs(i)<0.1; approximated here as all 57..71.
    ToreInit() {
        for (int i = 1; i <= N_PARAM; ++i) tore[i] = 0.0;
        for (int i = 1; i <= N_PARAM; ++i) tore[i] = ios[i] + iop[i] + iod[i];
        for (int i = 57; i <= 71; ++i) tore[i] = 3.0;   // lanthanide special (zs<0.1)
    }
} _tore_init;
std::vector<double> ams = {
    0.0,  // unused index 0
    1.0079, 4.0026, 6.94, 9.01218, 10.81, 12.011, 14.0067, 15.9994, 
    18.9984, 20.179, 22.98977, 24.305, 26.98154, 28.0855, 30.97376, 32.06, 
    35.453, 39.948, 39.0983, 40.08, 44.9559, 47.9, 50.9415, 51.996, 
    54.938, 55.847, 58.9332, 58.71, 63.546, 65.38, 69.735, 72.59, 
    74.9216, 78.96, 79.904, 83.8, 85.4678, 87.62, 88.9059, 91.22, 
    92.9064, 95.94, 98.9062, 101.07, 102.9055, 106.4, 107.868, 112.41, 
    114.82, 118.69, 121.75, 127.6, 126.9045, 131.3, 132.9054, 137.33, 
    138.906, 140.116, 140.9077, 144.24, 145, 150.36, 151.964, 157.25, 
    158.9253, 162.5, 164.9303, 167.26, 168.9342, 173.04, 174.967, 178.49, 
    180.9479, 183.85, 186.207, 190.2, 192.22, 195.09, 196.9665, 200.59, 
    204.37, 207.2, 208.9804, 209, 210, 222, 223, 226, 
    227, 232.0381, 231.0359, 238.0289, 0.0, 0.0, 0.0, 0.0, 
    0.0, 5e-05, 0.0, 0.0, 0.0, 1.0079, 0.0, 0.0, 
    0.0, 0.0, 0.0, 
};

double xfac[N_XFAC + 1][N_XFAC + 1] = {};
double alpb[N_XFAC + 1][N_XFAC + 1] = {};
double alp[N_PARAM + 1] = {};
int natorb[N_PARAM + 1] = {};
double tore[N_PARAM + 1] = {};
double guess1[N_PARAM + 1][5] = {};
double guess2[N_PARAM + 1][5] = {};
double guess3[N_PARAM + 1][5] = {};
double betas[N_PARAM + 1] = {};
double betap[N_PARAM + 1] = {};
double betad[N_PARAM + 1] = {};
double am[N_PARAM + 1] = {};
double ad[N_PARAM + 1] = {};
double aq[N_PARAM + 1] = {};
double dd[N_PARAM + 1] = {};
double qq[N_PARAM + 1] = {};
int ios[N_PARAM + 1] = { 0,  // [0] unused
      1,   // ios(1) F90 parameters_C data
      2,   // ios(2) F90 parameters_C data
      1,   // ios(3) F90 parameters_C data
      2,   // ios(4) F90 parameters_C data
      2,   // ios(5) F90 parameters_C data
      2,   // ios(6) F90 parameters_C data
      2,   // ios(7) F90 parameters_C data
      2,   // ios(8) F90 parameters_C data
      2,   // ios(9) F90 parameters_C data
      0,   // ios(10) F90 parameters_C data
      1,   // ios(11) F90 parameters_C data
      2,   // ios(12) F90 parameters_C data
      2,   // ios(13) F90 parameters_C data
      2,   // ios(14) F90 parameters_C data
      2,   // ios(15) F90 parameters_C data
      2,   // ios(16) F90 parameters_C data
      2,   // ios(17) F90 parameters_C data
      0,   // ios(18) F90 parameters_C data
      1,   // ios(19) F90 parameters_C data
      2,   // ios(20) F90 parameters_C data
      2,   // ios(21) F90 parameters_C data
      2,   // ios(22) F90 parameters_C data
      2,   // ios(23) F90 parameters_C data
      1,   // ios(24) F90 parameters_C data
      2,   // ios(25) F90 parameters_C data
      2,   // ios(26) F90 parameters_C data
      2,   // ios(27) F90 parameters_C data
      2,   // ios(28) F90 parameters_C data
      1,   // ios(29) F90 parameters_C data
      2,   // ios(30) F90 parameters_C data
      2,   // ios(31) F90 parameters_C data
      2,   // ios(32) F90 parameters_C data
      2,   // ios(33) F90 parameters_C data
      2,   // ios(34) F90 parameters_C data
      2,   // ios(35) F90 parameters_C data
      0,   // ios(36) F90 parameters_C data
      1,   // ios(37) F90 parameters_C data
      2,   // ios(38) F90 parameters_C data
      2,   // ios(39) F90 parameters_C data
      2,   // ios(40) F90 parameters_C data
      1,   // ios(41) F90 parameters_C data
      1,   // ios(42) F90 parameters_C data
      2,   // ios(43) F90 parameters_C data
      1,   // ios(44) F90 parameters_C data
      1,   // ios(45) F90 parameters_C data
      0,   // ios(46) F90 parameters_C data
      1,   // ios(47) F90 parameters_C data
      2,   // ios(48) F90 parameters_C data
      2,   // ios(49) F90 parameters_C data
      2,   // ios(50) F90 parameters_C data
      2,   // ios(51) F90 parameters_C data
      2,   // ios(52) F90 parameters_C data
      2,   // ios(53) F90 parameters_C data
      0,   // ios(54) F90 parameters_C data
      1,   // ios(55) F90 parameters_C data
      2,   // ios(56) F90 parameters_C data
      2,   // ios(57) F90 parameters_C data
      0,   // ios(58) F90 parameters_C data
      0,   // ios(59) F90 parameters_C data
      0,   // ios(60) F90 parameters_C data
      0,   // ios(61) F90 parameters_C data
      0,   // ios(62) F90 parameters_C data
      2,   // ios(63) F90 parameters_C data
      2,   // ios(64) F90 parameters_C data
      2,   // ios(65) F90 parameters_C data
      2,   // ios(66) F90 parameters_C data
      2,   // ios(67) F90 parameters_C data
      2,   // ios(68) F90 parameters_C data
      2,   // ios(69) F90 parameters_C data
      2,   // ios(70) F90 parameters_C data
      2,   // ios(71) F90 parameters_C data
      2,   // ios(72) F90 parameters_C data
      2,   // ios(73) F90 parameters_C data
      1,   // ios(74) F90 parameters_C data
      2,   // ios(75) F90 parameters_C data
      2,   // ios(76) F90 parameters_C data
      2,   // ios(77) F90 parameters_C data
      1,   // ios(78) F90 parameters_C data
      1,   // ios(79) F90 parameters_C data
      2,   // ios(80) F90 parameters_C data
      2,   // ios(81) F90 parameters_C data
      2,   // ios(82) F90 parameters_C data
      2,   // ios(83) F90 parameters_C data
      2,   // ios(84) F90 parameters_C data
      2,   // ios(85) F90 parameters_C data
      0,   // ios(86) F90 parameters_C data
      1,   // ios(87) F90 parameters_C data
      1,   // ios(88) F90 parameters_C data
      2,   // ios(89) F90 parameters_C data
      4,   // ios(90) F90 parameters_C data
      2,   // ios(91) F90 parameters_C data
      2,   // ios(92) F90 parameters_C data
      2,   // ios(93) F90 parameters_C data
      2,   // ios(94) F90 parameters_C data
      2,   // ios(95) F90 parameters_C data
      2,   // ios(96) F90 parameters_C data
      2,   // ios(97) F90 parameters_C data
      1,   // ios(98) F90 parameters_C data
      0,   // ios(99) F90 parameters_C data
      3,   // ios(100) F90 parameters_C data
     -3,   // ios(101) F90 parameters_C data
      1,   // ios(102) F90 parameters_C data
      2,   // ios(103) F90 parameters_C data
      1,   // ios(104) F90 parameters_C data
     -2,   // ios(105) F90 parameters_C data
     -1,   // ios(106) F90 parameters_C data
      0,  // ios(107) F90 parameters_C data
};
int iop[N_PARAM + 1] = { 0,  // [0] unused
      0,   // iop(1) F90 parameters_C data
      0,   // iop(2) F90 parameters_C data
      0,   // iop(3) F90 parameters_C data
      0,   // iop(4) F90 parameters_C data
      1,   // iop(5) F90 parameters_C data
      2,   // iop(6) F90 parameters_C data
      3,   // iop(7) F90 parameters_C data
      4,   // iop(8) F90 parameters_C data
      5,   // iop(9) F90 parameters_C data
      6,   // iop(10) F90 parameters_C data
      0,   // iop(11) F90 parameters_C data
      0,   // iop(12) F90 parameters_C data
      1,   // iop(13) F90 parameters_C data
      2,   // iop(14) F90 parameters_C data
      3,   // iop(15) F90 parameters_C data
      4,   // iop(16) F90 parameters_C data
      5,   // iop(17) F90 parameters_C data
      6,   // iop(18) F90 parameters_C data
      0,   // iop(19) F90 parameters_C data
      0,   // iop(20) F90 parameters_C data
      0,   // iop(21) F90 parameters_C data
      0,   // iop(22) F90 parameters_C data
      0,   // iop(23) F90 parameters_C data
      0,   // iop(24) F90 parameters_C data
      0,   // iop(25) F90 parameters_C data
      0,   // iop(26) F90 parameters_C data
      0,   // iop(27) F90 parameters_C data
      0,   // iop(28) F90 parameters_C data
      0,   // iop(29) F90 parameters_C data
      0,   // iop(30) F90 parameters_C data
      1,   // iop(31) F90 parameters_C data
      2,   // iop(32) F90 parameters_C data
      3,   // iop(33) F90 parameters_C data
      4,   // iop(34) F90 parameters_C data
      5,   // iop(35) F90 parameters_C data
      6,   // iop(36) F90 parameters_C data
      0,   // iop(37) F90 parameters_C data
      0,   // iop(38) F90 parameters_C data
      0,   // iop(39) F90 parameters_C data
      0,   // iop(40) F90 parameters_C data
      0,   // iop(41) F90 parameters_C data
      0,   // iop(42) F90 parameters_C data
      0,   // iop(43) F90 parameters_C data
      0,   // iop(44) F90 parameters_C data
      0,   // iop(45) F90 parameters_C data
      0,   // iop(46) F90 parameters_C data
      0,   // iop(47) F90 parameters_C data
      0,   // iop(48) F90 parameters_C data
      1,   // iop(49) F90 parameters_C data
      2,   // iop(50) F90 parameters_C data
      3,   // iop(51) F90 parameters_C data
      4,   // iop(52) F90 parameters_C data
      5,   // iop(53) F90 parameters_C data
      6,   // iop(54) F90 parameters_C data
      0,   // iop(55) F90 parameters_C data
      0,   // iop(56) F90 parameters_C data
      0,   // iop(57) F90 parameters_C data
      0,   // iop(58) F90 parameters_C data
      0,   // iop(59) F90 parameters_C data
      0,   // iop(60) F90 parameters_C data
      0,   // iop(61) F90 parameters_C data
      0,   // iop(62) F90 parameters_C data
      0,   // iop(63) F90 parameters_C data
      0,   // iop(64) F90 parameters_C data
      0,   // iop(65) F90 parameters_C data
      0,   // iop(66) F90 parameters_C data
      0,   // iop(67) F90 parameters_C data
      0,   // iop(68) F90 parameters_C data
      0,   // iop(69) F90 parameters_C data
      0,   // iop(70) F90 parameters_C data
      0,   // iop(71) F90 parameters_C data
      0,   // iop(72) F90 parameters_C data
      0,   // iop(73) F90 parameters_C data
      0,   // iop(74) F90 parameters_C data
      0,   // iop(75) F90 parameters_C data
      0,   // iop(76) F90 parameters_C data
      0,   // iop(77) F90 parameters_C data
      0,   // iop(78) F90 parameters_C data
      0,   // iop(79) F90 parameters_C data
      0,   // iop(80) F90 parameters_C data
      1,   // iop(81) F90 parameters_C data
      2,   // iop(82) F90 parameters_C data
      3,   // iop(83) F90 parameters_C data
      4,   // iop(84) F90 parameters_C data
      5,   // iop(85) F90 parameters_C data
      6,   // iop(86) F90 parameters_C data
      0,   // iop(87) F90 parameters_C data
      0,   // iop(88) F90 parameters_C data
      0,   // iop(89) F90 parameters_C data
      0,   // iop(90) F90 parameters_C data
      0,   // iop(91) F90 parameters_C data
      0,   // iop(92) F90 parameters_C data
      0,   // iop(93) F90 parameters_C data
      0,   // iop(94) F90 parameters_C data
      0,   // iop(95) F90 parameters_C data
      0,   // iop(96) F90 parameters_C data
      0,   // iop(97) F90 parameters_C data
      0,   // iop(98) F90 parameters_C data
      0,   // iop(99) F90 parameters_C data
      0,   // iop(100) F90 parameters_C data
      0,   // iop(101) F90 parameters_C data
      0,   // iop(102) F90 parameters_C data
      0,   // iop(103) F90 parameters_C data
      0,   // iop(104) F90 parameters_C data
      0,   // iop(105) F90 parameters_C data
      0,   // iop(106) F90 parameters_C data
      0,  // iop(107) F90 parameters_C data
};
int iod[N_PARAM + 1] = { 0,  // [0] unused
      0,   // iod(1) F90 parameters_C data
      0,   // iod(2) F90 parameters_C data
      0,   // iod(3) F90 parameters_C data
      0,   // iod(4) F90 parameters_C data
      0,   // iod(5) F90 parameters_C data
      0,   // iod(6) F90 parameters_C data
      0,   // iod(7) F90 parameters_C data
      0,   // iod(8) F90 parameters_C data
      0,   // iod(9) F90 parameters_C data
      0,   // iod(10) F90 parameters_C data
      0,   // iod(11) F90 parameters_C data
      0,   // iod(12) F90 parameters_C data
      0,   // iod(13) F90 parameters_C data
      0,   // iod(14) F90 parameters_C data
      0,   // iod(15) F90 parameters_C data
      0,   // iod(16) F90 parameters_C data
      0,   // iod(17) F90 parameters_C data
      0,   // iod(18) F90 parameters_C data
      0,   // iod(19) F90 parameters_C data
      0,   // iod(20) F90 parameters_C data
      1,   // iod(21) F90 parameters_C data
      2,   // iod(22) F90 parameters_C data
      3,   // iod(23) F90 parameters_C data
      5,   // iod(24) F90 parameters_C data
      5,   // iod(25) F90 parameters_C data
      6,   // iod(26) F90 parameters_C data
      7,   // iod(27) F90 parameters_C data
      8,   // iod(28) F90 parameters_C data
     10,   // iod(29) F90 parameters_C data
      0,   // iod(30) F90 parameters_C data
      0,   // iod(31) F90 parameters_C data
      0,   // iod(32) F90 parameters_C data
      0,   // iod(33) F90 parameters_C data
      0,   // iod(34) F90 parameters_C data
      0,   // iod(35) F90 parameters_C data
      0,   // iod(36) F90 parameters_C data
      0,   // iod(37) F90 parameters_C data
      0,   // iod(38) F90 parameters_C data
      1,   // iod(39) F90 parameters_C data
      2,   // iod(40) F90 parameters_C data
      4,   // iod(41) F90 parameters_C data
      5,   // iod(42) F90 parameters_C data
      5,   // iod(43) F90 parameters_C data
      7,   // iod(44) F90 parameters_C data
      8,   // iod(45) F90 parameters_C data
     10,   // iod(46) F90 parameters_C data
     10,   // iod(47) F90 parameters_C data
      0,   // iod(48) F90 parameters_C data
      0,   // iod(49) F90 parameters_C data
      0,   // iod(50) F90 parameters_C data
      0,   // iod(51) F90 parameters_C data
      0,   // iod(52) F90 parameters_C data
      0,   // iod(53) F90 parameters_C data
      0,   // iod(54) F90 parameters_C data
      0,   // iod(55) F90 parameters_C data
      0,   // iod(56) F90 parameters_C data
      1,   // iod(57) F90 parameters_C data
      0,   // iod(58) F90 parameters_C data
      0,   // iod(59) F90 parameters_C data
      0,   // iod(60) F90 parameters_C data
      0,   // iod(61) F90 parameters_C data
      0,   // iod(62) F90 parameters_C data
      0,   // iod(63) F90 parameters_C data
      0,   // iod(64) F90 parameters_C data
      0,   // iod(65) F90 parameters_C data
      0,   // iod(66) F90 parameters_C data
      0,   // iod(67) F90 parameters_C data
      0,   // iod(68) F90 parameters_C data
      0,   // iod(69) F90 parameters_C data
      0,   // iod(70) F90 parameters_C data
      1,   // iod(71) F90 parameters_C data
      2,   // iod(72) F90 parameters_C data
      3,   // iod(73) F90 parameters_C data
      5,   // iod(74) F90 parameters_C data
      5,   // iod(75) F90 parameters_C data
      6,   // iod(76) F90 parameters_C data
      7,   // iod(77) F90 parameters_C data
      9,   // iod(78) F90 parameters_C data
     10,   // iod(79) F90 parameters_C data
      0,   // iod(80) F90 parameters_C data
      0,   // iod(81) F90 parameters_C data
      0,   // iod(82) F90 parameters_C data
      0,   // iod(83) F90 parameters_C data
      0,   // iod(84) F90 parameters_C data
      0,   // iod(85) F90 parameters_C data
      0,   // iod(86) F90 parameters_C data
      0,   // iod(87) F90 parameters_C data
      0,   // iod(88) F90 parameters_C data
      1,   // iod(89) F90 parameters_C data
      0,   // iod(90) F90 parameters_C data
      0,   // iod(91) F90 parameters_C data
      0,   // iod(92) F90 parameters_C data
      0,   // iod(93) F90 parameters_C data
      0,   // iod(94) F90 parameters_C data
      0,   // iod(95) F90 parameters_C data
      0,   // iod(96) F90 parameters_C data
      0,   // iod(97) F90 parameters_C data
      0,   // iod(98) F90 parameters_C data
      0,   // iod(99) F90 parameters_C data
      0,   // iod(100) F90 parameters_C data
      0,   // iod(101) F90 parameters_C data
      0,   // iod(102) F90 parameters_C data
      0,   // iod(103) F90 parameters_C data
      0,   // iod(104) F90 parameters_C data
      0,   // iod(105) F90 parameters_C data
      0,   // iod(106) F90 parameters_C data
      0,  // iod(107) F90 parameters_C data
};
double gpp[N_PARAM + 1] = {};
double gp2[N_PARAM + 1] = {};
double hsp[N_PARAM + 1] = {};
double gss[N_PARAM + 1] = {};
double gsp[N_PARAM + 1] = {};
double zs[N_PARAM + 1] = {};
double zp[N_PARAM + 1] = {};
double zd[N_PARAM + 1] = {};
int npq[N_PARAM + 1][4] = {};
namespace {
struct NpqInit {
    NpqInit() {
        for (int i = 1; i <= N_PARAM; ++i)
            for (int j = 1; j <= 3; ++j) npq[i][j] = 0;
        static const int ns[108] = {
            0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3,
            3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 3, 0, 0, 0, 0, 0,
        };
        static const int np[108] = {
            0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3,
            3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        };
        static const int nd[108] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3,
            3, 3, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6,
            6, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        };
        for (int i = 1; i <= N_PARAM; ++i) {
            npq[i][1] = ns[i];
            npq[i][2] = np[i];
            npq[i][3] = nd[i];
        }
    }
} _npq_init;
}  // namespace
double uss[N_PARAM + 1] = {};
double upp[N_PARAM + 1] = {};
double udd[N_PARAM + 1] = {};
double polvol[N_PARAM + 1] = {};
double pocord[N_PARAM + 1] = {};
double eisol[N_PARAM + 1] = {};
double eheat[N_PARAM + 1] = { 0.0,  // [0] unused
    52.102   ,  // eheat(1) F90 parameters_C data
    0.000    ,  // eheat(2) F90 parameters_C data
    38.410   ,  // eheat(3) F90 parameters_C data
    76.960   ,  // eheat(4) F90 parameters_C data
    135.700  ,  // eheat(5) F90 parameters_C data
    170.890  ,  // eheat(6) F90 parameters_C data
    113.000  ,  // eheat(7) F90 parameters_C data
    59.559   ,  // eheat(8) F90 parameters_C data
    18.890   ,  // eheat(9) F90 parameters_C data
    0.000    ,  // eheat(10) F90 parameters_C data
    25.650   ,  // eheat(11) F90 parameters_C data
    35.000   ,  // eheat(12) F90 parameters_C data
    79.490   ,  // eheat(13) F90 parameters_C data
    108.390  ,  // eheat(14) F90 parameters_C data
    75.570   ,  // eheat(15) F90 parameters_C data
    66.400   ,  // eheat(16) F90 parameters_C data
    28.990   ,  // eheat(17) F90 parameters_C data
    0.000    ,  // eheat(18) F90 parameters_C data
    21.420   ,  // eheat(19) F90 parameters_C data
    42.600   ,  // eheat(20) F90 parameters_C data
    90.300   ,  // eheat(21) F90 parameters_C data
    112.300  ,  // eheat(22) F90 parameters_C data
    122.900  ,  // eheat(23) F90 parameters_C data
    95.000   ,  // eheat(24) F90 parameters_C data
    67.700   ,  // eheat(25) F90 parameters_C data
    99.300   ,  // eheat(26) F90 parameters_C data
    102.400  ,  // eheat(27) F90 parameters_C data
    102.800  ,  // eheat(28) F90 parameters_C data
    80.700   ,  // eheat(29) F90 parameters_C data
    31.170   ,  // eheat(30) F90 parameters_C data
    65.400   ,  // eheat(31) F90 parameters_C data
    89.500   ,  // eheat(32) F90 parameters_C data
    72.300   ,  // eheat(33) F90 parameters_C data
    54.300   ,  // eheat(34) F90 parameters_C data
    26.740   ,  // eheat(35) F90 parameters_C data
    0.000    ,  // eheat(36) F90 parameters_C data
    19.600   ,  // eheat(37) F90 parameters_C data
    39.100   ,  // eheat(38) F90 parameters_C data
    101.500  ,  // eheat(39) F90 parameters_C data
    145.500  ,  // eheat(40) F90 parameters_C data
    172.400  ,  // eheat(41) F90 parameters_C data
    157.300  ,  // eheat(42) F90 parameters_C data
    162.000  ,  // eheat(43) F90 parameters_C data
    155.500  ,  // eheat(44) F90 parameters_C data
    133.000  ,  // eheat(45) F90 parameters_C data
    90.000   ,  // eheat(46) F90 parameters_C data
    68.100   ,  // eheat(47) F90 parameters_C data
    26.720   ,  // eheat(48) F90 parameters_C data
    58.000   ,  // eheat(49) F90 parameters_C data
    72.200   ,  // eheat(50) F90 parameters_C data
    63.200   ,  // eheat(51) F90 parameters_C data
    47.000   ,  // eheat(52) F90 parameters_C data
    25.517   ,  // eheat(53) F90 parameters_C data
    0.000    ,  // eheat(54) F90 parameters_C data
    18.700   ,  // eheat(55) F90 parameters_C data
    42.500   ,  // eheat(56) F90 parameters_C data
    103.100  ,  // eheat(57) F90 parameters_C data
    101.300  ,  // eheat(58) F90 parameters_C data
    0.000    ,  // eheat(59) F90 parameters_C data
    0.000    ,  // eheat(60) F90 parameters_C data
    0.000    ,  // eheat(61) F90 parameters_C data
    49.400   ,  // eheat(62) F90 parameters_C data
    0.000    ,  // eheat(63) F90 parameters_C data
    0.000    ,  // eheat(64) F90 parameters_C data
    0.000    ,  // eheat(65) F90 parameters_C data
    0.000    ,  // eheat(66) F90 parameters_C data
    0.000    ,  // eheat(67) F90 parameters_C data
    75.800   ,  // eheat(68) F90 parameters_C data
    0.000    ,  // eheat(69) F90 parameters_C data
    36.360   ,  // eheat(70) F90 parameters_C data
    102.100  ,  // eheat(71) F90 parameters_C data
    148.000  ,  // eheat(72) F90 parameters_C data
    186.900  ,  // eheat(73) F90 parameters_C data
    203.100  ,  // eheat(74) F90 parameters_C data
    185.000  ,  // eheat(75) F90 parameters_C data
    188.000  ,  // eheat(76) F90 parameters_C data
    160.000  ,  // eheat(77) F90 parameters_C data
    135.200  ,  // eheat(78) F90 parameters_C data
    88.000   ,  // eheat(79) F90 parameters_C data
    14.690   ,  // eheat(80) F90 parameters_C data
    43.550   ,  // eheat(81) F90 parameters_C data
    46.620   ,  // eheat(82) F90 parameters_C data
    50.100   ,  // eheat(83) F90 parameters_C data
    0.000    ,  // eheat(84) F90 parameters_C data
    0.000    ,  // eheat(85) F90 parameters_C data
    0.000    ,  // eheat(86) F90 parameters_C data
    0.000    ,  // eheat(87) F90 parameters_C data
    0.000    ,  // eheat(88) F90 parameters_C data
    0.000    ,  // eheat(89) F90 parameters_C data
    1674.640 ,  // eheat(90) F90 parameters_C data
    0.000    ,  // eheat(91) F90 parameters_C data
    0.000    ,  // eheat(92) F90 parameters_C data
    0.000    ,  // eheat(93) F90 parameters_C data
    0.000    ,  // eheat(94) F90 parameters_C data
    0.000    ,  // eheat(95) F90 parameters_C data
    0.000    ,  // eheat(96) F90 parameters_C data
    0.000    ,  // eheat(97) F90 parameters_C data
    0.000    ,  // eheat(98) F90 parameters_C data
    0.000    ,  // eheat(99) F90 parameters_C data
    0.000    ,  // eheat(100) F90 parameters_C data
    0.000    ,  // eheat(101) F90 parameters_C data
    207.000  ,  // eheat(102) F90 parameters_C data
    0.000    ,  // eheat(103) F90 parameters_C data
    0.000    ,  // eheat(104) F90 parameters_C data
    0.000    ,  // eheat(105) F90 parameters_C data
    0.000    ,  // eheat(106) F90 parameters_C data
    0.000      // eheat(107) F90 parameters_C data
};
double eheat_sparkles[N_PARAM + 1] = { 0.0,  // [0] unused
    0.0000   ,  // eheat_sparkles(1)
    0.0000   ,  // eheat_sparkles(2)
    0.0000   ,  // eheat_sparkles(3)
    0.0000   ,  // eheat_sparkles(4)
    0.0000   ,  // eheat_sparkles(5)
    0.0000   ,  // eheat_sparkles(6)
    0.0000   ,  // eheat_sparkles(7)
    0.0000   ,  // eheat_sparkles(8)
    0.0000   ,  // eheat_sparkles(9)
    0.0000   ,  // eheat_sparkles(10)
    0.0000   ,  // eheat_sparkles(11)
    0.0000   ,  // eheat_sparkles(12)
    0.0000   ,  // eheat_sparkles(13)
    0.0000   ,  // eheat_sparkles(14)
    0.0000   ,  // eheat_sparkles(15)
    0.0000   ,  // eheat_sparkles(16)
    0.0000   ,  // eheat_sparkles(17)
    0.0000   ,  // eheat_sparkles(18)
    0.0000   ,  // eheat_sparkles(19)
    0.0000   ,  // eheat_sparkles(20)
    0.0000   ,  // eheat_sparkles(21)
    0.0000   ,  // eheat_sparkles(22)
    0.0000   ,  // eheat_sparkles(23)
    0.0000   ,  // eheat_sparkles(24)
    0.0000   ,  // eheat_sparkles(25)
    0.0000   ,  // eheat_sparkles(26)
    0.0000   ,  // eheat_sparkles(27)
    0.0000   ,  // eheat_sparkles(28)
    0.0000   ,  // eheat_sparkles(29)
    0.0000   ,  // eheat_sparkles(30)
    0.0000   ,  // eheat_sparkles(31)
    0.0000   ,  // eheat_sparkles(32)
    0.0000   ,  // eheat_sparkles(33)
    0.0000   ,  // eheat_sparkles(34)
    0.0000   ,  // eheat_sparkles(35)
    0.0000   ,  // eheat_sparkles(36)
    0.0000   ,  // eheat_sparkles(37)
    0.0000   ,  // eheat_sparkles(38)
    0.0000   ,  // eheat_sparkles(39)
    0.0000   ,  // eheat_sparkles(40)
    0.0000   ,  // eheat_sparkles(41)
    0.0000   ,  // eheat_sparkles(42)
    0.0000   ,  // eheat_sparkles(43)
    0.0000   ,  // eheat_sparkles(44)
    0.0000   ,  // eheat_sparkles(45)
    0.0000   ,  // eheat_sparkles(46)
    0.0000   ,  // eheat_sparkles(47)
    0.0000   ,  // eheat_sparkles(48)
    0.0000   ,  // eheat_sparkles(49)
    0.0000   ,  // eheat_sparkles(50)
    0.0000   ,  // eheat_sparkles(51)
    0.0000   ,  // eheat_sparkles(52)
    0.0000   ,  // eheat_sparkles(53)
    0.0000   ,  // eheat_sparkles(54)
    0.0000   ,  // eheat_sparkles(55)
    0.0000   ,  // eheat_sparkles(56)
    928.9000 ,  // eheat_sparkles(57)
    944.7000 ,  // eheat_sparkles(58)
    952.9000 ,  // eheat_sparkles(59)
    962.8000 ,  // eheat_sparkles(60)
    976.9000 ,  // eheat_sparkles(61)
    974.4000 ,  // eheat_sparkles(62)
    1006.6000,  // eheat_sparkles(63)
    991.3700 ,  // eheat_sparkles(64)
    999.0000 ,  // eheat_sparkles(65)
    1001.3000,  // eheat_sparkles(66)
    1009.6000,  // eheat_sparkles(67)
    1016.1500,  // eheat_sparkles(68)
    1022.0600,  // eheat_sparkles(69)
    1039.0300,  // eheat_sparkles(70)
    1031.2000,  // eheat_sparkles(71)
    0.0000   ,  // eheat_sparkles(72)
    0.0000   ,  // eheat_sparkles(73)
    0.0000   ,  // eheat_sparkles(74)
    0.0000   ,  // eheat_sparkles(75)
    0.0000   ,  // eheat_sparkles(76)
    0.0000   ,  // eheat_sparkles(77)
    0.0000   ,  // eheat_sparkles(78)
    0.0000   ,  // eheat_sparkles(79)
    0.0000   ,  // eheat_sparkles(80)
    0.0000   ,  // eheat_sparkles(81)
    0.0000   ,  // eheat_sparkles(82)
    0.0000   ,  // eheat_sparkles(83)
    0.0000   ,  // eheat_sparkles(84)
    0.0000   ,  // eheat_sparkles(85)
    0.0000   ,  // eheat_sparkles(86)
    0.0000   ,  // eheat_sparkles(87)
    0.0000   ,  // eheat_sparkles(88)
    0.0000   ,  // eheat_sparkles(89)
    0.0000   ,  // eheat_sparkles(90)
    0.0000   ,  // eheat_sparkles(91)
    0.0000   ,  // eheat_sparkles(92)
    0.0000   ,  // eheat_sparkles(93)
    0.0000   ,  // eheat_sparkles(94)
    0.0000   ,  // eheat_sparkles(95)
    0.0000   ,  // eheat_sparkles(96)
    0.0000   ,  // eheat_sparkles(97)
    0.0000   ,  // eheat_sparkles(98)
    0.0000   ,  // eheat_sparkles(99)
    0.0000   ,  // eheat_sparkles(100)
    0.0000   ,  // eheat_sparkles(101)
    0.0000   ,  // eheat_sparkles(102)
    0.0000   ,  // eheat_sparkles(103)
    0.0000   ,  // eheat_sparkles(104)
    0.0000   ,  // eheat_sparkles(105)
    0.0000   ,  // eheat_sparkles(106)
    0.0000     // eheat_sparkles(107)
};
double f0sd[N_PARAM + 1] = {};
double g2sd[N_PARAM + 1] = {};
double f0sd_store[N_PARAM + 1] = {};
double g2sd_store[N_PARAM + 1] = {};
struct MGInit {
    MGInit() {
        // Fortran data main_group: T for H..Ca(1-20), Zn..Kr(30-36),
        // Cd..Xe(48-54), Hg..Tv(80-107); F elsewhere.
        for (int i = 1; i <= N_PARAM; ++i) main_group[i] = false;
        for (int i = 1; i <= 20; ++i) main_group[i] = true;
        for (int i = 30; i <= 36; ++i) main_group[i] = true;
        for (int i = 48; i <= 54; ++i) main_group[i] = true;
        for (int i = 80; i <= N_PARAM; ++i) main_group[i] = true;
    }
} _mg_init;
bool main_group[N_PARAM + 1] = {};
double zsn[N_PARAM + 1] = {};
double zpn[N_PARAM + 1] = {};
double zdn[N_PARAM + 1] = {};
bool dorbs[N_PARAM + 1] = {};
// PM7 par1..par4 = gues7*(99,*) from parameters_for_PM7_C.F90 (PAR1..PAR4)
double par1 = 8.947612, par2 = 6.024265, par3 = -0.012037, par4 = 0.701333;
int ndelec[N_PARAM + 1] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 4, 4, 6,
    8, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
    2, 4, 4, 6, 8, 10, 10, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 2, 2, 4, 4, 6, 8, 10, 10, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

double ddp[7][108];
double gss6sp[N_PARAM + 1] = {};
double gssam1sp[N_PARAM + 1] = {};
double gssPM3sp[N_PARAM + 1] = {};
double dh2_a_parameters[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double f0dd[N_PARAM + 1] = {};
double f2dd[N_PARAM + 1] = {};
double f4dd[N_PARAM + 1] = {};
double f0pd[N_PARAM + 1] = {};
double f2pd[N_PARAM + 1] = {};
double g1pd[N_PARAM + 1] = {};
double g3pd[N_PARAM + 1] = {};

}

// REF (PARAM) parameter metadata, parameters_C.F90:301-337 (used by PARAM/REF only).
const char* partyp[37 + 1] = {
    nullptr,
    "USS",
    "UPP",
    "UDD",
    "ZS",
    "ZP",
    "ZD",
    "BETAS",
    "BETAP",
    "BETAD",
    "GSS",
    "GSP",
    "GPP",
    "GP2",
    "HSP",
    "F0SD",
    "G2SD",
    "POC",
    "ALP",
    "ZSN",
    "ZPN",
    "ZDN",
    "FN11",
    "FN21",
    "FN31",
    "FN12",
    "FN22",
    "FN32",
    "FN13",
    "FN23",
    "FN33",
    "FN14",
    "FN24",
    "FN34",
    "ALPB_",
    "XFAC_",
    "FN33",
    "NORBS",
};
double defmin[37 + 1] = {
    0.0,
    -200,
    -200,
    -200,
    0.6,
    0.6,
    0.6,
    -70,
    -70,
    -70,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    0.5,
    0.5,
    0.5,
    0.5,
    -1,
    0.5,
    0.5,
    -1,
    0.15,
    0.5,
    -1,
    1,
    0.5,
    -1,
    1,
    0.5,
    0.9,
    0.5,
    0.5,
    0.9,
};
double defmax[37 + 1] = {
    0.0,
    60,
    60,
    60,
    6,
    6,
    6,
    10,
    10,
    10,
    20,
    20,
    20,
    19,
    5,
    10,
    10,
    10,
    6,
    25,
    25,
    25,
    1,
    3,
    3,
    1,
    3,
    3.5,
    1,
    3,
    4,
    1,
    3,
    3,
    3,
    30,
    6,
    9.1,
};

