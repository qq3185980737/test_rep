// polar_helpers.h - independent helpers of polar.F90 (MOPAC2016 -> C++)
// Fortran 1-based indexing retained (vectors padded by one element).
#ifndef MOPAC_POLAR_HELPERS_H
#define MOPAC_POLAR_HELPERS_H

#include <string>
#include <vector>

void zerom(std::vector<std::vector<double>>& x, int m);
void tf(std::vector<std::vector<double>>& ua, std::vector<std::vector<double>>& ga,
        std::vector<std::vector<double>>& ub, std::vector<std::vector<double>>& gb,
        std::vector<std::vector<double>>& t, int norbs);
void transf(std::vector<std::vector<double>>& f, std::vector<std::vector<double>>& g,
            std::vector<std::vector<double>>& c, int norb);
double trsub(std::vector<std::vector<double>>& ul, std::vector<std::vector<double>>& x,
             std::vector<std::vector<double>>& ur, int l1, int lm, int ndim);
double trudgu(std::vector<std::vector<double>>& ul, std::vector<std::vector<double>>& x,
              std::vector<std::vector<double>>& ur, int l1, int lm, int ndim);
double trugdu(std::vector<std::vector<double>>& ul, std::vector<std::vector<double>>& x,
              std::vector<std::vector<double>>& ur, int l1, int lm, int ndim);
double trugud(std::vector<std::vector<double>>& ul, std::vector<std::vector<double>>& x,
              std::vector<std::vector<double>>& ur, int l1, int lm, int ndim);
double wrdkey(const std::string& keywrd_, const std::string& key, int nk,
              const std::string& refkey, int nr, double def);
double aval(std::vector<std::vector<double>>& h, std::vector<std::vector<double>>& d, int norbs);
double pol_vol(double average);

void copym(std::vector<std::vector<double>>& h, std::vector<std::vector<double>>& f, int m);
void hplusf(std::vector<std::vector<double>>& f, std::vector<std::vector<double>>& h, int norbs);
void hmuf(std::vector<std::vector<double>>& h1, int id, std::vector<std::vector<double>>& coord_,
          std::vector<int>& nfirst_, std::vector<int>& nlast_, std::vector<int>& nat_,
          int norbs, int numat);
void makeuf(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& uold,
            std::vector<std::vector<double>>& g, std::vector<double>& eigs_, bool& last,
            int norbs, int nclose, double& diff, double atol);
void densf(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& c,
           std::vector<std::vector<double>>& d, std::vector<std::vector<double>>& da,
           int norbs, int nclose, std::vector<double>& w2);
void ffreq2(std::vector<std::vector<double>>& f, std::vector<std::vector<double>>& ptot,
            std::vector<double>& w);
void ffreq1(std::vector<std::vector<double>>& f, std::vector<std::vector<double>>& ptot,
            std::vector<double>& pa, std::vector<double>& pb, int norbs);
void dawrt1(std::vector<double>& v, int len, int idaf_, int ns);
void dawrit(std::vector<std::vector<double>>& v, int len, int nrec);
void darea1(std::vector<double>& v, int len, int idaf_, int ns);
void daread(std::vector<std::vector<double>>& v, int len, int nrec);
void openda(int irest);
void betal1(std::vector<std::vector<double>>& u0a, std::vector<std::vector<double>>& g0a,
            std::vector<std::vector<double>>& u1b, std::vector<std::vector<double>>& g1b,
            std::vector<std::vector<double>>& u1c, std::vector<std::vector<double>>& g1c,
            int nclose, int norbs, double& term);
void betall(std::vector<std::vector<double>>& u2a, std::vector<std::vector<double>>& g2a,
            std::vector<std::vector<double>>& u1b, std::vector<std::vector<double>>& g1b,
            std::vector<std::vector<double>>& u1c, std::vector<std::vector<double>>& g1c,
            int nclose, int norbs, double& term);
void betcom(std::vector<std::vector<double>>& u1, std::vector<std::vector<double>>& g1,
            std::vector<std::vector<double>>& u2, std::vector<std::vector<double>>& g2,
            int nclose, int norbs, double& term);

#endif
