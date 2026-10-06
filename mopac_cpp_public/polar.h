// polar.h — C++ translation.
#pragma once
#include <vector>
#include <string>
void polar();
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