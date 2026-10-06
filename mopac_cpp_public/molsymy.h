// molsymy.h — C++ translation.
#pragma once
void molsym(double* coord, int& ierror, double* r);
void chi(double toler, double* coord, int ioper, int& iqual);
void makopr(int numat, double* coord, int& ierror, double* r);
void orient(int numat, double* coord, double* r);
void plato(double* coord, double* r, bool& cubic);
bool symdec(int n1, const int* ielem);
void cartab();
