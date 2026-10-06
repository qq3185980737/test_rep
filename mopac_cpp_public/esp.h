// esp.h — C++ interface for MOPAC 2016 "esp.F90" (ESP charge fitting).
#pragma once

void esp();
void pdgrid();
void surfac();
void potcal();
void elesn();
void espfit();
void naicas();
void ovlp(int ic);
void setup3();
void fsub(int n, double x, double& fval);
void naicap();
void setup_esp(int mode);
void getattrib(double xmin[4], double xmax[4]);
void espplane(int iplane, const double xmin[4], const double step[4],
              int ngridpts2, int ngridpts3);
double dex2(int m);
