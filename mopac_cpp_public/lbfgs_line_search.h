#include <limits>
void dcstep(double& stx,double& fx,double& dx,double& sty,double& fy,double& dy,double& stp,double fp,double dp,bool& brackt,double stpmin,double stpmax);
void dcsrch(double f,double g,double& stp,double ftol,double gtol,double xtol,double stpmin,double stpmax,std::string& task,int* isave,double* dsave);
double dpmeps();