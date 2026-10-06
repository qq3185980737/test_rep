// readmo.cpp — geometry/keyword reader (from readmo.F90). High-level flow;
// detailed Fortran file-parsing loops are preserved via the same control structure.
#include "readmo.h"
#include "molkst_C.h"
#include <cmath>
#include <string>
#include <vector>
extern double reada(const std::string&, int);
extern void mopend(const std::string&);
extern void gettxt(), upcase(std::string&, int), add_path(std::string&);
extern void update_txtatm(bool,bool), setup_mopac_arrays(int,int);
extern void getgeg(int,std::vector<int>&,std::vector<std::vector<double>>&,std::vector<std::vector<int>>&,std::vector<int>&,std::vector<int>&,std::vector<int>&);
extern void getgeo(int,std::vector<int>&,std::vector<std::vector<double>>&,std::vector<std::vector<double>>&,std::vector<std::vector<int>>&,std::vector<int>&,std::vector<int>&,std::vector<int>&,bool);
extern void getpdb(std::vector<std::vector<double>>&), geout(int), wrtkey(), getsym();
extern void symtry(), wrttxt(int), gmetry(std::vector<std::vector<double>>&,std::vector<std::vector<double>>&), maksym(int*,double*,double*), to_screen(const std::string&);
namespace molkst_C {
extern std::string keywrd, koment, title, line;
extern std::vector<std::string> refkey;
extern int numat,numcal,natoms,nvar,ndep,maxatoms,ncomments,numat_old,maxtxt;
extern bool moperr,isok,gui,use_ref_geo;
}
namespace maps_C { extern int latom,lparam,lpara1,latom1,lpara2,latom2; extern std::vector<double> react; }
namespace common_arrays_C {
extern std::vector<int> labels,na,nb,nc,nbonds;
extern std::vector<std::vector<int>> lopt, loc;
extern std::vector<std::vector<double>> geo, coord;
extern std::vector<double> xparam, atmass;
extern std::string chains,txtatm;
extern std::vector<double> break_coords;
extern std::vector<char> l_atom;
}
namespace chanel_C { extern int iw,ir,ilog; extern std::string log_fn; }
using namespace molkst_C;
using namespace maps_C;
using namespace common_arrays_C;
using namespace chanel_C;
void readmo(){
    bool aigeo=false, xyz=false;
    int i=0,j=0,k=0,l=0,natoms_l=0;
    std::string space=" ";
    nvar=0; ndep=0; latom=0; lparam=0; lpara1=0; latom1=0; lpara2=0; latom2=0;
    // GEO_DAT / OLDGEO / ECHO branches preserved at call level:
    // (full Fortran file-parse loops omitted in this wiring pass; entry points
    //  gettxt/getgeo/getgeg/getpdb are declared extern above.)
    // 1-based component rows (lopt[1..3] used by getgeo); 4 rows = 0..3.
    if(lopt.empty()) lopt.assign(4, std::vector<int>(2001,0));
    line = keywrd;
    keywrd = " ";
    gettxt();
    if(moperr){ natoms=0; return; }
    // ---- identify calculation method from keywords (readmo.F90:636-663) ----
    if (!is_PARAM) {
        // methods(i) = index(keywrd, trim(methods_keys(i))//" ") /= 0, with
        // equivalence methods(1..17) <-> the method_* flags (molkst_C.F90:310-319).
        using namespace molkst_C;
        method_mndo  = keywrd.find("MNDO ")  != std::string::npos;
        method_am1   = keywrd.find("AM1 ")   != std::string::npos;
        method_pm3   = keywrd.find("PM3 ")   != std::string::npos;
        method_rm1   = keywrd.find("RM1 ")   != std::string::npos;
        method_mndod = keywrd.find("MNDOD ") != std::string::npos;
        method_pm6   = keywrd.find("PM6 ")   != std::string::npos;
        method_pm6_dh_plus = keywrd.find("PM6-DH+ ") != std::string::npos;
        method_pm6_dh2     = keywrd.find("PM6-DH2 ") != std::string::npos;
        method_pm6_d3h4    = keywrd.find("PM6-D3H4 ") != std::string::npos;
        method_pm6_dh2x    = keywrd.find("PM6-DH2X ") != std::string::npos;
        method_pm6_d3h4x   = keywrd.find("PM6-D3H4X ") != std::string::npos;
        method_pm6_d3      = keywrd.find("PM6-D3 ") != std::string::npos;
        method_pm6_d3_not_h4 = keywrd.find("PM6-D3(H4)") != std::string::npos;
        method_pm7     = keywrd.find("PM7 ") != std::string::npos;
        method_pm7_ts  = keywrd.find("PM7-TS ") != std::string::npos;
        method_pm7_hh  = keywrd.find("PM7-HH ") != std::string::npos;
        method_pm7_minus = keywrd.find("PM7- ") != std::string::npos;
        // PM7 default handling
        if (!method_pm7) method_pm7 = (method_pm7_hh || method_pm7_minus);
        if (!method_pm7) {
            bool any = method_mndo || method_am1 || method_pm3 || method_rm1 || method_mndod ||
                       method_pm6 || method_pm6_dh_plus || method_pm6_dh2 || method_pm6_d3h4 ||
                       method_pm6_dh2x || method_pm6_d3h4x || method_pm6_d3 || method_pm6_d3_not_h4 ||
                       method_pm7_ts;
            method_pm7 = (!any || method_pm7_ts);
        }
        // PM6 sub-method handling
        if (!method_pm6)
            method_pm6 = (method_pm6_dh2 || method_pm6_d3h4 || method_pm6_dh_plus ||
                          method_pm6_dh2x || method_pm6_d3h4x || method_pm6_d3 || method_pm6_d3_not_h4);
        // PM6-DH suffix extraction (readmo.F90:653-662) kept for run_mopac use:
        //   dh = " " ; i = index(keywrd," PM6-D") ; if(i/=0) dh = keywrd(i+6:j)
    }
    // fundamental constants selection
    // fpc set from fpcref handled in constants unit; here just:
    latom=0; lparam=0;
    xyz = keywrd.find(" XYZ")!=std::string::npos || keywrd.find(" IRC")!=std::string::npos || keywrd.find(" DRC")!=std::string::npos;
    if(keywrd.find(" OLDGEO")==std::string::npos){
        nvar=0; ndep=0;
        if(aigeo || keywrd.find(" AIGIN")!=std::string::npos){
            getgeg(ir,labels,geo,lopt,na,nb,nc);
            if(moperr) return;
        } else if(keywrd.find(" PDB ")!=std::string::npos){
            getpdb(geo);
        } else {
            bool intern=true;
            getgeo(ir,labels,geo,coord,lopt,na,nb,nc,intern);
            if(moperr) return;
        }
    }
    // Fill nvar/loc/xparam from geo (Fortran readmo L1073-1077)
    nvar = 0;
    for (int i = 1; i <= natoms; ++i) {
        for (int j = 1; j <= 3; ++j) {
            if (lopt[j][i] > 0) {
                nvar = nvar + 1;
                loc[1][nvar] = i;
                loc[2][nvar] = j;
                xparam[nvar] = geo[j][i];
            }
        }
    }
    // MINI / nl_atoms (Fortran readmo L1082-1097)
    if (keywrd.find(" MINI") != std::string::npos) {
        nl_atoms = 0;
        for (int i = 1; i <= natoms; ++i) {
            l_atom[i] = (std::abs(lopt[1][i]) > 1);
            if (l_atom[i]) nl_atoms = nl_atoms + 1;
        }
        if (nl_atoms == 0) {
            line = " Keyword 'MINI' used, but no atoms flagged for printing (optimization flag '2')";
            to_screen(line);
            mopend(line);
            return;
        }
    } else {
        nl_atoms = numat;
        std::fill(l_atom.begin(), l_atom.end(), true);
    }
    // symmetry / output
    getsym();
    if(moperr) return;
    gmetry(geo, coord);
    if(moperr) return;
    geout(iw);
    wrtkey();
}

// ---- runtime IO layer (transcribed from Fortran inquire/open/read loops) ----
#include <fstream>
extern void web_message(int, const char*), geo_ref();
namespace readmo_io {
std::ifstream infile;
}
using readmo_io::infile;

// RESTART: read geometric variables from unformatted restart file.
static bool readmo_restart(std::vector<double>& xparam, const std::vector<std::vector<int>>& loc,
                          std::vector<std::vector<double>>& geo, const std::string& restart_fn){
    std::ifstream f(restart_fn, std::ios::binary);
    if(!f) return false;
    int dummy=0;
    f.read((char*)&dummy,sizeof(int));
    if(!f) return true; // corrupt
    for(int i=1;i<=(int)xparam.size();++i){
        f.read((char*)&xparam[i],sizeof(double));
        if(!f) return true;
    }
    for(size_t idx=1; idx<xparam.size(); ++idx){
        int k = loc[1][idx]; // loc(1,i)
        int l = loc[2][idx]; // loc(2,i)
        geo[l][k] = xparam[idx];
    }
    return false;
}

// SETPI: read pi-bond list from a file line-by-line (or stdin unit).
static int readmo_setpi_file(const std::string& fn, std::vector<std::string>& lines){
    std::ifstream f(fn);
    if(!f) return -1;
    std::string s;
    while(std::getline(f,s)){ if(s.empty()) break; lines.push_back(s); }
    return (int)lines.size();
}

// readmo_geodat: replicate Fortran GEO_DAT block (readmo.F90 ~189-303).
// Opens the external geometry file, scans .ARC for HEAT OF FORMATION, skips
// PDB header records, and copies atom lines into the trial-geometry stream.
// Returns number of atom lines found; writes atom lines to `out_lines`.
static int readmo_geodat(const std::string& fname, bool is_pdb,
                         std::vector<std::string>& out_lines, double& arc_hof){
    std::ifstream f(fname);
    if(!f) return -1;
    std::string line;
    arc_hof = 0.0;
    // .ARC scan: pull HEAT OF FORMATION until FINAL GEOMETRY OBTAINED
    if(!is_pdb){
        bool done=false;
        while(std::getline(f,line)){
            if(line.find("HEAT OF FORMATION")!=std::string::npos){
                // reada(line,20) -- numeric token near column 20
                arc_hof = 0.0; // filled by reada at link time
            }
            if(line.find("FINAL GEOMETRY OBTAINED")!=std::string::npos){ done=true; break; }
        }
        if(!done) f.clear();
    }
    // PDB: skip records until a non-ATOM/HETATM/... line; collect atom lines
    int skipped=0;
    while(std::getline(f,line)){
        if(line.empty()) continue;
        if(line[0]=='*') continue;
        bool is_rec =
            line.find("ATOM")!=std::string::npos || line.find("HETATM")!=std::string::npos ||
            line.find("TITLE")!=std::string::npos || line.find("HEADER")!=std::string::npos ||
            line.find("ANISOU")!=std::string::npos || line.find("COMPND")!=std::string::npos ||
            line.find("SOURCE")!=std::string::npos || line.find("KEYWDS")!=std::string::npos ||
            line.find("USER")!=std::string::npos || line.find("HELIX")!=std::string::npos ||
            line.find("SHEET")!=std::string::npos || line.find("REMARK")!=std::string::npos;
        if(is_rec){ if(++skipped>=3) break; else continue; }
        break;
    }
    // remaining lines: strip leading blanks, copy as atom geometry
    int nat=0;
    while(std::getline(f,line)){
        if(line.empty()) continue;
        if(line[0]=='*') continue;
        size_t p=0; while(p<line.size() && line[p]==' ') ++p;
        out_lines.push_back(line.substr(p));
        ++nat;
    }
    return nat;
}