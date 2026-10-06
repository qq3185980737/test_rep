// den_in_out.cpp — C++ translation of "den_in_out" (iter.F90, contained routine).
// Read/write the density matrix; MOZYME jobs delegate to pinout().
#include "den_in_out.h"
#include "molkst_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include <cstdio>
#include <fstream>
#include <string>

namespace molkst_C {
extern int norbs, numat;
extern bool uhf, mozyme;
extern std::string keywrd;
}
namespace chanel_C {
extern int iden;
extern std::string density_fn;
}
namespace common_arrays_C {
extern std::vector<double> p, pa, pb;
}

extern void pinout(int);
extern void mopend(const std::string&);
extern void to_screen(const std::string&);

namespace {
// Sequential binary read of vector elements 1..size-1.
bool read_vec(std::fstream& fs, std::vector<double>& v) {
    for (size_t i = 1; i < v.size(); ++i) fs.read((char*)&v[i], sizeof(double));
    return !fs.fail();
}
bool write_vec(std::fstream& fs, const std::vector<double>& v) {
    for (size_t i = 1; i < v.size(); ++i) fs.write((const char*)&v[i], sizeof(double));
    return !fs.fail();
}
}  // namespace

void den_in_out(int mode) {
    if (molkst_C::mozyme) {
        pinout(mode);
        return;
    }
    bool formatted = (molkst_C::keywrd.find(" DENOUTF") != std::string::npos);
    int io_stat = 0, old_norbs = 0, old_numat = 0;
    for (int icount = 1; icount <= 2; ++icount) {
        std::fstream fs;
        std::ios::openmode om = std::ios::binary;
        if (mode == 0) om |= std::ios::in;
        else om |= std::ios::in | std::ios::out | std::ios::trunc;
        fs.open(chanel_C::density_fn, om);
        fs.clear();
        fs.seekg(0);
        fs.seekp(0);
        if (mode == 0) {
            if (formatted) {
                fs >> old_norbs >> old_numat;
                double v;
                for (size_t i = 1; i < common_arrays_C::pa.size(); ++i) { fs >> v; common_arrays_C::pa[i] = v; }
            } else {
                fs.read((char*)&old_norbs, sizeof(int));
                fs.read((char*)&old_numat, sizeof(int));
                read_vec(fs, common_arrays_C::pa);
            }
            io_stat = fs.fail() ? 1 : 0;
            if (old_norbs > 0 && old_norbs < 100000) {
                if (molkst_C::norbs != old_norbs || molkst_C::numat != old_numat) {
                    mopend("Density file read in does not match current data set");
                    return;
                }
            }
            if (icount < 2 && io_stat != 0) {
                formatted = !formatted;
                fs.close();
                continue;
            }
            if (io_stat != 0) {
                to_screen(" Density Restart File missing or corrupt");
                mopend("Density Restart File missing or corrupt");
                return;
            }
            if (molkst_C::uhf) {
                if (formatted) {
                    double v;
                    for (size_t i = 1; i < common_arrays_C::pb.size(); ++i) { fs >> v; common_arrays_C::pb[i] = v; }
                } else {
                    read_vec(fs, common_arrays_C::pb);
                }
                io_stat = fs.fail() ? 1 : 0;
                if (io_stat != 0) {
                    to_screen(" Beta Density Restart File missing or corrupt");
                    to_screen(" (Most likely the previous job did not use UHF)");
                    mopend("Beta Density Restart File missing or corrupt");
                    return;
                }
                for (size_t i = 1; i < common_arrays_C::p.size(); ++i)
                    common_arrays_C::p[i] = common_arrays_C::pa[i] + common_arrays_C::pb[i];
            } else {
                for (size_t i = 1; i < common_arrays_C::p.size(); ++i)
                    common_arrays_C::p[i] = 2.0 * common_arrays_C::pa[i];
            }
        } else {
            if (formatted) {
                fs << molkst_C::norbs << " " << molkst_C::numat;
                for (size_t i = 1; i < common_arrays_C::pa.size(); ++i) fs << " " << common_arrays_C::pa[i];
            } else {
                fs.write((const char*)&molkst_C::norbs, sizeof(int));
                fs.write((const char*)&molkst_C::numat, sizeof(int));
                write_vec(fs, common_arrays_C::pa);
            }
            io_stat = fs.fail() ? 1 : 0;
            if (molkst_C::uhf) {
                if (formatted) for (size_t i = 1; i < common_arrays_C::pb.size(); ++i) fs << " " << common_arrays_C::pb[i];
                else write_vec(fs, common_arrays_C::pb);
            }
            fs.close();
            return;
        }
        fs.close();
    }
}
