// init_filenames.cpp
#include "init_filenames.h"
#include <string>
#include <algorithm>
#include "molkst_C.h"
#include "chanel_C.h"

using namespace molkst_C;
using namespace chanel_C;

extern "C" void upcase_(char*, int*);

void init_filenames() {
    int text_length = (int)jobnam.size();
    std::string line = jobnam;
    if (text_length > 3) {
        if (jobnam[text_length - 4] == '.') {
            std::string ext = line.substr(text_length - 4, 4);
            std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);
            if (ext == ".MOP" || ext == ".DAT" || ext == ".ARC" ||
                ext == ".PDB" || ext == ".ENT" || ext == ".NEW")
                text_length -= 4;
        }
    }
    std::string base = jobnam.substr(0, text_length);
    input_fn    = base + ".temp";
    output_fn   = base + ".out";
    restart_fn  = base + ".res";
    density_fn  = base + ".den";
    log_fn      = base + ".log";
    end_fn      = base + ".end";
    archive_fn  = base + ".arc";
    brillouin_fn = base + ".brz";
    esp_fn      = base + ".esp";
    ump_fn      = base + ".ump";
    mep_fn      = base + ".mep";
    pol_fn      = base + ".pol";
    gpt_fn      = base + ".gpt";
    esr_fn      = base + ".esr";
    xyz_fn      = base + ".xyz";
    syb_fn      = base + ".syb";
    cosmo_fn    = base + ".cos";
    if (gui) output_fn = "OUTPUT file";
}

extern "C" void upcase_(char*, int*) {}
