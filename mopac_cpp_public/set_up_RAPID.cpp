// set_up_RAPID.cpp — translation of MOPAC 2016 set_up_RAPID.F90.
#include "set_up_RAPID.h"

#include "common_arrays_C.h"
#include "compfg.h"
#include "hcore_for_MOZYME.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "picopt.h"
#include "pinout.h"

using namespace common_arrays_C;
using namespace molkst_C;
using namespace MOZYME_C;

// RAPID mode save/restore of MOZYME bookkeeping (mode = "ON"/"OFF"/"RES").
// Stub-level adaptation: rapid SCF switches handled at the F90 driver level;
// the entry keeps MOPAC's state transitions faithful.
void set_up_rapid(const char* txt) {
    static int store_numred=0, store_nelred=0, store_mode=0;
    static bool store_use_ref_geo=false;
    if (txt[0]=='O'&&txt[1]=='F') {
        store_numred=numred; store_mode=mode; store_nelred=nelred;
        numred=numat; mode=0; nelred=nelecs;
        return;
    }
    if (txt[0]=='R'&&txt[1]=='E') {
        numred=store_numred; mode=store_mode; nelred=store_nelred;
        return;
    }
    store_use_ref_geo=use_ref_geo; use_ref_geo=false;
    numred=numat; mode=0; nelred=nelecs;
    picopt(-1);
    mode=0;
    compfg(xparam, true, escf, true, grad, true);
    pinout(1);
    picopt(1);
    iflepo=1;
    step_num=step_num+1;
    mode=-1;
    hcore_for_MOZYME();
    if (moperr) return;
    mode=1;
    use_ref_geo=store_use_ref_geo;
}
