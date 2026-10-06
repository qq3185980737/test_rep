// bonds_for_MOZYME.h — C++ translation of MOPAC 2016
// "bonds_for_MOZYME.F90". MOZYME-style bond-order/valency report.
#pragma once

// ijbo: packed AO-block start in p between atoms i,j (external).
int ijbo(int ii, int jj);

void bonds_for_MOZYME();
