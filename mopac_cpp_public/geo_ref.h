// geo_ref.h — C++ translation of geo_ref.F90 (MOPAC 2016).
#pragma once

// geo_ref reads in the geometry reference file and makes it as similar as
// possible to the current geometry: PDB labels are re-sequenced, the
// geometries are rotated/translated into maximum coincidence (RMS minimum),
// and pairs of atoms are swapped when that lowers the RMS difference.
void geo_ref();
