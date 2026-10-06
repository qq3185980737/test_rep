// addhb.h — C++ translation of MOPAC 2016 "addhb.F90"
//
// CHECK FOCK MATRIX AND DENSITY MATRIX TO ENSURE THAT ALL FOCK MATRIX
// ELEMENTS ARE BEING CORRECTLY HANDLED BY THE DIAGONALISER.
//
// nocc1  Number of occupied M.O.s (might be a subset of noccupied)
// nvir1  Number of virtual M.O.s   (might be a subset of nvirtual)
// nij    Out: number of hydrogen bonds that need to be created (0 if none)
// nhb    Threshold selector: 1..4, selects hblims(nhb).
void addhb(int nocc1, int nvir1, int idiagg, int& nij, int nhb);
