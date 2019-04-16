#ifndef _IOLIBISI_
#define _IOLIBISI_
#include "latticelib.h"
#include <string>

using namespace std;

void WriteLatticeToBinary(string fname, lattice &alattice);
lattice ReadLattticeFromBinary(string fname);

void WriteLatticeToBinary_xyModel(string fname, lattice alattice,bool lattice_skip=false);
lattice ReadLattticeFromBinary_xyModel(string fname);
#endif

