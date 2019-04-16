#ifndef _IOLIBISI_
#define _IOLIBISI_
#include "isinglib.h"
#include <string>

using namespace std;

void WriteLatticeToBinary(string fname, lattice &alattice);
lattice ReadLattticeFromBinary(string fname);
#endif

