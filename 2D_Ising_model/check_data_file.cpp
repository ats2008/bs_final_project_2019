#include <iostream>
#include <ctime>
#include "isinglib.h"
#include "iolib.h"
#include "stattools.h"
#include <string>
#include <stdlib.h>

using namespace std;

int main()
{
	string file_name="multi_histo_data/N10/N_10_T_2.400_ising.dat";
	lattice alat=ReadLattticeFromBinary(file_name);
	alat.display_lattice_detailed();	
	alat.clear_lattice();
	file_name="multi_histo_data/N10/N_10_T_2.400_ising.datR";
	alat=ReadLattticeFromBinary(file_name);
	alat.display_lattice_detailed();	
	alat.clear_lattice();
	return 0;
}
