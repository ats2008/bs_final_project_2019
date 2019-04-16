#include <iostream>
#include <ctime>
#include "latticelib.h"
#include "iolib.h"
#include <string>
#include <sys/stat.h>


int main(int argc,char ** argv)
{
	lattice alat;
	alat=ReadLattticeFromBinary_xyModel("/media/aravind/data/ACADEMIC/IISc/sem_7/FINAL_PROJECT/workbooks/cpp/xymodel/codes/simulation_data/compiled_thermalized_lattices/N16/data/1/1_1_L_16_T_1.000_J_4.400_h4_0.060_h8_0.150_xy.dat");
	alat.display_lattice_detailed();
}



















