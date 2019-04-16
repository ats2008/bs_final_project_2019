#include <iostream>
#include "isinglib.h"
#include "iolib.h"
#include <string>
using namespace std;

//int main()
int main() 
{

cout<<"\n THE ISING MODEL ROKS @ C++ !!\n";
//lattice(int dim,int lattice_len,int N,double T_equilib,double J=1.0,int spin_count=2,char initialization='Z');
int dim=2;
int L=3;
int N=L*L;
double T=8.2;
double J=1;
int sp_count=2;
char initialization='Z';
string fname;
fname="multi_histo_data/N_5_T_1.800_potts.dat";
lattice blat;
blat=ReadLattticeFromBinary(fname);

cout<<"\n THE READ ONE \n";
blat.display_lattice_detailed();
cout<<"\nEN  = "<<blat.calc_energy();
cout<<"\nMAGN  = "<<blat.calc_magnetization();
cout<<"\nEN  = "<<blat.energy.front();
cout<<"\nMAGN  = "<<blat.magnetization.front();
blat.wolf_steps(5); 
//blat.display_lattice_detailed();
//cout<<"\nEN  = "<<blat.calc_energy();
//cout<<"\nMAGN  = "<<blat.calc_magnetization();

return 0;
}
