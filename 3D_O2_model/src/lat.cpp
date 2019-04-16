#include <iostream>
#include <ctime>
#include "latticelib.h"
#include "iolib.h"
#include <string>
#include <sys/stat.h>
#include <chrono>
#include <thread>

using namespace std::this_thread; 
using namespace std::chrono;
using namespace std;

int msf_write_interval=40;
int wlf_write_interval=200;

int THERMALIZTION_SWEEP=2500;
int MSF_thermalize_step=100;
int WLF_sweep_max=2000;
int WLF_sweep_step=100;
int MSF_sweep_step=10;
int MSF_WLF_gap=800;


char initialization='H';
char model_type='X';
string BASE_FOLDER="simulation_data/thermalized_lattices/";


int main(int argc,char ** argv)
{
		double H[2]={0.0,0.0};
		int dim = 2;
		int n=2;
		int L = 5;
		double T = 1;
		double J = 1;
		int nib_count =4 ;
		int N = int(pow(L,dim));
		char model_type = 'O';
		char initialization='Z';
		
		lattice_O_n alat;
		alat.setup_model(dim,n,L,N,nib_count,T,J,H,initialization,model_type);
		alat.display_lattice_detailed();
		cout<<"\n ";
		alat.display_neighbours();
		
		alat.mc_steps(5000);
		alat.wolf_steps(5000);
//		alat.mc_steps(5000);
//		alat.wolf_steps(5000);
//		alat.mc_steps(5000);
		
		cout<<"\n\nwriting txt\n";
		alat.writeToText("alat.txt");
		cout<<"writing binary\n";
		alat.writeToBinary("alat.dat");
		cout<<"reading the binary\n";
		lattice_O_n blat;
		blat.readFromBinary("alat.dat");
		cout<<"newlat \n";
		blat.display_lattice_detailed();

		alat.clear_lattice();
		blat.clear_lattice();
		
		return 0;
}
