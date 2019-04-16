#include <iostream>
#include <ctime>
#include "isinglib.h"
#include "iolib.h"
#include <string>

using namespace std;

//int main()
int main() 
{

int dim=2;
//int NCOUNT=12;
//int LEN[20] = {22,24,26,28,30,32};
int LEN[20] = {16};
int NCOUNT=1;

double TSTART=1.2;
double TMAX=3.2;
double dt =0.05;

double J=1;

int THERMALIZTION_SWEEP=10000;
int MC_sweep_max=25000;
int WLF_sweep_max=25000;

int sp_count=2;

char initialization='Z';
char model_type='I';
string fname;

int L,N;
double T;

clock_t time_req;
lattice alat;

float t;
float total_time=0;

for (int j=0;j<NCOUNT;j++)
  {
	T=TSTART;
	L=LEN[j];
	N=L*L;
	string BASE_FOLDER="multi_histo_data/N16/";
	alat=lattice(dim,L,N,T,J,sp_count,initialization,model_type);
	T=TSTART;
	cout<<"\n thermalizing the initial temperature.... \n";
	alat.mc_steps(THERMALIZTION_SWEEP);
	alat.clear_data();
	cout<<"\n DONE THE THERMALIZATION !!";
	
	while(T<1.9)
		{
			if(T>TMAX)
				break;
			time_req = clock();
			alat.set_temparature(T);
			cout<<"\n_______________________________________________";
			cout<<"\n\tUpdate Algo\t=\t"<<"METRO";
			cout<<"\n\tUpdate Steps\t=\t"<<MC_sweep_max;
			cout<<"\n\tL\t\t=\t"<<L<<"\n";
			cout<<"\n\tsp_count\t=\t"<<sp_count<<"\n";
			cout<<"\tT\t\t=\t"<<T<<"  "<<alat.get_temparature()<<"\n";
			alat.mc_steps(MC_sweep_max);
			cout<<"DONE !! writing to file .. ";
			fname="N_"+to_string(L)+"_T_"+to_string(round(T*1000)/1000).substr(0,5)+"_ising";
			alat.writeToText(BASE_FOLDER+fname+".txt");
			WriteLatticeToBinary(BASE_FOLDER+fname+".dat",alat);
			alat.clear_data();
			time_req=clock()- time_req;
			t=float(time_req/CLOCKS_PER_SEC);
			total_time+=t;
			cout<<" DONE in "<<t<<" secs  .. estimated time left = "<<(TMAX-T)*t/dt<<" secs"<<"   Total time elapsed = "<<total_time;
			T+=dt;
		}
	while( T<2.6 )
		{
			if(T>TMAX)
				break;
			time_req = clock();
			alat.set_temparature(T);
			cout<<"\n_______________________________________________";
			cout<<"\n\tUpdate Algo\t=\t"<<"WOLF";
			cout<<"\n\tUpdate Steps\t=\t"<<WLF_sweep_max;
			cout<<"\n\tL\t\t=\t"<<L<<"\n";
			cout<<"\n\tsp_count\t=\t"<<sp_count<<"\n";
			cout<<"\tT\t\t=\t"<<T<<"  "<<alat.get_temparature()<<"\n";
			alat.wolf_steps(WLF_sweep_max);
			cout<<"DONE !! writing to file .. ";
			fname="N_"+to_string(L)+"_T_"+to_string(round(T*1000)/1000).substr(0,5)+"_ising";
			alat.writeToText(BASE_FOLDER+fname+".txt");
			WriteLatticeToBinary(BASE_FOLDER+fname+".dat",alat);
			alat.clear_data();
			time_req=clock()- time_req;
			t=float(time_req/CLOCKS_PER_SEC);
			total_time+=t;
			cout<<" DONE in "<<t<<" secs  .. estimated time left = "<<(TMAX-T)*t/dt<<" secs"<<"   Total time elapsed = "<<total_time;
			T+=dt;
		}
	
	while(T<TMAX)
		{
			time_req = clock();
			alat.set_temparature(T);
			cout<<"\n_______________________________________________";
			cout<<"\n\tUpdate Algo\t=\t"<<"Metropolis";
			cout<<"\n\tUpdate Steps\t=\t"<<MC_sweep_max;
			cout<<"\n\tL\t\t=\t"<<L<<"\n";
			cout<<"\n\tsp_count\t=\t"<<sp_count<<"\n";
			cout<<"\tT\t\t=\t"<<T<<"  "<<alat.get_temparature()<<"\n";
			alat.mc_steps(MC_sweep_max);
			cout<<"DONE !! writing to file .. ";
			fname="N_"+to_string(L)+"_T_"+to_string(round(T*1000)/1000).substr(0,5)+"_ising";
			alat.writeToText(BASE_FOLDER+fname+".txt");
			WriteLatticeToBinary(BASE_FOLDER+fname+".dat",alat);
			alat.clear_data();
			time_req=clock()- time_req;
			t=float(time_req/CLOCKS_PER_SEC);
			total_time+=t;
			cout<<" DONE in "<<t<<" secs  .. estimated time left = "<<(TMAX-T)*t/dt<<" secs"<<"   Total time elapsed = "<<total_time;
			T+=dt;
		}
	alat.clear_lattice( );
  }	
cout<<"\n\n";
return 0;
}
