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
//int LEN[20] = {10,12,14,16,18,22,24,26,28,30,32};
int LEN[20] = {5};
int NCOUNT=1;

double TSTART=1.5;
double TMAX=2;
double J=1;

int MC_sweep_max=25000;
int WLF_sweep_max=10000;

int sp_count=2;

char initialization='Z';
char model_type='P';
string fname;

int L,N;
double T;
double dt =0.05;

clock_t time_req;
lattice alat;

float t;
float total_time=0;
for (int j=0;j<NCOUNT;j++)
  {
	T=TSTART;
	L=LEN[j];
	N=L*L;
	alat=lattice(dim,L,N,T,J,sp_count,initialization,model_type);
	T=TSTART;
	while(1)
		{
			if(T>TMAX)
				break;
			time_req = clock();
			alat.set_temparature(T);
			cout<<"\n_______________________________________________";
			cout<<"\n\tUpdate Algo\t=\t"<<"WLF";
			cout<<"\n\tUpdate Steps\t=\t"<<WLF_sweep_max;
			cout<<"\n\tL\t\t=\t"<<L<<"\n";
			cout<<"\n\tsp_count\t=\t"<<sp_count<<"\n";
			cout<<"\tT\t\t=\t"<<T<<"  "<<alat.get_temparature()<<"\n";
			alat.wolf_potts_steps(WLF_sweep_max);
			cout<<"DONE !! writing to file .. ";
			fname="N_"+to_string(L)+"_T_"+to_string(round(T*1000)/1000).substr(0,5)+"_potts"+"WTRY";
			alat.writeToText("wolf_verify/data/potts/"+fname+".txt");
//			WriteLatticeToBinary("ising_data/"+fname+".dat",alat);
			alat.clear_data();
			time_req=clock()- time_req;
			t=float(time_req/CLOCKS_PER_SEC);
			total_time+=t;
			cout<<" DONE in "<<t<<" secs  .. estimated time left = "<<(TMAX-T)*t/dt<<" secs"<<"   Total time elapsed = "<<total_time;
			T+=dt;
		}
	alat.clear_lattice();
		
	/*alat=lattice(dim,L,N,T,J,sp_count,initialization,model_type);
	T=TSTART;
	while(1)
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
			alat.mc_potts_steps(MC_sweep_max);
			cout<<"DONE !! writing to file .. ";
			fname="N_"+to_string(L)+"_T_"+to_string(round(T*1000)/1000).substr(0,5)+"_potts"+"MCTRY";
			alat.writeToText("wolf_verify/data/potts/"+fname+".txt");
//			WriteLatticeToBinary("ising_data/"+fname+".dat",alat);
			alat.clear_data();
			time_req=clock()- time_req;
			t=float(time_req/CLOCKS_PER_SEC);
			total_time+=t;
			cout<<" DONE in "<<t<<" secs  .. estimated time left = "<<(TMAX-T)*t/dt<<" secs"<<"   Total time elapsed = "<<total_time;
			T+=dt;
		}*/
	/*while(T<1.6)
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
			alat.wolf_potts_steps(WLF_sweep_max);
			cout<<"DONE !! writing to file .. ";
			fname="N_"+to_string(L)+"_T_"+to_string(round(T*1000)/1000).substr(0,5)+"_potts";
			alat.writeToText("potts_data/"+fname+".txt");
			WriteLatticeToBinary("potts_data/"+fname+".dat",alat);
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
			alat.mc_potts_steps(MC_sweep_max);
			cout<<"DONE !! writing to file .. ";
			fname="N_"+to_string(L)+"_T_"+to_string(round(T*1000)/1000).substr(0,5)+"_potts";
			alat.writeToText("potts_data/"+fname+".txt");
			WriteLatticeToBinary("potts_data/"+fname+".dat",alat);
			alat.clear_data();
			time_req=clock()- time_req;
			t=float(time_req/CLOCKS_PER_SEC);
			total_time+=t;
			cout<<" DONE in "<<t<<" secs  .. estimated time left = "<<(TMAX-T)*t/dt<<" secs"<<"   Total time elapsed = "<<total_time;
			T+=dt;
		}*/
	alat.clear_lattice();
  }	
cout<<"\n\n";
return 0;
}
