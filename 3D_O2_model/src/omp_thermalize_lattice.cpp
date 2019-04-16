#include <iostream>
#include <ctime>
#include "latticelib.h"
#include "iolib.h"
#include <string>
#include <sys/stat.h>
#include <chrono>
#include <thread>
#include<omp.h>

using namespace std::this_thread; 
using namespace std::chrono;
using namespace std;

int msf_write_interval=500;
int wlf_write_interval=2000;

int THERMALIZTION_SWEEP=2500;
int MSF_thermalize_step=500;
int WLF_sweep_max=10000;
int WLF_sweep_step=500;
int MSF_sweep_step=10;
int MSF_WLF_gap=1000;


char initialization='Z';
char model_type='O';
string BASE_FOLDER="simulation_data/raw_thermalized_lattices/";


int main(int argc,char ** argv)
{
		int dim=3;
		int LEN[20] = {3,12,16,20,24,28,32,36,40};
		int n=2;
		int nib_count=2*dim;

		double TEMP[10]={1};
		double H[30]={0.0,0.0,0.0,0.0,0.0};
		double Hx[30];
		int HXCOUNT=16;
		for(int i=0;i<HXCOUNT;i++)
		{
			Hx[i]=(0.01+i*0.01)*(0.01+i*0.01);
		}

		vector<string> lat_details_string,hist_names;

		double K[50]={0.4};
		double origJ=0.462;
		if(argc>6)
			origJ=atof(argv[6]);
		cout<<"\n ORGJ = "<<origJ<<"\n";
		int KCOUNT=1;
		for(int i=0;i<KCOUNT;i++)
			{
				K[i]=origJ;
			}

		string temp_str,LatSeries="";
		int LatID_GLO=0;

		if(argc>2)
		{
			LatSeries+=argv[2];
			BASE_FOLDER=argv[1];
			BASE_FOLDER+=argv[2]+string("/");
		}
		else
		{
			cout<<" please enter *BASE_FOLDER,*series_name  , *Lattice Len,  WLF_sweep_max, max number of threads, J coupling "<<"\n";
			cout<<"files will be saved under "<< BASE_FOLDER<<"BASE_FOLDER\n";
			return 0;
		}
		if(argc>3)
			LEN[0]=int(atof(argv[3]));
		if(argc>4)
			WLF_sweep_max=int(atof(argv[4]));
		


		fstream  alog;
		alog.open("log.log",ios::out|ios::app);
		cout<<"\nstarting to append new log to log.log ``__``\n";
		alog<<"\nstarting to append new log to log.log ``__``\n";
		alog.close();

		fstream  Gfnames;
		Gfnames.open(BASE_FOLDER+"series_id.txt",ios::in);
		if(!Gfnames)
		{
			cout<<temp_str<<"\nSTARTING NEW REGISTRY of Lattices in LatSeries : -> "<<LatSeries<<"<-\n";
			mkdir(BASE_FOLDER.c_str(),0777);
		}
		else
		{
			getline(Gfnames,LatSeries,'\n');
			getline(Gfnames,LatSeries,'\n');
			getline(Gfnames,temp_str,'\n');
			getline(Gfnames,temp_str,'\n');
			LatID_GLO=int(atof(temp_str.c_str()));
			LatID_GLO+=1;
		}
		Gfnames.close();
		

	int *lat_ids =new int [HXCOUNT];
	for(int i=0;i<HXCOUNT;i++)
		{
			lat_ids[i]=LatID_GLO+i;
		}
		
		int hx_id=0;
		int m=0;
		int omp_max_thread_count=30;
		if(argc>5)
		 omp_max_thread_count=int(atof(argv[5]));
		omp_set_num_threads(omp_max_thread_count);
		cout<<"\n OMP THREAD MXOUNT MAX = "<<omp_max_thread_count<<"\n\n";
		
		Gfnames.open("log.log",ios::out|ios::app);
                Gfnames<<" starting with omp_thread_coun =  "<<omp_max_thread_count<<"\n";
                Gfnames.close();

		#pragma omp parallel for schedule(dynamic)
		for (hx_id=0;hx_id<HXCOUNT;hx_id++)
		{
			string fname;
			int L,N;
			double J,T;
			T=TEMP[0];
			L=LEN[0];
			N=pow(L,dim);
			J=K[0];
			
			fstream fnames,log;
			clock_t time_req =clock();
			int wlf_step_count=0,gap_count=0;
			time_t now = time(0);
			char* dtime = ctime(&now);
			int LatID=lat_ids[hx_id];
			sleep_for(nanoseconds(100000*m));
			lattice_O_n alat;
			alat.setup_model(dim,n,L,N,nib_count,T,J,H,initialization,model_type);
			alat.seed_generator(LatID,-1);
			alat.set_temparature(T);
			H[0]=Hx[hx_id];
			alat.set_ext_field(H);
			alat.set_J(J);
			cout<<" AT "<<J<<" , "<<m<<" thermalizing for  L "<<L<<", T "<<alat.get_temparature()<<", J "<<alat.get_J();
			cout<<", H = "<<alat.get_H_as_string()<<"\n";
			
			// THERMALIZATION
			log.open("log.log",ios::out|ios::app);
			log<<"initaial MC_SWEEP[ "<<THERMALIZTION_SWEEP<<" ] @ dim = "<<alat.get_diamension()<<"  L = "<<L<<" , Temp = "<<T<<" J  = "<<J<<" H = "<<alat.get_H_as_string();
			log<<" starting\n";
			log.close();

			fname=LatSeries+"_"+to_string(LatID)+"_L_"+to_string(L)+"_T_"+to_string(round(T*1000)/1000).substr(0,5)+
															"_J_"+to_string(round(alat.get_J()*1000)/1000).substr(0,5)+
															"_H_"+alat.get_H_as_string()+"_O"+to_string(n);
			clock_t Temptime_req =clock();
			int part_count=0;
			wlf_step_count=0;
			alat.display_lattice();
			for(int i=0;i<(THERMALIZTION_SWEEP-1)/MSF_thermalize_step+1;i++)
			{	
				alat.mc_steps(MSF_thermalize_step);
//				alat.wolf_steps(MSF_thermalize_step); //
				wlf_step_count+=MSF_thermalize_step;
				if(wlf_step_count>=msf_write_interval)
				{
					part_count++;
					now = time(0);
					dtime = ctime(&now);
					alat.update_history=string(dtime);
					//# pragma omp critical
					 {
						alat.writeToText(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".txt");
						alat.writeToBinary(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".dat");
					 }
					fnames.open(string(BASE_FOLDER+"PartHfnames_"+fname+".txt"),ios::out|ios::app);
					fnames<<to_string(part_count)<<"_"+fname<<"\n";
					fnames.close();
					alat.clear_data();
					wlf_step_count=0;
				}
				cout<<"TSweep["<<THERMALIZTION_SWEEP<<"] sweep="<<MSF_thermalize_step*(i+1)<<" , @ dim = "<<alat.get_diamension()<<"  L = "<<L<<" , Temp = "<<T<<" J  = ";
				cout<<J<<" H = "<<alat.get_H_as_string()<<"\n";
				log.open("log.log",ios::out|ios::app);
				log<<"MC_SWEEP["<<THERMALIZTION_SWEEP<<"] @ dim = "<<alat.get_diamension()<<"  L = "<<L<<" , Temp = "<<T<<" J  = "<<J<<" H = "<<alat.get_H_as_string()<<" ";
				log<<" : "<<i<<" / "<<THERMALIZTION_SWEEP/MSF_thermalize_step<<"\n";
				log.close();
			}
			Temptime_req =clock()-time_req;
			log.open("log.log",ios::out|ios::app);
			log<<"MC_SWEEP["<<THERMALIZTION_SWEEP<<"] @ dim = "<<alat.get_diamension()<<"  L = "<<L<<" , Temp = "<<T<<" J  = "<<J<<" H = "<<alat.get_H_as_string()<<" ";
			log<<"&& time :: "<<float(Temptime_req/CLOCKS_PER_SEC)<<"\n";
			log.close();
			wlf_step_count=0;
			gap_count=0;
			cout<<"STARTING WLF_MSF CYCLE for "<<"@ dim = "<<alat.get_diamension()<<"  L = "<<L<<" , Temp = "<<T<<" J  = "<<J<<" H = "<<alat.get_H_as_string()<<"\n";;
			for(int qw=0;qw<int((WLF_sweep_max-1)/WLF_sweep_step)+1;qw++)
			{
//				alat.wolf_steps(WLF_sweep_step);
				alat.mc_steps(WLF_sweep_step);
				wlf_step_count+=WLF_sweep_step;
				gap_count+=WLF_sweep_step;
				if(wlf_step_count>=wlf_write_interval)
				{
					part_count++;
					now = time(0);
					dtime = ctime(&now);
					alat.update_history=string(dtime);
					//# pragma omp critical
					 {
						alat.writeToText(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".txt");
						alat.writeToBinary(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".dat");
					 }
					fnames.open(BASE_FOLDER+"PartHfnames_"+fname+".txt",ios::out|ios::app);
					fnames<<part_count<<"_"+fname<<"\n";
					fnames.close();
					alat.clear_data();
					wlf_step_count=0;
				}
				if(gap_count>=MSF_WLF_gap)
				{
					alat.mc_steps( MSF_sweep_step );
//					alat.wolf_steps(MSF_thermalize_step); //
					wlf_step_count+=MSF_sweep_step;
					if (wlf_step_count >= wlf_write_interval)
					{
						part_count++;
						now = time(0);
						dtime = ctime(&now);
						alat.update_history=string(dtime);
//						# pragma omp critical
						 {
							alat.writeToText(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".txt");
							alat.writeToBinary(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".dat");
						 }
 						fnames.open(BASE_FOLDER+"PartHfnames_"+fname+".txt",ios::out|ios::app);
						fnames<<part_count<<"_"+fname<<"\n";
						fnames.close();
						alat.clear_data();
						wlf_step_count=0;
					}
					gap_count=0;
				}
				
				Temptime_req = clock() - time_req;
				log.open("log.log",ios::out|ios::app);
				log<<"WLF_SWEEP["<<WLF_sweep_max<<"] sweeps = "<<WLF_sweep_step*(qw+1)<<", @ dim = "<<alat.get_diamension();
				log<<"  L = "<<L<<" , Temp = "<<T<<" J  = "<<J<<" H = "<<alat.get_H_as_string()<<"["<<qw+1<<"/";
				log<<int((WLF_sweep_max-1)/WLF_sweep_step)<<"]  && time :: "<<float(Temptime_req/CLOCKS_PER_SEC)<<"\n";
				log.close();
				cout<<"WLF_SWEEP["<<WLF_sweep_step<<"] sweeps = "<<WLF_sweep_step*(qw+1)<<", @ dim = "<<alat.get_diamension();
				cout<<"  L = "<<L<<" , Temp = "<<T<<" J  = "<<J<<" H = "<<alat.get_H_as_string()<<"["<<qw+1<<"/";
				cout<<int((WLF_sweep_max-1)/WLF_sweep_step)<<"]  && time :: "<<float(Temptime_req/CLOCKS_PER_SEC)<<"\n";
				
			}
			
			
			log.open("log.log",ios::out|ios::app);
			log<<"@ dim = "<<alat.get_diamension()<<"  L = "<<L<<" , Temp = "<<T<<" .. -> done wlf cycles \n";
			log<<"WLF_SWEEP DONE FOR  @ dim = "<<alat.get_diamension()<<"  L = "<<L<<" , Temp = "<<T<<" J  = "<<J<<" H = "<<alat.get_H_as_string()<<"\n";
			log.close();
			cout<<"DONE !! now writing to file .. ";
			part_count++;
			now = time(0);
			dtime = ctime(&now);
			alat.update_history=string(dtime);
			//# pragma omp critical
			 {
				alat.writeToText(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".txt");
				alat.writeToBinary(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".dat");
			 }
			fnames.open(BASE_FOLDER+"PartHfnames_"+fname+".txt",ios::out|ios::app);
			fnames<<part_count<<"_"+fname<<"\n";
			fnames.close();
			log.open("log.log",ios::out|ios::app);
			log.close();
			alat.clear_lattice( );
			# pragma omp critical
			{
				fnames.open(BASE_FOLDER+"lat_names.txt",ios::out|ios::app);
				fnames<<fname<<"\n";
				fnames.close();
			}
		}
		  
		Gfnames.open(BASE_FOLDER+"series_id.txt",ios::out);
		Gfnames<<"#seriesName\n"<<LatSeries<<"\n";
		Gfnames<<"#LAT_ID\n"<<LatID_GLO+HXCOUNT<<"\n";
		Gfnames.close();
		cout<<"\n_______^^^__________^^__________\n";
		
		delete [] lat_ids;
		return 0;
}
