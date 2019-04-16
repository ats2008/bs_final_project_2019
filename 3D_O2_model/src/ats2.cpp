
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
		int dim=2;
		int LEN[20] = {5,12,16,20,24,28,32,36,40};

		double TEMP[100]={1};
		double H4[30]={0.06};
		double H8[10]={0.15};
		vector<string> lat_details_string,hist_names;

		double K[50]={1};
		int KCOUNT=2;
		for(int i=0;i<KCOUNT;i++)
			{
				K[i]=4.2+0.2*i;
			}

		string temp_str,LatSeries="";
		int LatID_GLO=0;

		if(argc>1)
		{
			LatSeries+=argv[1];
			BASE_FOLDER+=argv[1]+string("/");
		}
		if(argc>2)
			LEN[0]=int(atof(argv[2]));
		if(argc>3)
			WLF_sweep_max=int(atof(argv[3]));
		
		int L,N;
		double T,h4,h8;


		fstream  alog;
		alog.open("log.log",ios::out|ios::app);
		cout<<"\nstarting to append new log to log.log ``__``\n";
		alog<<"\nstarting to append new log to log.log ``__``\n";
		alog.close();

		fstream  Gfnames;
		Gfnames.open(BASE_FOLDER+"series_id.txt",ios::in);
		if(!Gfnames)
		{
			cout<<temp_str<<"\nSTARTING NEW REGISTRY of Lattices in LatSeries : ->"<<LatSeries<<"<-\n";
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
		T=TEMP[0];
		L=LEN[0];
		N=L*L;
		h4=H4[0];
		h8=H8[0];


	int *lat_ids =new int [KCOUNT];
	for(int i=0;i<KCOUNT;i++)
		{
			lat_ids[i]=LatID_GLO+i;
		}

	#pragma omp parallel for 
		for (int m=0;m<KCOUNT;m++)
		{
			string fname;
			fstream fnames,log;
			clock_t time_req =clock();
			int wlf_step_count=0,gap_count=0;
			time_t now = time(0);
			char* dtime = ctime(&now);
			double J=K[m];
			int LatID=lat_ids[m];
			sleep_for(nanoseconds(100000*m));
			lattice alat;
			alat.setup_xy_model(dim,L,N,T,J,h4,h8,initialization,model_type);
			alat.seed_generator(LatID,-1);
			alat.set_temparature(T);
			alat.set_h4(h4);
			alat.set_h8(h8);
			alat.set_J(K[m]);
			cout<<" AT "<<J<<" , "<<m<<" thermalizing for  L "<<L<<", T "<<alat.get_temparature()<<", J "<<alat.get_J();
			cout<<", h4 "<<alat.get_h4()<<", h8 "<<alat.get_h8()<<"\n";
			
			// THERMALIZATION
			log.open("log.log",ios::out|ios::app);
			log<<"initaial MC_SWEEP[ "<<THERMALIZTION_SWEEP<<" ] @ L = "<<L<<" , Temp = "<<T<<" J  = "<<J<<" h4 = "<<h4<<" h8="<<h8<<" starting\n";
			log.close();

			fname=LatSeries+"_"+to_string(LatID)+"_L_"+to_string(L)+"_T_"+to_string(round(T*1000)/1000).substr(0,5)+
															"_J_"+to_string(round(alat.get_J()*1000)/1000).substr(0,5)+
															"_h4_"+to_string(round(alat.get_h4()*1000)/1000).substr(0,5)+
															"_h8_"+to_string(round(alat.get_h8()*1000)/1000).substr(0,5)+"_xy";
			clock_t Temptime_req =clock();
			int part_count=0;
			wlf_step_count=0;
			for(int i=0;i<(THERMALIZTION_SWEEP-1)/MSF_thermalize_step+1;i++)
			{	
				alat.mc_steps_xyModel(MSF_thermalize_step);
				wlf_step_count+=MSF_thermalize_step;
				if(wlf_step_count>=msf_write_interval)
				{
					part_count++;
					now = time(0);
					dtime = ctime(&now);
					alat.update_history=string(dtime);
					# pragma omp critical
					 {
						alat.writeToText_xyModel(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".txt");
						WriteLatticeToBinary_xyModel(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".dat",alat);
					 }
//					fnames.open(string(BASE_FOLDER+"PartHfnames_"+fname+".txt"),ios::out|ios::app);
//					fnames<<to_string(part_count)<<"_"+fname<<"\n";
//					fnames.close();
					alat.clear_data();
					wlf_step_count=0;
				}
				cout<<"TSweep["<<THERMALIZTION_SWEEP<<"] sweep="<<MSF_thermalize_step*(i+1)<<" , "<<"] @ L = "<<L<<" , Temp = "<<T<<" J  = ";
				cout<<J<<" h4 = "<<h4<<" h8="<<h8<<"\n";
				log.open("log.log",ios::out|ios::app);
				log<<"MC_SWEEP["<<THERMALIZTION_SWEEP<<"] @ L = "<<L<<" , Temp = "<<T<<" J  = "<<J<<" h4 = "<<h4<<" h8="<<h8<<" ";
				log<<" : "<<i<<" / "<<THERMALIZTION_SWEEP/MSF_thermalize_step<<"\n";
				log.close();
			}
			Temptime_req =clock()-time_req;
			log.open("log.log",ios::out|ios::app);
			log<<"MC_SWEEP["<<THERMALIZTION_SWEEP<<"] @ L = "<<L<<" , Temp = "<<T<<" J  = "<<J<<" h4 = "<<h4<<" h8="<<h8<<" ";
			log<<"&& time :: "<<float(Temptime_req/CLOCKS_PER_SEC)<<"\n";
			log.close();
			wlf_step_count=0;
			gap_count=0;
			cout<<"STARTING WLF_MSF CYCLE for "<<"@ L = "<<L<<" , Temp = "<<T<<" J  = "<<J<<" h4 = "<<h4<<" h8="<<h8<<" ";;
			for(int qw=0;qw<int((WLF_sweep_max-1)/WLF_sweep_step)+1;qw++)
			{
				alat.wolf_steps_xyModel(WLF_sweep_step);
				wlf_step_count+=WLF_sweep_step;
				gap_count+=WLF_sweep_step;
				if(wlf_step_count>=wlf_write_interval)
				{
					part_count++;
					now = time(0);
					dtime = ctime(&now);
					alat.update_history=string(dtime);
					# pragma omp critical
					 {
						alat.writeToText_xyModel(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".txt");
						WriteLatticeToBinary_xyModel(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".dat",alat);
					 }
					fnames.open(BASE_FOLDER+"PartHfnames_"+fname+".txt",ios::out|ios::app);
					fnames<<part_count<<"_"+fname<<"\n";
					fnames.close();
					alat.clear_data();
					wlf_step_count=0;
				}
				if(gap_count>=MSF_WLF_gap)
				{
					alat.mc_steps_xyModel( MSF_sweep_step );
					wlf_step_count+=MSF_sweep_step;
					if (wlf_step_count >= wlf_write_interval)
					{
						part_count++;
						now = time(0);
						dtime = ctime(&now);
						alat.update_history=string(dtime);
						# pragma omp critical
						 {
							alat.writeToText_xyModel(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".txt");
							WriteLatticeToBinary_xyModel(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".dat",alat);
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
				log<<"WLF_SWEEP["<<WLF_sweep_step<<"] @ L = "<<L<<" , Temp = "<<T<<" J  = "<<J<<" h4 = "<<h4<<" h8="<<h8<<"["<<qw+1<<"/";
				log<<int((WLF_sweep_max-1)/WLF_sweep_step)<<"]  && time :: "<<float(Temptime_req/CLOCKS_PER_SEC)<<"\n";
				log.close();
				cout<<"WLF_SWEEP["<<WLF_sweep_step<<"] @ L = "<<L<<" , Temp = "<<T<<" J  = "<<J<<" h4 = "<<h4<<" h8="<<h8<<"["<<qw+1<<"/";
				cout<<int((WLF_sweep_max-1)/WLF_sweep_step)<<"]  && time :: "<<float(Temptime_req/CLOCKS_PER_SEC)<<"\n";
				
			}
			
			
			log.open("log.log",ios::out|ios::app);
			log<<"@ L = "<<L<<" , Temp = "<<T<<" .. -> done wlf cycles \n";
			log<<"WLF_SWEEP DONE FOR  @ L = "<<L<<" , Temp = "<<T<<" J  = "<<J<<" h4 = "<<h4<<" h8="<<h8<<"\n";
			log.close();
			cout<<"DONE !! now writing to file .. ";
			part_count++;
			now = time(0);
			dtime = ctime(&now);
			alat.update_history=string(dtime);
			# pragma omp critical
			 {
				alat.writeToText_xyModel(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".txt");
				WriteLatticeToBinary_xyModel(BASE_FOLDER+"history_"+to_string(part_count)+"_"+fname+".dat",alat);
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
		Gfnames<<"#LAT_ID\n"<<LatID_GLO+KCOUNT<<"\n";
		Gfnames.close();
		cout<<"\n_______^^^__________^^__________\n";
		
		delete [] lat_ids;
		return 0;
}
