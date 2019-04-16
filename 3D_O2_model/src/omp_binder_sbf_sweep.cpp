
#include <iostream>
#include <ctime>
#include "latticelib.h"
#include "iolib.h"
#include <string>
#include <sys/stat.h>
#include<omp.h>

using namespace std;


void get_NthOrder_spin_expectaion(lattice alat ,int n,double *expec );
void spawn_child_lattices(lattice plat,int id,int lvl,int BASE_SEED,int spawn_max,int thread_count_max,string BASE_FOLDER);
int dim=2;
int LEN[20] = {16,64,12,16,20,24,28,32,36,40};
int NCOUNT=1;

double TEMP[100]={1};
int TCOUNT=1;

double K[50]={5.10,5.0,4.95,4.93,4.90,4.80,4.97,4.70};
int KCOUNT=1;
int K_BEGI=0;

double H4[30]={0.06};
int H4COUNT=1;

double H8[10]={0.15};
int H8COUNT =1 ;


int RG_REDUCTION_FACTOR[20]={1,2,4,8,12,16};
int RG_REDUCTION_COUNT=6;

//MOMENTS SCHEME SETUP
int MOMENTS[30]={1,2,4,6,8,12,16};
int MOMENTS_COUNT=4;

int THERMALIZTION_SWEEP=0;
int MC_sweep_step=10;
int WLF_sweep_step=100;
int WLF_sweep_max=4000;
int WLF_MSF_GAP =1000;

char initialization='Z';
char model_type='X';

//Spawning realted vars
int SPAWN_MAX = 3;

int main(int argc,char ** argv)
{
		string fname;
		string BASE_FOLDER="simulation_data/compiled_thermalized_lattices/";
		string DESTN_FOLDER="";

	    if(argc>1)
		{
			BASE_FOLDER+=argv[1]+string("/");
//			mkdir(BASE_FOLDER.c_str(),0777);
		}
		else
		{
			cout<<"Please provide destination name\n";
			return 0;
		}
		if(argc>2)
		{
			WLF_sweep_max=int(atof(argv[2]));
		}
		if(argc>3)
		{
			DESTN_FOLDER=argv[3];
			DESTN_FOLDER+="data/";
			mkdir((DESTN_FOLDER).c_str(),0777);
		}
		else
		{
			DESTN_FOLDER=BASE_FOLDER;
			DESTN_FOLDER+="data/";
			mkdir(DESTN_FOLDER.c_str(),0777);
		}


		fstream  log;
		log.open("log.log",ios::out|ios::app);
		log<<"starting to append new log to log.log  |^_^| \n";
		log.close();
		cout<<"starting to append new log to log.log  |^_^| \n";
		cout<<" base folder : "<<BASE_FOLDER<<"\n";
		vector<string> seed_names;
		fstream  ifile;
		ifile.open(BASE_FOLDER+"lattices/latnames.txt",ios::in);
		string l="";
		while( getline (ifile,l,',') )
		{
			seed_names.push_back(l+string(".dat"));
			cout<<" adding to seed_list : "<<l<<"\n";
			getline (ifile,l,'\n'); 
		
		}
		ifile.close();
	
	
		int BASE_SEED=0;
		int data_len_max=20000;
		ifile.open(BASE_FOLDER+"lattices/config.txt",ios::in);
		if(ifile)
		{	
			getline (ifile,l,'\n') ;
			getline (ifile,l,'\n') ;
			BASE_SEED=int(atof(l.c_str()));
			getline (ifile,l,'\n') ;
			getline (ifile,l,'\n') ;
			data_len_max=int(atof(l.c_str()));
			cout<<"\n read base seed and data_len_max "<<BASE_SEED<<" , "<<data_len_max<<"\n";
		}
		else
		{
			ifile.open(BASE_FOLDER+"lattices/config.txt",ios::out);
			ifile<<"#BASE_SEED"<<"\n";;
			ifile<<BASE_SEED<<"\n";
			ifile<<"#data_len_max"<<"\n";
			ifile<<data_len_max<<"\n";
		}
		ifile.close();
		// READ_INITIAL_LATTICES
		vector<string>::iterator name = seed_names.begin();
		vector<lattice> seed_lats;
		for(;name!=seed_names.end();name++)
		{
			seed_lats.push_back(ReadLattticeFromBinary_xyModel(BASE_FOLDER+"lattices/"+*name));
		}

		//set the thread count
		int thread_count_max = int((data_len_max-1)/WLF_sweep_max)+1+1;
		//Spawn the tree of lattices
		int seed_len=seed_lats.size();
		cout<<"BASE_SEED        = "<<BASE_SEED<<"\n";
		cout<<"data_len_max     = "<<data_len_max<<"\n";
		cout<<"seed_len         = "<<seed_len<<"\n";
		cout<<"thread_count_max = "<<thread_count_max<<"\n";
		
		omp_set_nested(1);
		#pragma omp parallel for schedule(dynamic)
		for (int i=0;i<seed_len;i++)
		{
			//seed_lats[i].display_lattice_detailed();
			cout<<"\n seeding base lattices  : seed = "<<BASE_SEED+thread_count_max*i;
			cout<<" , SPAWN_MAX = "<<SPAWN_MAX<<" ,thread_count_max : "<<thread_count_max<<"\n";
			spawn_child_lattices(seed_lats[i],1,0,BASE_SEED+thread_count_max*i,SPAWN_MAX,thread_count_max,DESTN_FOLDER);
			//seed_lats[i].clear_lattice();
			cout<<"\n\n";
		}
		cout<<"\n\n";
		return 0;
}

void get_NthOrder_spin_expectaion(lattice alat ,int n,double *expec )
{
	site * state =alat.get_state();
	expec[0]=0;
	expec[1]=0;
	expec[2]=0;
	for(int i=0;i<alat.get_site_count();i++)
	{
		expec[0]+=pow((state+i)->_components[1]*(state+i)->_components[1]+(state+i)->_components[0]*(state+i)->_components[1],n/2);
		expec[1]+=pow((state+i)->_components[0],n);
		expec[2]+=pow((state+i)->_components[1],n);
	}
}

void spawn_child_lattices(lattice plat,int id,int lvl,int BASE_SEED,int spawn_max,int thread_count_max,string BASE_FOLDER)
{
	int curr_gen=int(pow(spawn_max,lvl));
	int total_count=0;
	for(int i=0;i<=lvl;i++) total_count+=int(pow(spawn_max,i));
	int shift = curr_gen-(total_count-id)-1;
	int child_id_base = total_count+shift*spawn_max+1;
	
	
//	cout<<"seed = "<<id<<"  ,";
//	cout<<"level = "<<lvl<<"  ,";
//	cout<<"total_count = "<<total_count<<"  ,";
//	cout<<"curr_gen = "<<curr_gen<<"  ,";
//	cout<<"shift = "<<shift<<"  ,";
//	cout<<"child_base = "<<child_id_base<<"  ,";
//	cout<<"SEED = "<<id+BASE_SEED<<"\n";

	fstream fnames;
//	if(id>thread_count_max) return;
  #pragma omp parallel for schedule(dynamic)
	for(int i=0;i<spawn_max;i++)
	{	
		int idx,idy;
		int child_id=child_id_base+i;
		time_t now = time(0);
		char* dtime = ctime(&now);
		if(child_id>thread_count_max) 
		{
			cout<<"Skipping spawning as child_id>thread_count_max c id [seed]:"<<child_id<<" [ "<<child_id+BASE_SEED<<" ] ";
			cout<<" MAX : "<<thread_count_max<<"\n";
			continue;
		}
		cout<<"seeding child with id : "<<child_id<<" [ "<<child_id+BASE_SEED<<" ]  MAX : "<<thread_count_max<<"\n";
	//cout<<"\n";
		fstream log,ofile;
	//cout<<"bare\n";
		lattice alat=plat.get_bare_copy();
	//cout<<"bare got\n";
		string fname;
		vector<string> MomentNames;
		alat.seed_generator(child_id+BASE_SEED,-1);
		clock_t time_req;
		time_req = clock( );
		
		// Moment Scheme setup
		vector< vector<double * > >* SPIN_moments;
		vector< vector<double * > >::iterator moment_itr;
		
		//RG SCHEME SETUP
		vector< lattice > RG_lattices;
		vector< lattice >::iterator lat_itr;

		// SETING UP THE RG LATTICES from Parent lattice
		SPIN_moments = new vector< vector<double *> > [RG_REDUCTION_COUNT];
		
		double * exec=NULL;
		for(int rgc=0;rgc<RG_REDUCTION_COUNT;rgc++)
		{
			RG_lattices.push_back(alat.get_RG_blocked_lattice(RG_REDUCTION_FACTOR[rgc]));
			for(int mc=0;mc<MOMENTS_COUNT;mc++)
				SPIN_moments[rgc].push_back(*(new vector<double *>));
		}
//		cout<<"Starting for  L "<<alat.get_lattice_len()<<", T "<<alat.get_temparature()<<", J "<<alat.get_J();
//		cout<<", h4 "<<alat.get_h4()<<", h8 "<<alat.get_h8()<<"\n";
		log.open("log.log",ios::out|ios::app);
		log<<child_id<<" : Starting for  L "<<alat.get_lattice_len()<<", T "<<alat.get_temparature()<<", J "<<alat.get_J();
		log<<", h4 "<<alat.get_h4()<<", h8 "<<alat.get_h8()<<" time "<<0<<"\n";		
		log.close();
//		alat.display_lattice_detailed();
		// DOING THE MONTE CARLO SIMULATION CALLS with Wolf Steps ( with intermediate Metro. Step)
		int wlf_gap=0;
		int write_gap=0;
		for(int qw=0;qw<int((WLF_sweep_max-1)/WLF_sweep_step)+1;qw++)
		{
		 wlf_gap+=WLF_sweep_step;
		 write_gap+=WLF_sweep_step;
			for(int qe=0;qe<WLF_sweep_step;qe++)
			{
				alat.wolf_steps_xyModel(1);
				idx=0;
				for( lat_itr = RG_lattices.begin()
									; lat_itr!=RG_lattices.end(); lat_itr ++)
						{
							RG_block_xy_model(alat,*lat_itr,RG_REDUCTION_FACTOR[idx]);
							idy=0;
							for(moment_itr=SPIN_moments[idx].begin();
											moment_itr!=SPIN_moments[idx].end(); moment_itr++)
								{
									exec = new double [3];
									get_NthOrder_spin_expectaion(*lat_itr,MOMENTS[idy],exec);
									moment_itr->push_back(exec);
									idy++;
								}
							idx++;
						}
			}
			// DOING THE MONTE CARLO SIMULATION CALLS with Metropolis Single flip Steps
			if(wlf_gap> WLF_MSF_GAP)
			{
				wlf_gap=0;
				for(int qe=0;qe<MC_sweep_step;qe++)
				{
					alat.mc_steps_xyModel(1);
					idx=0;
					for( lat_itr = RG_lattices.begin()
										; lat_itr!=RG_lattices.end(); lat_itr ++)
							{
								RG_block_xy_model(alat,*lat_itr,RG_REDUCTION_FACTOR[idx]);
								idy=0;
								for(moment_itr=SPIN_moments[idx].begin();
												moment_itr!=SPIN_moments[idx].end(); moment_itr++)
									{
										exec = new double [3];
										get_NthOrder_spin_expectaion(*lat_itr,MOMENTS[idy],exec);
										moment_itr->push_back(exec);
										idy++;
									}
								idx++;
							}
				}
			}
			log.open("log.log",ios::out|ios::app);
			log<<child_id<<" @ L = "<<alat.get_lattice_len()<<", T "<<alat.get_temparature()<<", J "<<alat.get_J();
			log<<", h4 "<<alat.get_h4()<<", h8 "<<alat.get_h8()<<" on ";
			log<<qw+1<<" / "<<int(WLF_sweep_max/WLF_sweep_step)+1<<"  [ ";
			log<<float((clock()- time_req)/CLOCKS_PER_SEC)<<" ] \n";
			log.close();

		}


		// WRITING LATTICE TO FILE along with data
//		cout<<"DONE !! now writing to file .. ";
		fname=to_string(child_id)+"_L_"+to_string(alat.get_lattice_len())+
														"_T_"+to_string(round(alat.get_temparature()*1000)/1000).substr(0,5)+
														"_J_"+to_string(round(alat.get_J()*1000)/1000).substr(0,5)+
														"_h4_"+to_string(round(alat.get_h4()*1000)/1000).substr(0,5)+
														"_h8_"+to_string(round(alat.get_h8()*1000)/1000).substr(0,5)+"_xy";
		now = time(0);
		dtime = ctime(&now);
		alat.update_history=string(dtime);
		alat.writeToText_xyModel(BASE_FOLDER+fname+".txt");
		WriteLatticeToBinary_xyModel(BASE_FOLDER+fname+".dat",alat);

		//WRITING RG MOMENTS TO FILES
//		cout<<"writing rg things now \n";
		idx=0;
		for( lat_itr = RG_lattices.begin()
									; lat_itr!=RG_lattices.end(); lat_itr ++)
			{
				ofile.open(BASE_FOLDER+fname+"_RG_"+to_string(RG_REDUCTION_FACTOR[idx])+".txt",ios::out);
				MomentNames.push_back(fname+"_RG_"+to_string(RG_REDUCTION_FACTOR[idx]));
				ofile<<"#L , "<<alat.get_lattice_len()<<"\n";
				ofile<<"#T , "<<alat.get_temparature()<<"\n";
				ofile<<"#J , "<<alat.get_J()<<"\n";
				ofile<<"#h4 , "<<alat.get_h4()<<"\n";
				ofile<<"#h8 , "<<alat.get_h8()<<"\n";
				ofile<<"#RG , "<<RG_REDUCTION_FACTOR[idx]<<"\n";
				ofile<<"id";
				for(int mc=0;mc<MOMENTS_COUNT;mc++)
					ofile<<",M_s"<<MOMENTS[mc]<<",Mx_s"<<MOMENTS[mc]<<",My_s"<<MOMENTS[mc];
				ofile<<"\n";
				moment_itr=SPIN_moments[idx].begin();
				int data_len=moment_itr->size();
				idy=0;
				while(idy<data_len)
				{
					ofile<<idy;
					for(moment_itr=SPIN_moments[idx].begin();
								moment_itr!=SPIN_moments[idx].end(); moment_itr++)
					{
						ofile<<","<<(*moment_itr)[idy][0]<<","<<(*moment_itr)[idy][1]<<","<<(*moment_itr)[idy][2];
						delete [] (*moment_itr)[idy];
					}
					ofile<<"\n";
				idy++;
				}
				idx++;
				ofile.close();
			}

		//clearing RG sweep data
		idx=0;
		for( idx = 0; idx< RG_REDUCTION_COUNT; idx ++)
			for(moment_itr=SPIN_moments[idx].begin();
								moment_itr!=SPIN_moments[idx].end(); moment_itr++)
					{
						moment_itr->clear();
					}
		delete [] SPIN_moments;
		for( lat_itr = RG_lattices.begin()
										; lat_itr!=RG_lattices.end(); lat_itr ++)
			{
				lat_itr->clear_lattice();
			}
			
		log.open("log.log",ios::out|ios::app);
		log<<"completed "<<child_id<<"  L "<<alat.get_lattice_len()<<", T "<<alat.get_temparature()<<", J "<<alat.get_J();
		log<<", h4 "<<alat.get_h4()<<", h8 "<<alat.get_h8()<<"  [ ";
		log<<float((clock()- time_req)/CLOCKS_PER_SEC)<<" ] \n";
		log.close();
		
	//clearing Main Lattice Sweep data
	alat.clear_data( );
//void spawn_lattices(lattice plat,int id,int lvl,int BASE_SEED,int spawn_max,int thread_count_max,string BASE_FOLDER);
	spawn_child_lattices(alat,child_id,lvl+1,BASE_SEED,spawn_max,thread_count_max,BASE_FOLDER);
	alat.clear_lattice();
		#pragma omp critical
		{
			fnames.open(BASE_FOLDER+"fnames.txt",ios::out|ios::app);
			fnames<<fname<<".txt\n";
			fnames.close();
		}
		#pragma omp critical
		{
			fnames.open(BASE_FOLDER+"pdetails.txt",ios::out|ios::app);
			fnames<<fname<<",";
			fname=to_string(id);
			fname+="_L_"+to_string(alat.get_lattice_len())+
														"_T_"+to_string(round(alat.get_temparature()*1000)/1000).substr(0,5)+
														"_J_"+to_string(round(alat.get_J()*1000)/1000).substr(0,5)+
														"_h4_"+to_string(round(alat.get_h4()*1000)/1000).substr(0,5)+
														"_h8_"+to_string(round(alat.get_h8()*1000)/1000).substr(0,5)+"_xy";	
			fnames<<fname<<"\n";
			fnames.close();
		}
		#pragma omp critical
		{
			fnames.open(BASE_FOLDER+"Mfnames.txt",ios::out|ios::app);
			for(vector<string>::iterator itr=MomentNames.begin();itr!=MomentNames.end();itr++)
				fnames<<*itr<<"\n";
			fnames.close();
		}
	
	}

}

