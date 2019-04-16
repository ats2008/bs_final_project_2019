#include <iostream>
#include <ctime>
#include "latticelib.h"
//#include "iolib.h"
#include <string>
#include<omp.h>
#include <sys/stat.h>

void get_NthOrder_spin_expectaion(lattice alat ,int n,double *expec );

double K[50]={5.10,5.0,4.95,4.93,4.90,4.80,4.97,4.70};
int KCOUNT=1;
int K_BEGI=0;

char initialization='Z';
char model_type='X';


//MOMENTS SCHEME SETUP
int MOMENTS[30]={1,2,4,6,8,12,16};
int MOMENTS_COUNT=4;

int THERMALIZTION_SWEEP=0;
int MC_sweep_step=10;
int WLF_sweep_step=200;
int WLF_sweep_max=4000;
int WLF_MSF_GAP =1000;

// DATA MANAGEMENT SETUP
int  FILE_WRITE_GAP = 400;
int LATTICE_WRITE_GAP =2000;


int main(int argc,char ** argv)
{

		string fname;
		string BASE_FOLDER="simulation_data/compiled_thermalized_lattices/";
		string DESTN_FOLDER="";

	    if(argc>1)
		{
			BASE_FOLDER=argv[1];
//			mkdir(BASE_FOLDER.c_str(),0777);
		}
		else
		{
			cout<<"Please provide *source_path ,destn_path,destn_folder_name,multiplicity_of_seed_lats,data_len_max,";
			cout<<"omp_max_thread_count,start,end\n";
			cout<<"if destn_path not set files are tobe saved under source_path/data\n";
			cout<<"if destn_folder_name not set files are under source_path/data\n\n";
			return 0;
		}
		if(argc>3)
		{
			DESTN_FOLDER=argv[2];
			DESTN_FOLDER+=argv[3];
			DESTN_FOLDER+="/";
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
		cout<<" destn folder : "<<DESTN_FOLDER<<"\n";

		// READ_INITIAL_LATTICE names
		vector<string> seed_names;
		fstream  ifile;
		ifile.open(BASE_FOLDER+"lattices/latnames.txt",ios::in);
		int multiplicity_of_seed_lats=1;
		if(argc>4)
		{
			multiplicity_of_seed_lats=int(atof(argv[4]));
		}
		string l="";
		while( getline (ifile,l,',') )
		{
			for(int i=0;i<multiplicity_of_seed_lats;i++)
				seed_names.push_back(l+string(".dat"));
			cout<<" adding to seed_list : "<<l<<" x "<<multiplicity_of_seed_lats<<"\n";
			getline (ifile,l,'\n'); 
		
		}
		ifile.close();
		
		// READING BASIC CONFIGS
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
			ifile.close();
			ifile.open(BASE_FOLDER+"lattices/config.txt",ios::out);
			ifile<<"#BASE_SEED"<<"\n";;
			ifile<<BASE_SEED<<"\n";
			ifile<<"#data_len_max"<<"\n";
			ifile<<data_len_max<<"\n";
		}
		ifile.close();
		
		int seed_len=seed_names.size();
		// SETUP LATIDs
		int *lat_ids = new int [seed_len];
		for(int i=0;i<seed_len;i++)
			lat_ids[i]=BASE_SEED+i;
			
		ifile.open(BASE_FOLDER+"lattices/config.txt",ios::out);
		ifile<<"#BASE_SEED"<<"\n";
		ifile<<BASE_SEED+seed_len<<"\n";
		ifile<<"#data_len_max"<<"\n";
		ifile<<data_len_max<<"\n";
		ifile.close();
		
		if(argc>5)
		{
			data_len_max=int(atof(argv[5]));
		}

		int start=0,end=seed_names.size();
		if(argc>7)
		{
			start=int(atof(argv[7]));
		}
		if(argc>8)
		{
			end=int(atof(argv[8]));
		}
		
		//Basic simulation stats : 
		cout<<"\n\n";
		cout<<"BASE_SEED        = "<<BASE_SEED<<"\n";
		cout<<"data_len_max     = "<<data_len_max<<"\n";
		cout<<"seed_len         = "<<seed_len<<"\n";
		
		
		WLF_sweep_max=data_len_max;
		int omp_max_thread_count=3;
		if(argc>6)
		{
			omp_max_thread_count=int(atof(argv[6]));
		}
		omp_set_num_threads(omp_max_thread_count);
//		#pragma omp parallel for schedule(dynamic)
		for (int m=start;m<end;m++)
		{
			cout<<"\n"<<m<<"\n";
			lattice_O_n current_lattice;
			current_lattice.readFromBinary(BASE_FOLDER+"lattices/"+seed_names[m]);
			string DESTN_FOLDER_local,fname,partfname;
			fstream fnames,log,ofile;
			clock_t time_req =clock();
			clock_t Temptime_req;
			time_t now = time(0);
			char* dtime = ctime(&now);
			int LatID=lat_ids[m];
			int lattice_seed=LatID;
			current_lattice.seed_generator(LatID);
			
			DESTN_FOLDER_local=DESTN_FOLDER+"LatID_"+to_string(LatID)+"/";
			mkdir(DESTN_FOLDER_local.c_str(),0777);
			cout<<"made folder = "<<DESTN_FOLDER_local.c_str()<<"\n";
			mkdir((DESTN_FOLDER_local+"lattices/").c_str(),0777);
			cout<<"made folder = "<<(DESTN_FOLDER_local+"lattices/").c_str()<<"\n";
			
			# pragma omp critical
			{
				fnames.open(DESTN_FOLDER+"pdetails.txt",ios::out|ios::app);
				fnames<<"LatID_"<<LatID<<","<<seed_names[m]<<","<<lattice_seed<<"\n";
				fnames.close();
			}
			// base filename from the seed_lat
			fname="LId_"+to_string(LatID)+"_L_"+to_string(current_lattice.get_lattice_len())+
												"_T_"+to_string(round(current_lattice.get_temparature()*1000)/1000).substr(0,5)+
												"_J_"+to_string(round(current_lattice.get_J()*1000)/1000).substr(0,5)+
												"_H_"+current_lattice.get_H_as_string()+"O"+to_string(current_lattice.get_spin_dim());
												
			// registering the lattice to be simulated
			//# pragma omp critical
			{
				fnames.open(DESTN_FOLDER_local+"lat_names.txt",ios::out|ios::app);
				fnames<<fname<<"\n";
				fnames.close();
			}
		
			
			cout<<"setup done for : "<<LatID<<"  L "<<current_lattice.get_lattice_len()<<", T "<<current_lattice.get_temparature();
			cout<<", J "<<current_lattice.get_J()<<", H = "<<current_lattice.get_H_as_string()<<" starting sweeps \n";

			// bookkeeping variable updates
			int part_count=0;
			int wlf_gap=0;
			int write_gap=0;
			int lat_write_gap=0;
			int lat_count=0;

			//the progressive sweep
			
			// DOING THE MONTE CARLO SIMULATION CALLS with Wolf Steps ( with intermediate Metro. Step)
			for(int qw=0;qw<int((WLF_sweep_max-1)/WLF_sweep_step)+1;qw++)
				{
				// log & prints
//					current_lattice.display_lattice();
					Temptime_req = clock() - time_req;
					log.open("log.log",ios::out|ios::app);
					log<<"WLF_SWEEP["<<WLF_sweep_step<<"] @ dim = "<<current_lattice.get_diamension();
					log<<"  L = "<<current_lattice.get_lattice_len()<<", T "<<current_lattice.get_temparature();
					log<<", J "<<current_lattice.get_J()<<", H = "<<current_lattice.get_H_as_string()<<" [ "<<qw+1<<"/";
					log<<int((WLF_sweep_max-1)/WLF_sweep_step)<<" ]  && time :: "<<float(Temptime_req/CLOCKS_PER_SEC)<<"\n";
					log.close();
					cout<<"WLF_SWEEP["<<WLF_sweep_step<<"] @ dim = "<<current_lattice.get_diamension();
					cout<<"  L = "<<current_lattice.get_lattice_len()<<", T "<<current_lattice.get_temparature();
					cout<<", J "<<current_lattice.get_J()<<", H = "<<current_lattice.get_H_as_string()<<" [ "<<qw+1<<"/";
					cout<<int((WLF_sweep_max-1)/WLF_sweep_step)<<" ]  && time :: "<<float(Temptime_req/CLOCKS_PER_SEC)<<"\n";
				
					
//					current_lattice.wolf_steps(WLF_sweep_step);
					current_lattice.mc_steps(WLF_sweep_step);
					// bookkeeping variable updates
					 wlf_gap+=WLF_sweep_step;
					 write_gap+=WLF_sweep_step;
					 lat_write_gap+=WLF_sweep_step;
					 
					
					// DOING THE MONTE CARLO SIMULATION CALLS with Metropolis Single flip Steps
					if(wlf_gap> WLF_MSF_GAP)
					{
						wlf_gap=0;
						current_lattice.mc_steps(MC_sweep_step);
						write_gap+=MC_sweep_step;
						lat_write_gap+=MC_sweep_step;
					}
					
					if(write_gap> FILE_WRITE_GAP or qw>int((WLF_sweep_max-1)/WLF_sweep_step))
					{
						cout<<"WRITING to file \n";
						part_count++;
						now = time(0);
						dtime = ctime(&now);
						current_lattice.update_history=string(dtime);
						partfname=to_string(part_count)+"_"+fname;
						//WRITING LATTICE-DATA TO FILES
						current_lattice.writeToText(DESTN_FOLDER_local+partfname+".txt",true);
						current_lattice.writeToBinary(DESTN_FOLDER_local+partfname+".dat",true);
						//# pragma omp critical
						{
							fnames.open(DESTN_FOLDER_local+"PartHfnames_"+fname+".txt",ios::out|ios::app);
							fnames<<partfname<<"\n";
							fnames.close();
						}
						
						//clearing Main Lattice Sweep data
						current_lattice.clear_data( );

						log.open("log.log",ios::out|ios::app);
						log<<"completed "<<LatID<<"  L "<<current_lattice.get_lattice_len()<<", T "<<current_lattice.get_temparature();
						log<<", J "<<current_lattice.get_J();
						log<<", H = "<<current_lattice.get_H_as_string()<<"  [ ";
						log<<float(( clock()- time_req )/CLOCKS_PER_SEC)<<" ] \n";
						log.close();
						write_gap=0;
					}
					
				if(lat_write_gap>=LATTICE_WRITE_GAP or qw>int((WLF_sweep_max-1)/WLF_sweep_step))
				{
					lat_count++;
					cout<<"WRITING LATTICE \n";
					lattice_O_n alat=current_lattice.get_bare_copy();
					now = time(0);
					dtime = ctime(&now);
					alat.update_history=string(dtime);
					alat.set_sweep_count(current_lattice.get_sweep_count());
					alat.writeToText(DESTN_FOLDER_local+"lattices/lattice_"+partfname+".txt");
//					cout<<"\t txt to "<<DESTN_FOLDER_local+"lattices/lattice_"+partfname+".txt"<<"\n";
					alat.writeToBinary(DESTN_FOLDER_local+"lattices/lattice_"+partfname+".dat");
//					cout<<"\t txt to "<<DESTN_FOLDER_local+"lattices/lattice_"+partfname+".txt"<<"\n";
					alat.clear_lattice();
					//# pragma omp critical
					{
						fnames.open(DESTN_FOLDER_local+"lattices/latnames.txt",ios::out|ios::app);
						fnames<<"lattice_"<<partfname<<","<<string(dtime);
						fnames.close();
					}
					lat_write_gap=0;
				}
			}
		current_lattice.clear_lattice();
		}
		
		delete [] lat_ids;
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



























