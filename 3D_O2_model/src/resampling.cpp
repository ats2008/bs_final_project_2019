#include <iostream>
#include <ctime>
#include "latticelib.h"
#include "stattools.h"
#include <string>
#include <stdlib.h>
#include <sys/stat.h>

using namespace std;

int main(int argc,char ** argv)
{
	string BASE_FOLDER  = "simulation_data/raw/";
	string DESTN_FOLDER = "simulation_data/resampled/";
	int thermalization_skip = 0;
	if(argc>3)
		{
			BASE_FOLDER=argv[1];
			DESTN_FOLDER=argv[2];
			DESTN_FOLDER+=argv[3];
			DESTN_FOLDER+="/";
			mkdir(DESTN_FOLDER.c_str(),0777);
		}
	else
	{
		cout<<"Please provide *source_path ,*destn_path,*destn_folder_name,thermalization_skip\n";
		cout<<"files are to be saved under destn_path/destn_folder_name \n";
		return 0;
	}
	if(argc>4)
		thermalization_skip=int(atof(argv[4]));
	
	string temp_str;
	
	//Reading names & tau of the simulation files
	temp_str=BASE_FOLDER+"resampling.config";
	fstream file(temp_str.c_str(),ios::in);
	if(!file)
	{
		cout<<temp_str<<"\n ERROR OPENING FILE \n";
	}
	lattice_O_n alat,blat;
	vector<string> fnames;
	vector<int> tau;
	vector<double> tau_err;
	string l="";
	while( getline (file,l,',') )
	{
		fnames.push_back(l);
		getline(file,l,',');
		tau.push_back(int(atof(l.c_str())));
		getline(file,l,',');
		tau_err.push_back(atof(l.c_str()));
		getline(file,l,'\n');
	}
	file.close();
	
	
	
	// Resampling  the main files
	int i=0;
	vector<string>::iterator name = fnames.begin();
	vector<int>::iterator tau_it=tau.begin();
	vector<double>::iterator tau_err_it=tau_err.begin();

	fstream ofile;
	ofile.open(DESTN_FOLDER+"Hfnames.txt",ios::out);
	cout<<"_______________________________________________\n";
	cout<<"_______________________________________________\n";
	i=0;
	for(;name!=fnames.end();name++)
	{
		i++;
		cout<<i<<" / "<<fnames.size()<<"\n";
		cout<<"READING file\t:"<< i <<" : "<<*name<<"\n";
		alat.readFromBinary(BASE_FOLDER+*name+".dat");
		cout<<"Resampling "<<i<<" with tau,tau_err : [ "<<*tau_it<<" , "<<*tau_err_it<<" ] \n";
//		blat=alat.resample_thermalized(*tau_it,*tau_err_it,thermalization_skip);
		blat=alat.resample_thermalized_averages(*tau_it,*tau_err_it,thermalization_skip);
		cout<<"new sample space has "<<blat.energy.size()<<" items \n";
		blat.writeToBinary(DESTN_FOLDER+*name+".dat");
		blat.writeToText(DESTN_FOLDER+*name+".txt");
		ofile<<*name+"\n";
		cout<<"Writing the resampled file at  : "<<DESTN_FOLDER+*name+".***"<<"\n";
		tau_it++;
		tau_err_it++;
	}
	ofile.close();
	
	
return 0;
}
