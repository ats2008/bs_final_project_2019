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
	int iid_spacing = 0,off_set=0;
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
		cout<<"Please provide *source_path ,*destn_path,*destn_folder_name,iid_spacing,off_set\n";
		cout<<"files are to be saved under destn_path/destn_folder_name \n";
		return 0;
	}
	if(argc>4)
		iid_spacing=int(atof(argv[4]));
	if(argc>5)
		off_set=int(atof(argv[5]));
	
	string temp_str;
	
	//Reading names & tau of the simulation files
	temp_str=BASE_FOLDER+"Hfnames.txt";
	fstream file(temp_str.c_str(),ios::in);
	if(!file)
	{
		cout<<temp_str<<"\n ERROR OPENING FILE \n";
	}
	lattice_O_n alat,blat;
	string l="";
	vector<string> fnames;
	while( getline (file,l,'\n') )
	{
		fnames.push_back(l);
	}
	file.close();
	
	
	// Resampling  the main files
	int i=0;
	vector<string>::iterator name = fnames.begin();

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
		cout<<"Resampling "<<i<<" with iid_spacing :  "<<iid_spacing<<"  and offset : "<<off_set <<"\n";
		blat=alat.resample_data(iid_spacing,off_set);
		cout<<"new sample space has "<<blat.energy.size()<<" items \n";
		blat.writeToBinary(DESTN_FOLDER+*name+".dat");
		blat.writeToText(DESTN_FOLDER+*name+".txt");
		ofile<<*name+"\n";
		cout<<"Writing the resampled file at  : "<<DESTN_FOLDER+*name+".***"<<"\n";
	}
	ofile.close();
return 0;
}
