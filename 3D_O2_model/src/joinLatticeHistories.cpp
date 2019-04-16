#include <iostream>
#include <ctime>
#include "latticelib.h"
//#include "iolib.h"
#include "stattools.h"
#include <string>
#include <stdlib.h>
#include <sys/stat.h>

using namespace std;


string BASE_FOLDER  = "simulation_data/raw_thermalized_lattices/";
string DESTN_FOLDER = "simulation_data/compiled_thermalized_lattices/";


int main(int argc,char ** argv)
{
	if(argc>1)
		{
			BASE_FOLDER=argv[1];
		}
	else
		{
			cout<<"\n plz enter the *source_path , destn path , destn folder name\n ";
			cout<<"IF not set, destn of files = source_path/compiled/";
			return 0;
		}
	if(argc>3)
		{
			DESTN_FOLDER=argv[2];
			DESTN_FOLDER+=argv[3];
			DESTN_FOLDER+="/";
			mkdir(DESTN_FOLDER.c_str(),0777);
			mkdir((DESTN_FOLDER+"lattices/").c_str(),0777);
			mkdir((DESTN_FOLDER+"history/").c_str(),0777);

		}
	else
		{
			DESTN_FOLDER=argv[1];
			DESTN_FOLDER+="compiled/";
			mkdir(DESTN_FOLDER.c_str(),0777);
			mkdir((DESTN_FOLDER+"lattices/").c_str(),0777);
			mkdir((DESTN_FOLDER+"history/").c_str(),0777);
		}
	
	string temp_str;
	
	//Reading names of the simulation part files
	temp_str=BASE_FOLDER+"lat_names.txt";
	fstream file(temp_str.c_str(),ios::in);
	if(!file)
	{
		cout<<temp_str<<"\n ERROR OPENING FILE \n";
	}
	lattice_O_n alat,blat;
	vector<string> fnames,partFnames;
	string l="";
	while( getline (file,l,'\n') )
	{
		fnames.push_back(l);
		cout<<"adding partfile name : "<<l<<"\n";
	}
	file.close();
	
	// Resampling  the main files
	int i=0;
	vector<string>::iterator name = fnames.begin();
	vector<string>::iterator partname;
	vector<double>::iterator eng_itr;;
	vector<double>::iterator mag_itr;
	vector<double * >::iterator magn_itr;
	double * mg_n;

	fstream ofile,latnames;
	ofile.open(DESTN_FOLDER+"history/Hfnames.txt",ios::out);
	latnames.open(DESTN_FOLDER+"lattices/latnames.txt",ios::out);
	cout<<"_______________________________________________\n";
	cout<<"             PROCESSING PART FILES\n";
	cout<<"_______________________________________________\n";
	i=0;
	int _n;
	for(;name!=fnames.end();name++)
	{
		i++;
//		if (i>1) break;
		cout<<i<<" / "<<fnames.size()<<"\n";
		cout<<"READING PartHfnames file for \t:"<< i <<" : "<<BASE_FOLDER+"PartHfnames_"+*name<<"\n";
		file.open(BASE_FOLDER+"PartHfnames_"+*name+".txt");
		
		if(!file)
		{
			cout<<BASE_FOLDER+"PartHfnames_"+*name<<"\n ERROR OPENING FILE \n";
			continue;
		}
		while( getline (file,l,'\n') )
		{
			partFnames.push_back(l);
		}
		file.close();
		cout<<"lengh of part file list : "<<partFnames.size()<<"\n";
		int j=0;
		partname=partFnames.begin();
		alat.readFromBinary(BASE_FOLDER+"history_"+*partname+".dat");
//		alat.display_lattice_detailed();
		partname++;
		j++;
		int count=alat.energy.size();
//		cout<<"\n adding the additional lattices and now count = "<<count<<"\n";
		for(;partname!=partFnames.end();partname++)
		{
			j++;
//			if(j>2) break;
			cout<<j<<" / "<<partFnames.size()<<".."<<"history_"+*partname;
			blat.readFromBinary(BASE_FOLDER+"history_"+*partname+".dat");
			eng_itr=blat.energy.begin();
			mag_itr=blat.magnetization.begin();
			magn_itr=blat.magnetization_O_n.begin();
			_n=blat.get_spin_dim();
			int size =blat.energy.size();
			cout<<" size of current file  = "<<size<<" count = "<<count<< ",  -> total size : "<<alat.magnetization.size()<<" \n";
			eng_itr++;
			mag_itr++;
			magn_itr++;
			while(size-1 >0)
			{	
				count++;
//				cout<<count<<",";
				alat.energy.push_back(*eng_itr);
				alat.magnetization.push_back(*mag_itr);
				mg_n =new double [_n];
				for(int js=0;js<_n;js++)
					mg_n[js]=(*magn_itr)[js];
				alat.magnetization_O_n.push_back(mg_n);
				size--;
				eng_itr++;
				mag_itr++;
				magn_itr++;
			}
//			cout<<"energy and mags set \n";
			alat.set_sweep_count(count);
			alat.update_history=string(blat.update_history);
//			cout<<"clearing blat now \n";
			blat.clear_lattice();
//			cout<<"cleared blat now \n";
		}
		cout<<"reading for 'lattice save' : "<<BASE_FOLDER+"history_"+partFnames.back()+".dat"<<"\n";
		blat.readFromBinary(BASE_FOLDER+"history_"+partFnames.back()+".dat");
		blat.clear_data();
		latnames<<"lattice_"<<*name<<",";
		latnames<<blat.update_history;
		blat.writeToBinary(DESTN_FOLDER+"lattices/lattice_"+*name+".dat");
		blat.writeToText(DESTN_FOLDER+"lattices/lattice_"+*name+".txt");

		partFnames.clear();
//		alat.display_lattice_detailed();
		cout<<"Done with the joining of the lattice histories  @ "<<alat.get_sweep_count()<<"\n";
		alat.writeToBinary(DESTN_FOLDER+"history/"+*name+".dat");
		alat.writeToText(DESTN_FOLDER+"history/"+*name+".txt");
		cout<<"Writing the resampled file at  : "<<DESTN_FOLDER+"history/"+*name+".***"<<"\n\n";
		ofile<<*name+"\n";
//			cout<<"\n final clear\n";
		blat.clear_lattice();
//		cout<<"\n final B clear\n";
		alat.clear_lattice();
//		cout<<"\n final A clear\n";
	}
	cout<<"\n DONE AND DONE !! \n";
	ofile.close();
	latnames.close();
return 0;
}
