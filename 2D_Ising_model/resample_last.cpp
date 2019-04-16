#include <iostream>
#include <ctime>
#include "isinglib.h"
#include "iolib.h"
#include "stattools.h"
#include <string>
#include <stdlib.h>

using namespace std;
int main()
{
	string BASE_FOLDER="multi_histo_data/";
	string temp_str;
	temp_str=BASE_FOLDER+"fnames.dat";
	fstream file(temp_str.c_str(),ios::in);
	temp_str=BASE_FOLDER+"fnames.txt";
	fstream file_txt_names(temp_str.c_str(),ios::out);
	int end_read_count = 10000;
	
	lattice alat,blat;
	vector<string> fnames;
	string l="";
	while( getline (file,l,'\n') )
	{
		fnames.push_back(l);
		cout<<l<<"\n";
	}
	file.close();
	
	vector<string>::iterator name = fnames.begin();
	int i=0;
	for(;name!=fnames.end();name++)
	{	i++;
		cout<<"READING "<< i <<"th file : "<<*name<<"\n";
		alat = ReadLattticeFromBinary(BASE_FOLDER+*name);
		cout<<"Resampling ...";
		blat=resample_end(alat,end_read_count);
		cout<<"Writing the resampled file at  : "<<BASE_FOLDER+*name+"RL.*";
		//WriteLatticeToBinary(BASE_FOLDER+*name+"RL",blat);
		blat.writeToText(BASE_FOLDER+*name+"RL.txt");
		file_txt_names<<(*name+"RL.txt")<<"\n";
		cout<<"  !! DONE !!\n";
	}
	
	file_txt_names.close();
	return 0;

}
