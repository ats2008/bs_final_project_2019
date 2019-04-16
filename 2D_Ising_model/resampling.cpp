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
	temp_str=BASE_FOLDER+"hist.config";
	fstream file(temp_str.c_str(),ios::in);
	int thermalization_skip = 10000;
	
	lattice alat,blat;
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
	

	int i=0;
	vector<string>::iterator name = fnames.begin();
	vector<int>::iterator tau_it=tau.begin();
	vector<double>::iterator tau_err_it=tau_err.begin();

	for(;name!=fnames.end();name++)
	{
		i++;
		cout<<"READING "<< i <<"th file : "<<*name<<"\n";
		alat = ReadLattticeFromBinary(BASE_FOLDER+*name);
		cout<<"Resampling "<<*name<<" with tau,tau_err : [ "<<*tau_it<<" , "<<*tau_err_it<<" ] ";
		blat=resample_thermalized(alat,*tau_it,*tau_err_it,thermalization_skip);
		WriteLatticeToBinary(BASE_FOLDER+*name+"R",blat);
		blat.writeToText(BASE_FOLDER+*name+"R.txt");
		cout<<"Writing the resampled file at  : "<<BASE_FOLDER+*name+"R"<<"\n";
		tau_it++;
		tau_err_it++;
	}

return 0;
}
