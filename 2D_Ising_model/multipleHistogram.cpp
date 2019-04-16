#include <iostream>
#include <ctime>
#include "isinglib.h"
#include "iolib.h"
#include "stattools.h"
#include <string>

using namespace std;
double e1(double e,double M)
	{
		return e;
	}

double e2(double e,double M)
	{
		return e*e;
	}
double m1(double e,double M)
	{
		return M;
	}

double m2(double e,double M)
	{
		return M*M;
	}
double m4(double e,double M)
	{
		return M*M*M*M;
	}

int main(int argc, char** argv) 
{
	string out_file;
	if(argc<2)
	{
		cout<<"\n OUT FILE NAME NOT SPECIFIED .. enter ./multt.. [outfile]\n";
		return 0;
	}
	out_file=argv[1];
	string BASE_FOLDER="multi_histo_data/";
	string temp_str;
	temp_str=BASE_FOLDER+"fnames.dat";
	fstream file(temp_str.c_str(),ios::in);
	
	if(!file.is_open())
	{
		cout<<"\n ERROR HAPPEND !! FILE DOES NOT EXIST !! \n";
		cout<<"fname : "<<temp_str;
		exit(0);
	}
	
	lattice alat,blat;
	vector<string> fnames;
	string l="";
	while( getline (file,l,'\n') )
	{
		fnames.push_back(l);
		cout<<l<<"\n";
	}
	file.close();
	temp_str=out_file+".config";
	file.open(temp_str.c_str(),ios::out);
	if(!file.is_open())
	{
		cout<<"\n ERROR HAPPEND !! FILE DOES NOT EXIST !! \n";
		cout<<"fname : "<<temp_str;
		exit(0);
	}
	
	vector<lattice> data_vector;
	vector<lattice>::iterator lat_it;

	int i=0;
	vector<string>::iterator name = fnames.begin();
	file<<"#read_files\n";
	for(;name!=fnames.end();name++)
	{	
		i++;
		cout<<"READING "<< i <<"th file : "<<*name<<"R\n";
		data_vector.push_back(ReadLattticeFromBinary(BASE_FOLDER+*name+"R"));
		file<<BASE_FOLDER+*name+"R"<<"\n";
	}
	cout<<"\n\n ON THE VECTOR NOW ";
	int N=data_vector.size();

	long double *Z;
	Z = new long double [N];
	for(i=0;i<N;i++)
		Z[i]=1;
//	long double Z0[4]={ 3.0663e+19, 9005.28, 8.82182e-10, 2.73177e-20};
//	Z=resolve_Z(data_vector,Z0);
	Z=resolve_Z(data_vector);
	cout<<"\n\n FINAL Z : = [ "<<Z[0];
	for(int p=1;p<N;p++)
		cout<<", "<<Z[p];
	cout<<"]\n";
	
	i=0;
	
	cout<<"\n__________________________________________\n";
	
	
	file<<"#Z\n";
	for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
		{
			cout<<" Z for t = "<<(*lat_it).get_temparature()<<"  :  "<<Z[i]<<"\n";
			file<<(*lat_it).get_temparature()<<","<<Z[i]<<"\n";
			i+=1;
		}
	
	cout<<"\n__________________________________________\n";
	
	double t,dt,temp,TMAX;
	
	TMAX=2.55;
	t=2.0;
	dt=0.002;
	file.close();
	file.open(out_file.c_str(),ios::out);
	
	file<<"Temparature,E1,E2,M,M2,M4\n";
	while(t<TMAX)
	{
		cout<<"FOR T = "<<t<<",";
		temp=get_expectation_value(t,e1,data_vector,Z);
		cout<<" U = "<<temp;
		file<<t<<","<<temp<<",";
		temp=get_expectation_value(t,e2,data_vector,Z);
		cout<<" U2 = "<<temp;
		file<<temp<<",";
		temp=get_expectation_value(t,m1,data_vector,Z);
		cout<<" M = "<<temp;
		file<<temp<<",";
		temp=get_expectation_value(t,m2,data_vector,Z);
		cout<<" M2 = "<<temp;
		file<<temp<<",";
		temp=get_expectation_value(t,m4,data_vector,Z);
		cout<<" M4 = "<<temp<<"\n";
		file<<temp;
		file<<"\n";
		t+=dt;
	}
	file.close();
return 0;
}
