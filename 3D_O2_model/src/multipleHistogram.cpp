#include <iostream>
#include <ctime>
#include "latticelib.h"
#include "stattools.h"
#include <string>
#include <boost/multiprecision/cpp_dec_float.hpp>

#define BSTRAP_COUNT 3

using namespace std;

double e1(double e,double * M)
	{
		return e;
	}

double e2(double e,double* M)
	{
		return e*e;
	}
double m1(double e,double* M)
	{
		return M[0];
	}

double m2(double e,double *M)
	{
		return M[0]*M[0];
	}
double m4(double e,double* M)
	{
		return M[0]*M[0]*M[0]*M[0];
	}

int main(int argc, char** argv) 
{
	string out_file,BASE_FOLDER,DSTN_FOLDER,Z_file;
	if(argc>1)
	{
		BASE_FOLDER=argv[1];
		mkdir((BASE_FOLDER+"analysis").c_str(),0777);
		DSTN_FOLDER=BASE_FOLDER+"analysis/";
		out_file=DSTN_FOLDER+"mh_data";
		Z_file=BASE_FOLDER+"analysis/partionFn.dat";
	}
	else
	{
		cout<<"\nplease enter *BASE_FOLDER,**DESTN_FOLDER,**output_filename\n";
		return 0;
	}
	if(argc>2)
	{
		out_file=argv[2];
	}
	if(argc>3)
	{
		DSTN_FOLDER=argv[3];
		out_file=DSTN_FOLDER+argv[2];
	}
	string temp_str;
	temp_str=BASE_FOLDER+"Hfnames.txt";
	fstream file(temp_str.c_str(),ios::in);
	
	if(!file.is_open())
	{
		cout<<"\n ERROR HAPPEND !! FILE DOES NOT EXIST!! \n";
		cout<<"fname : "<<temp_str;
		exit(0);
	}
	
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
	
	lattice_O_n *alat;
	vector<lattice_O_n> data_vector;
	vector<lattice_O_n> BS_data_vector;
	vector<lattice_O_n>::iterator lat_it;

	int i=0;
	vector<string>::iterator name = fnames.begin();
	file<<"#read_files\n";
	for(;name!=fnames.end();name++)
	{	
		i++;
		cout<<"READING "<< i <<"th file : "<<*name<<".dat ";
		alat=new lattice_O_n;
		alat->readFromBinary(BASE_FOLDER+*name+".dat");
		data_vector.push_back(*alat);
		cout<<" data len = "<<(alat->energy).size()<<"\n";
		file<<BASE_FOLDER+*name+"R"<<"\n";
//		alat->display_lattice();
	}
	
	cout<<"\n\n ON THE VECTOR NOW \n\n";
	int N=data_vector.size();
	cpp_dec_float_100 *Z=NULL;
//	Z= new cpp_dec_float_100 [N];
	int n;
	Z=read_Z(&n,Z_file);
	if(Z)
	{
		cout<<"\n\n loded Z : = [ "<<Z[0];
		for(int p=1;p<N;p++)
			cout<<", "<<Z[p];
		cout<<"]\n";
	}
	else
	{
		cout<<" No old file loaded \n";
	}
	
	Z= resolve_Z_large(data_vector,Z,Z_file);
	/*for(int p=0;p<N;p++)
		Z[p]=1.0;*/
		
	cout<<"\n\n FINAL Z : = [ "<<Z[0];
	for(int p=1;p<N;p++)
		cout<<", "<<Z[p];
	cout<<"]\n";
	
	cout<<"writeing the Z to : "<<Z_file<<"\n";
	save_Z(Z,N,Z_file);
	cout<<"wrote the Z\n";	
	//delete [] Zr;
	cout<<"\n__________________________________________\n";
	
	
	file<<"#Z\n";
	i=0;
	for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
		{
			cout<<" Z for t = "<<(*lat_it).get_temparature()<<"  :  "<<Z[i]<<"\n";
			file<<(*lat_it).get_temparature()<<","<<(*lat_it).get_J()<<","<<(*lat_it).get_H_as_string()<<","<<Z[i]<<"\n";
			i+=1;
		}
		
	i=0;
	file<<"\n["<<Z[0];
	for(i=1;i<N;i++)
		{
			file<<","<<Z[i];
			i+=1;
		}
	file<<"]";
	
	cout<<"\n__________________________________________\n";
	
	// Reading CONFIG FILE
	double J0=0.42,J;
	double dj=0.01;
	double Jmax=0.7;

	double h[3]={0.025,0,0};
	double dh0=10;
	double H0MAX=0.026;
	
	file.close();
	file.open("config/multipleHistogram.config",ios::in);
	if(file)
		{	
			getline (file,l,'\n') ;	getline (file,l,'\n') ; J0=atof(l.c_str());
			getline (file,l,'\n') ;	getline (file,l,'\n') ; Jmax=atof(l.c_str());
			getline (file,l,'\n') ;	getline (file,l,'\n') ; dj=atof(l.c_str());
			getline (file,l,'\n') ;	getline (file,l,'\n') ; h[0]=atof(l.c_str());
			getline (file,l,'\n') ;	getline (file,l,'\n') ; H0MAX=atof(l.c_str());
			getline (file,l,'\n') ;	getline (file,l,'\n') ; dh0=atof(l.c_str());
			file.close();
		}
	cout <<"J0, J0MAX, dj = "<<J0<<"  ,  "<<Jmax<<"  ,  "<<dj<<" = > "<<(Jmax-J0)/dj<<"\n";
	cout <<"H0, H0MAX, dh0 = "<<h[0]<<"  ,  "<<H0MAX<<"  ,  "<<dh0<<"  = > "<<(H0MAX-h[0])/dh0<<"\n";

	file.open(out_file.c_str(),ios::out);
	file<<"L,J,Hx,E,E_err,M,M_err\n";
	cout<<"\n H0MAX = "<<H0MAX<<"\n";

	cpp_dec_float_100 temp100;
	cpp_dec_float_100 *E_bs, *M_bs,E,M,E_err,M_err;
	E_bs = new cpp_dec_float_100 [BSTRAP_COUNT];
	M_bs = new cpp_dec_float_100 [BSTRAP_COUNT];
	int seed_val=0;
	int L=data_vector[0].get_lattice_len();
	while(h[0]<H0MAX )
	{
		for(J=J0;J<Jmax;J+=dj)
		{
		
			cout<<"Bootstraping for h0 = "<<h[0]<<" and J = "<<J<<"\n";
			#pragma omp parallel for private(lat_it,BS_data_vector)
			for(int i=0;i<BSTRAP_COUNT;i++)
			{
					lat_it=data_vector.begin();
					for(;lat_it!=data_vector.end();lat_it++)
					{
						#pragma omp critical
						seed_val++;
						BS_data_vector.push_back((*lat_it).get_bootstrap_copy(seed_val));
					}
					E_bs[i]=get_expectation_value_large(J,h,e1,BS_data_vector,Z);
					M_bs[i]=get_expectation_value_large(J,h,m1,BS_data_vector,Z);
	//				temp=get_expectation_value(j,h,e2,data_vector,Z);
	//				temp=get_expectation_value(j,h,m2,data_vector,Z);
	//				temp=get_expectation_value(j,h,m4,data_vector,Z);

	//				cout<<" i = "<<i<<" : E = "<<E_bs[i]<<" , M = "<<M_bs[i]<<"\n";
					lat_it=BS_data_vector.begin();
					for(;lat_it!=BS_data_vector.end();lat_it++)
					{
						(*lat_it).clear_lattice();
					}
					BS_data_vector.clear();
					
			}
			E=average(E_bs,BSTRAP_COUNT);
			E_err=sqrt(variance(E_bs,BSTRAP_COUNT));
			M=average(M_bs,BSTRAP_COUNT);
			M_err=sqrt(variance(M_bs,BSTRAP_COUNT));
			cout<<" L = "<<L<<" ,J = "<<J<<",H0 = "<<h[0]<<", E = "<<E<<", E_err = "<<E_err<<", M = "<<M<<", M_err ="<<M_err<<"\n";
			file<<L<<","<<J<<","<<h[0]<<","<<E<<","<<E_err<<","<<M<<","<<M_err<<"\n";
//			cout<<" L = "<<L<<" ,J = "<<J<<",H0 = "<<h[0]<<", E = "<<E_bs[0]<<", E_err = "<<E_err<<", M = "<<M_bs[0]<<", M_err ="<<M_err<<"\n";
//			file<<L<<","<<J<<","<<h[0]<<","<<E_bs[0]<<","<<10<<","<<M_bs[0]<<","<<0<<"\n";
		}
		h[0]+=dh0;
	}
	file.close();
	delete [] Z;
	delete [] E_bs;
	delete [] M_bs;
return 0;
}
















