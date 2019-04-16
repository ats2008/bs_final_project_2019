#include "stattools.h"

double average(double *array,int len)
{
	double x= 0;
	for(int i=0;i<len;i++)
		x+=array[i];
	return x/len;
}

double variance(double *array,int len)
{
	double x= 0;
	double y=0;
	for(int i=0;i<len;i++)
		{
			x+=array[i];
		}
		x/=len;
	for(int i=0;i<len;i++)
		{
			y+=(array[i]-x)*(array[i]-x);
		}
		y/=len;
//	cout<<y<<"\n";
	return y;
}

cpp_dec_float_100 average(cpp_dec_float_100 *array,int len)
{
	cpp_dec_float_100 x= 0;
	for(int i=0;i<len;i++)
		x+=array[i];
	return x/len;
}

cpp_dec_float_100 variance(cpp_dec_float_100 *array,int len)
{
	cpp_dec_float_100 x= 0;
	cpp_dec_float_100 y=0;
	for(int i=0;i<len;i++)
		{
			x+=array[i];
		}
		x/=len;
	for(int i=0;i<len;i++)
		{
			y+=(array[i]-x)*(array[i]-x);
		}
		y/=len;
//	cout<<y<<"\n";
	return y;
}



long double* resolve_Z(vector<lattice_O_n> data_vector,long double *Z0,int MAX_ITR,double TOLARENCE)
{
	int N=data_vector.size();
	long double t;
	
	long double *Z,*Ztemp;
	long double zmax,zmin;
	Z= new long double [N];
	Ztemp= new long double [N];
	
	double *J,*H_ext;
	vector<double*> H;
	J = new double [N];
	
	
	if(Z0==NULL)
		{
			cout<<"\n initialization inside resolve \n";
			for(int i=0;i<N;i++)
				Z[i]=1;
		}
	else
		{
			for(int i=0;i<N;i++)
				Z[i]=Z0[i];
		}
	vector<lattice_O_n>::iterator lat_it;
	vector<double>::iterator en_tot;
	vector<double *>::iterator mag_n_it;
	vector<double *>::iterator H_ext_it;
	int i=0;
	int k=0;
	
	int *count;
	count= new int [N];
	int spin_dim=data_vector[0].get_n();
	
	i=0;
	for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
	{
		J[i]=(*lat_it).get_J();
		H_ext = new double [spin_dim];
		(*lat_it).get_ext_field(H_ext); 
		H.push_back(H_ext);
		count[i]=(*lat_it).energy.size();	
		i+=1;
	}
	
	i=0;
	cout<<"for itr = "<<i << " Z : = [ "<<Z[0];
	for(int p=1;p<N;p++)
		cout<<", "<<Z[p];
	cout<<" ]\n";
	cin>>i;
	

	double eh1,eh2;
	for(i=0;i<MAX_ITR;i++)
	{
//		#pragma omp parallel for private(lat_it,en_tot,mag_n_it,eh1,eh2,t)
		for(k=0;k<N;k++)
		{
			Ztemp[k]=0;
			for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
			{
			
			
				en_tot=(*lat_it).energy.begin();
				mag_n_it=(*lat_it).magnetization_O_n.begin();
				for(  ;en_tot!=(*lat_it).energy.end() and mag_n_it!=(*lat_it).magnetization_O_n.end() ;en_tot++,mag_n_it++)
				{
					t=0;
					for(int j=0;j<N;j++)
					{
						eh2=0;eh1=0;
						for(int i=0;i<spin_dim;i++)
						{
							eh1+=H[j][i]*(*mag_n_it)[i];
							eh2+=H[k][i]*(*mag_n_it)[i];
						}
						t+=(count[j]*exp((J[k]-J[j])*(*en_tot+ eh1)/J[j] + (eh1-eh2) )/Z[j]);
//						cout<<"\t\t\t\t"<<(J[k]-J[j])*(*en_tot+ eh1)/J[j]<<"   deh = "<< (eh1-eh2)<<"\n";
//						if(Z[j]!=1) cout<<"\n \t\tat z=0 T = "<<t<<"\n";
					}
					Ztemp[k]+=(1/t);
				}
//				#pragma omp critical
//				cout<<" Ztemp "<<k<<" = "<<Ztemp[k]<<" t = "<<t<<"\n";
			}
		}
		t=0;
		zmin=Ztemp[0];
		zmax=Ztemp[0];
		for(k=0;k<N;k++)
		{
			zmin=fmin(zmin,Z[k]);
			zmax=fmax(zmax,Z[k]);
		}

		zmin=sqrt(zmin*zmax);
		for(k=0;k<N;k++)
			Ztemp[k]/=zmin;

		for(k=0;k<N;k++)
		{
			t+=fabs((Z[k]-Ztemp[k])/Z[k]);
			Z[k]=Ztemp[k];

		}

		if(t<TOLARENCE)
		{
			cout<<"  tolarance : "<<t<<"  breaking";
			break;	
		}
		cout<<"for itr = "<<i+1<< " Z : = [ "<<Z[0];
		for(int p=1;p<N;p++)
			cout<<", "<<Z[p	];
		cout<<" ] tol = "<<t<<"\n";
	}
	
	delete [] count;
	delete [] J;
	for(H_ext_it=H.begin();H_ext_it!=H.end();H_ext_it++)
		delete [] *H_ext_it;
	delete [] Ztemp; 
	return Z;
}


long double get_partition_function(double desiered_J,double* desiered_H,vector<lattice_O_n> data_vector,long double * Z)
{
	if (Z==NULL)
	{
		Z=resolve_Z(data_vector);
	}
	
//	double beta_d=1/desiered_T;	
	int N = data_vector.size();
	long double Z_desired=0;
	
	double *J,*H_ext;
	vector<double*> H;
	J = new double [N];
	
	int *count;
	count= new int [N];
	int spin_dim=data_vector[0].get_n();
	
	vector<lattice_O_n>::iterator lat_it;
	int i=0;
	for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
	{
		J[i]=(*lat_it).get_J();
		H_ext = new double [spin_dim];
		(*lat_it).get_ext_field(H_ext); 
		H.push_back(H_ext);
		count[i]=(*lat_it).energy.size();	
		i+=1;
	}
	
	long double t,eh1,eh2;
	vector<double>::iterator en_tot;
	vector<double *>::iterator mag_n_it;
	vector<double *>::iterator H_ext_it;
	
	lat_it=data_vector.begin();
	H_ext_it=H.begin();
	for(;lat_it!=data_vector.end() and H_ext_it!=H.end();lat_it++,H_ext_it++)
	{
		en_tot=(*lat_it).energy.begin();
		mag_n_it=(*lat_it).magnetization_O_n.begin();
		for(  ;en_tot!=(*lat_it).energy.end() and mag_n_it!=(*lat_it).magnetization_O_n.end() ;en_tot++,mag_n_it++)
		{
			t=0;
			for(int j=0;j<N;j++)
			{	
				eh2=0;eh1=0;
				for(int i=0;i<spin_dim;i++)
				{
					eh1+=(*H_ext_it)[i]*(*mag_n_it)[i];
					eh2+=desiered_H[i]*(*mag_n_it)[i];
				}
				t+=(count[j]*exp((desiered_J-J[j])*(*en_tot+ eh1)/J[j] + (eh1-eh2) )/Z[j]);
			}
			
			Z_desired+=(1/t);
		}
	}
	
	delete [] J;
	delete [] count;
	for(H_ext_it=H.begin();H_ext_it!=H.end();H_ext_it++)
		delete [] *H_ext_it;
	
	return Z_desired;
}

double get_expectation_value(double desiered_J,double* desiered_H,double (*func)(double E,double *M),vector<lattice_O_n> data_vector,long double * Z)
{
	if (Z==NULL)
	{
		Z=resolve_Z(data_vector);
	}
	
	int N = data_vector.size();
	long double Z_desired=get_partition_function(desiered_J,desiered_H,data_vector,Z);
	
//	
	double *J,*H_ext;
	vector<double*> H;
	J = new double [N];
	
	int *count;
	count= new int [N];
	int spin_dim=data_vector[0].get_n();
	
	vector<lattice_O_n>::iterator lat_it;
	int i=0;
	for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
	{
		J[i]=(*lat_it).get_J();
		H_ext = new double [spin_dim];
		(*lat_it).get_ext_field(H_ext); 
		H.push_back(H_ext);
		count[i]=(*lat_it).energy.size();	
		i+=1;
	}


//
	long double t,expetation=0;
	double eh1,eh2;
	
	vector<double>::iterator en_tot;
	vector<double *>::iterator mag_n_it;
	vector<double *>::iterator H_ext_it;	
	
	lat_it=data_vector.begin();
	H_ext_it=H.begin();
	
	for(;lat_it!=data_vector.end() and H_ext_it!=H.end();lat_it++,H_ext_it++)
	{

		for(en_tot=(*lat_it).energy.begin(),
			mag_n_it=(*lat_it).magnetization_O_n.begin();
			en_tot!=(*lat_it).energy.end() and mag_n_it!=(*lat_it).magnetization_O_n.end();
			en_tot++, mag_n_it++)
				{
					t=0;
					for(int j=0;j<N;j++)
					{	
						eh2=0;eh1=0;
						for(int i=0;i<spin_dim;i++)
						{
							eh1+=(*H_ext_it)[i]*(*mag_n_it)[i];
							eh2+=desiered_H[i]*(*mag_n_it)[i];
						}
						t+=(count[j]*exp((desiered_J-J[j])*(*en_tot+ eh1)/J[j] + (eh1-eh2) )*(Z_desired/Z[j]));
//						cout<<"\t\t\t\t\t"<<(desiered_J-J[j])*(*en_tot+ eh1)/J[j]<<(eh1-eh2) <<"\n";
					}
			
					expetation+=func(*en_tot,*mag_n_it)/t;
				}
	}
	
	delete [] J;
	delete [] count;
	for(H_ext_it=H.begin();H_ext_it!=H.end();H_ext_it++)
		delete [] *H_ext_it;
	
	return double(expetation);
}



cpp_dec_float_100 * resolve_Z_large(vector<lattice_O_n> data_vector,cpp_dec_float_100 *Z0,string Z_file,int MAX_ITR,double TOLARENCE)
{
	int N=data_vector.size();
	cpp_dec_float_100 t;
	
	cpp_dec_float_100 *Z,*Ztemp;
	cpp_dec_float_100 zmax,zmin;
	Z= new cpp_dec_float_100 [N];
	Ztemp= new cpp_dec_float_100 [N];
	cout<<"saving temp file to "<<Z_file<<"\n";
	double *J,*H_ext;
	vector<double*> H;
	J = new double [N];
	
	
	if(Z0==NULL)
		{
			cout<<"\n initialization inside resolve \n";
			for(int i=0;i<N;i++)
				Z[i]=1.0;
		}
	else
		{
			for(int i=0;i<N;i++)
				Z[i]=Z0[i];
		}
	vector<lattice_O_n>::iterator lat_it;
	vector<double>::iterator en_tot;
	vector<double *>::iterator mag_n_it;
	vector<double *>::iterator H_ext_it;
	int i=0;
	int k=0;
	
	int *count;
	count= new int [N];
	int spin_dim=data_vector[0].get_n();
	
	i=0;
	for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
	{
		J[i]=(*lat_it).get_J();
		H_ext = new double [spin_dim];
		(*lat_it).get_ext_field(H_ext); 
		H.push_back(H_ext);
		count[i]=(*lat_it).energy.size();	
		i+=1;
	}
	
	i=0;
	cout<<"for itr = "<<i << " Z : = [ "<<Z[0];
	for(int p=1;p<N;p++)
		cout<<", "<<Z[p	];
	cout<<" ]\n";
	if (i==0) return Z;
	double eh1,eh2;
	for(i=0;i<MAX_ITR;i++)
	{
		#pragma omp parallel for private(lat_it,en_tot,mag_n_it,eh1,eh2,t)
		for(k=0;k<N;k++)
		{
			Ztemp[k]=0.0;
			for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
			{
			
			
				en_tot=(*lat_it).energy.begin();
				mag_n_it=(*lat_it).magnetization_O_n.begin();
				for(  ;en_tot!=(*lat_it).energy.end() and mag_n_it!=(*lat_it).magnetization_O_n.end() ;en_tot++,mag_n_it++)
				{
					t=0;
					for(int j=0;j<N;j++)
					{
						eh2=0;eh1=0;
						for(int i=0;i<spin_dim;i++)
						{
							eh1+=H[j][i]*(*mag_n_it)[i];
							eh2+=H[k][i]*(*mag_n_it)[i];
						}
						t+=(count[j]*exp(cpp_dec_float_100((J[k]-J[j])*(*en_tot+ eh1)/J[j] + (eh1-eh2)))/Z[j]);
//						cout<<"\t\t\t\t"<<(J[k]-J[j])*(*en_tot+ eh1)/J[j]<<"   deh = "<< (eh1-eh2)<<"t = "<<t<<"\n";
//						if(Z[j]!=1) cout<<"\n \t\tat z=0 T = "<<t<<"\n";
					}
					Ztemp[k]+=(1/t);
				}
//				#pragma omp critical
//				cout<<"\t\tfor itr = "<<i<<" Z[ "<<k<<"] = "<<Ztemp[k]<<"\n";
			}
		}
		t=0;
		zmin=Ztemp[0];
		zmax=Ztemp[0];
		for(k=0;k<N;k++)
		{
			zmin=fmin(zmin,Z[k]);
			zmax=fmax(zmax,Z[k]);
		}

		zmin=sqrt(zmin*zmax);
		for(k=0;k<N;k++)
			Ztemp[k]/=zmin;

		for(k=0;k<N;k++)
		{
			t+=fabs((Z[k]-Ztemp[k])/Z[k]);
			Z[k]=Ztemp[k];

		}

		if(t<TOLARENCE)
		{
			cout<<"  tolarance : "<<t<<"  breaking";
			break;	
		}
		cout<<"for itr = "<<i+1<< " Z : = [ "<<Z[0];
		for(int p=1;p<N;p++)
			cout<<", "<<Z[p	];
		cout<<" ] tol = "<<t<<"\n";
		if(Z_file!="")
		{
			save_Z(Z,N,Z_file);
			cout<<"saving temp file to "<<Z_file<<"\n";
		}
	}
	
	delete [] count;
	delete [] J;
	for(H_ext_it=H.begin();H_ext_it!=H.end();H_ext_it++)
		delete [] *H_ext_it;
	delete [] Ztemp; 
	return Z;
}

cpp_dec_float_100 get_partition_function_large(double desiered_J,double* desiered_H,vector<lattice_O_n> data_vector,cpp_dec_float_100 * Z)
{
	if (Z==NULL)
	{
		cout<<"\n recived no Z at get_partition_function_large(), so resolving it now !!\n";
		Z=resolve_Z_large(data_vector);
	}
	
//	double beta_d=1/desiered_T;	
	int N = data_vector.size();
	cpp_dec_float_100 Z_desired=0;
	
	double *J,*H_ext;
	vector<double*> H;
	J = new double [N];
	
	int *count;
	count= new int [N];
	int spin_dim=data_vector[0].get_n();
	
	vector<lattice_O_n>::iterator lat_it;
	int i=0;
	for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
	{
		J[i]=(*lat_it).get_J();
		H_ext = new double [spin_dim];
		(*lat_it).get_ext_field(H_ext); 
		H.push_back(H_ext);
		count[i]=(*lat_it).energy.size();	
		i+=1;
	}
	
	cpp_dec_float_100 t;
	double eh1,eh2;
	vector<double>::iterator en_tot;
	vector<double *>::iterator mag_n_it;
	vector<double *>::iterator H_ext_it;
	
	lat_it=data_vector.begin();
	H_ext_it=H.begin();
	for(;lat_it!=data_vector.end() and H_ext_it!=H.end();lat_it++,H_ext_it++)
	{
		en_tot=(*lat_it).energy.begin();
		mag_n_it=(*lat_it).magnetization_O_n.begin();
		for(  ;en_tot!=(*lat_it).energy.end() and mag_n_it!=(*lat_it).magnetization_O_n.end() ;en_tot++,mag_n_it++)
		{
			t=0;
			for(int j=0;j<N;j++)
			{	
				eh2=0;eh1=0;
				for(int i=0;i<spin_dim;i++)
				{
					eh1+=(*H_ext_it)[i]*(*mag_n_it)[i];
					eh2+=desiered_H[i]*(*mag_n_it)[i];
				}
				t+=(count[j]*exp(cpp_dec_float_100((desiered_J-J[j])*(*en_tot+ eh1)/J[j] + (eh1-eh2) ))/Z[j]);
			}
			
			Z_desired+=(1/t);
		}
	}
	
	delete [] J;
	delete [] count;
	for(H_ext_it=H.begin();H_ext_it!=H.end();H_ext_it++)
		delete [] *H_ext_it;
	
	return Z_desired;
}

cpp_dec_float_100 get_expectation_value_large(double desiered_J,double* desiered_H,double (*func)(double E,double *M),vector<lattice_O_n> data_vector,
			cpp_dec_float_100 * Z)
{
	if (Z==NULL)
	{
		cout<<"\n recived no Z at get_expectation_value_large(), so resolving it now !!\n";
		Z=resolve_Z_large(data_vector);
	}
	
	int N = data_vector.size();
	cpp_dec_float_100 Z_desired=get_partition_function_large(desiered_J,desiered_H,data_vector,Z);
	
//	
	double *J,*H_ext;
	vector<double*> H;
	J = new double [N];
	
	int *count;
	count= new int [N];
	int spin_dim=data_vector[0].get_n();
	
	vector<lattice_O_n>::iterator lat_it;
	int i=0;
	for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
	{
		J[i]=(*lat_it).get_J();
		H_ext = new double [spin_dim];
		(*lat_it).get_ext_field(H_ext); 
		H.push_back(H_ext);
		count[i]=(*lat_it).energy.size();	
		i+=1;
	}


//
	cpp_dec_float_100 t,expetation=0;
	double eh1,eh2;
	
	vector<double>::iterator en_tot;
	vector<double *>::iterator mag_n_it;
	vector<double *>::iterator H_ext_it;	
	
	lat_it=data_vector.begin();
	H_ext_it=H.begin();
	
	for(;lat_it!=data_vector.end() and H_ext_it!=H.end();lat_it++,H_ext_it++)
	{

		for(en_tot=(*lat_it).energy.begin(),
			mag_n_it=(*lat_it).magnetization_O_n.begin();
			en_tot!=(*lat_it).energy.end() and mag_n_it!=(*lat_it).magnetization_O_n.end();
			en_tot++, mag_n_it++)
				{
					t=0;
					for(int j=0;j<N;j++)
					{	
						eh2=0;eh1=0;
						for(int i=0;i<spin_dim;i++)
						{
							eh1+=(*H_ext_it)[i]*(*mag_n_it)[i];
							eh2+=desiered_H[i]*(*mag_n_it)[i];
						}
						t+=(count[j]*exp(cpp_dec_float_100((desiered_J-J[j])*(*en_tot+ eh1)/J[j] + (eh1-eh2) ))*(Z_desired/Z[j]));
//						cout<<"\t\t\t\t\t"<<(desiered_J-J[j])*(*en_tot+ eh1)/J[j]<<(eh1-eh2) <<" and  t = "<<t<<"\n";
					}
			
					expetation+=cpp_dec_float_100(func(*en_tot,*mag_n_it))/t;
				}
//	cout<<"\t\t\t\t expetation  = "<<expetation <<" and  t = "<<t<<"\n";
	}
	
	delete [] J;
	delete [] count;
	for(H_ext_it=H.begin();H_ext_it!=H.end();H_ext_it++)
		delete [] *H_ext_it;
	
	return expetation;
}


void save_Z(cpp_dec_float_100 *Z,int N,string fname)
{
	fstream f;
	f.open(fname.c_str(),ios::out | ios::binary);
	f.write((char *)(&N),sizeof(int));
	for(int i=0;i<N;i++)
	{
		f.write((char *)(&Z[i]),sizeof(cpp_dec_float_100));
	}
	f.close();
}

cpp_dec_float_100 * read_Z(int *N,string fname)
{
	fstream f;
	cpp_dec_float_100 *Z;
	f.open(fname.c_str(),ios::in | ios::binary);
	if(!f.is_open())
	{
		cout<<"\n No exixiting file to read Z retuning NULL !! \n";
		Z=NULL;
		*N=0;
		return Z;
	}
	f.read((char *)(N),sizeof(int));
	Z=new cpp_dec_float_100 [*N];
	for(int i=0;i<*N;i++)
	{
		f.read((char *)(&Z[i]),sizeof(cpp_dec_float_100));
		cout<<i<<" -> "<<Z[i]<<"\n";
	}
	f.close();
	return Z;
}




