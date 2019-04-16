
#include "stattools.h"
long double _logadd(long double l1,long double l2)
{
	long double lM= fmax(l1,l2);
	long double lm= fmin(l1,l2);
	return lM + log1p(exp(lm-lM));
}

long double* resolve_Z(vector<lattice> data_vector,long double *Z0,int MAX_ITR,double TOLARENCE)
{
	int N=data_vector.size();
	long double t;
	
	long double *Z,*Ztemp;
	long double zmax,zmin;
	Z= new long double [N];
	Ztemp= new long double [N];
	
	double *beta;
	beta= new double [N];
	
	int *count;
	count= new int [N];
	
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
	vector<lattice>::iterator lat_it;
	vector<double>::iterator dub_it;
	int i=0;
	int j=0;
	int k=0;
	for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
	{
		beta[i]=1/(*lat_it).get_temparature();
		count[i]=(*lat_it).energy.size();	
		i+=1;
	}
	i=0;
	cout<<"for itr = "<<i << " Z : = [ "<<Z[0];
	for(int p=1;p<N;p++)
		cout<<", "<<Z[p	];
	cout<<" ]\n";

	for(i=0;i<MAX_ITR;i++)
	{
		for(k=0;k<N;k++)
		{
			Ztemp[k]=0;
			for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
			{
				int m=0;
				for(dub_it=(*lat_it).energy.begin();dub_it!=(*lat_it).energy.end();dub_it++)
				{
					m++;
					if(m<5000)
						continue;
					t=0;
					for(j=0;j<N;j++)
						t+=(count[j]*exp((beta[k]-beta[j])*(*dub_it))/Z[j]);
					Ztemp[k]+=(1/t);
				}	
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
	delete [] beta;
	delete [] Ztemp; 
	return Z;
}

long double get_partition_function(double desiered_T,vector<lattice> data_vector,long double * Z)
{
	if (Z==NULL)
	{
		Z=resolve_Z(data_vector);
	}
	
	double beta_d=1/desiered_T;	
	int N = data_vector.size();
	long double Z_desired=0;
	
	double *beta;
	beta= new double [N];
	
	int *count;
	count= new int [N];
	
	vector<lattice>::iterator lat_it;
	vector<double>::iterator dub_it;
	int i=0;
	for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
	{
		beta[i]=1/(*lat_it).get_temparature();
		count[i]=(*lat_it).energy.size();	
		i+=1;
	}
	
	long double t;
	
	for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
	{
		for(dub_it=(*lat_it).energy.begin();dub_it!=(*lat_it).energy.end();dub_it++)
		{
			t=0;
			for(int j=0;j<N;j++)
				t+=(count[j]*exp((beta_d-beta[j])*(*dub_it))/Z[j]);
			Z_desired+=(1/t);
		}
	}
	return Z_desired;
}

double get_expectation_value(double desiered_T,double (*func)(double E,double M),vector<lattice> data_vector,long double * Z)
{
	if (Z==NULL)
	{
		Z=resolve_Z(data_vector);
	}
	
	double beta_d=1/desiered_T;	
	int N = data_vector.size();
	long double Z_desired=get_partition_function(desiered_T,data_vector,Z);
	
	double *beta;
	beta= new double [N];
	
	int *count;
	count= new int [N];
	
	vector<lattice>::iterator lat_it;
	vector<double>::iterator dub_it_energy,dub_it_mag;
	int i=0;
	for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
	{
		beta[i]=1/(*lat_it).get_temparature();
		count[i]=(*lat_it).energy.size();	
		i+=1;
	}
	
	long double t,expetation=0;
	for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
	{
		for(dub_it_energy=(*lat_it).energy.begin(),
			dub_it_mag=(*lat_it).magnetization.begin();
			dub_it_energy!=(*lat_it).energy.end();
			dub_it_energy++, dub_it_mag++)
				{
					t=0;
					for(int j=0;j<N;j++)
						t+=(count[j]*exp((beta_d-beta[j])*(*dub_it_energy))*(Z_desired/Z[j]));
					expetation+=func(*dub_it_energy,*dub_it_mag)/t;
				}
	}
	return double(expetation);
}

lattice resample_thermalized(lattice alat,int tau,double tau_err,int thermalization_skip,float tau_err_tolarance)
{
		int indip_sample_period=(2*(tau+tau_err)+1);
		if( float(tau_err/tau)> tau_err_tolarance) 
			{
				cout<<"\n Tau Have err > Tol !! "<<tau_err*100/tau<<" %"<<"\n";
			}
		lattice_write parameters;
		parameters=alat.get_parameters();
		lattice blattice(
			parameters._dimension,
			parameters._L,
			parameters._site_count,
			parameters._T,
			parameters._J,
			parameters._q,
			parameters._initalization
			);
			
	blattice.update_history=parameters.update_history;
	blattice.seed_generator(parameters._urd_seed,parameters._initalization_seed);
	site * state,*state_old;
	state=blattice.get_state();
	state_old=alat.get_state();
	for(int i=0;i<parameters._site_count;i++)
		{
			state[i]=state_old[i];
		}
	blattice.set_neighboures();
	blattice.energy.clear();
	blattice.magnetization.clear();
	
	vector<double>::iterator en_it;
	vector<double>::iterator mg_it;
	en_it = alat.energy.begin();
	mg_it = alat.magnetization.begin();
	int m=0;
	int i=0;
	bool Not_thermalized=true;
	while(en_it!=alat.energy.end() and mg_it!=alat.magnetization.end())
	{	
		m++;
		if(Not_thermalized)
		{
			if(m>=thermalization_skip)
			{
				Not_thermalized=false;
				m=0;
			}
			mg_it++;
			en_it++;
			continue;
		}
		if(m==indip_sample_period)
		{
			i++;
			blattice.add_energy(*en_it);
			blattice.add_magnetization(*mg_it);
			m=0;
		}
		mg_it++;
		en_it++;
			
	}
	blattice.set_sweep_count(i-1);
	return blattice;
}

lattice resample_end(lattice alat,int count)
{
		lattice_write parameters;
		parameters=alat.get_parameters();
		lattice blattice(
			parameters._dimension,
			parameters._L,
			parameters._site_count,
			parameters._T,
			parameters._J,
			parameters._q,
			parameters._initalization
			);
			
	blattice.update_history=parameters.update_history;
	blattice.seed_generator(parameters._urd_seed,parameters._initalization_seed);
	site * state,*state_old;
	state=blattice.get_state();
	state_old=alat.get_state();
	for(int i=0;i<parameters._site_count;i++)
		{
			state[i]=state_old[i];
		}
	blattice.set_neighboures();
	blattice.energy.clear();
	blattice.magnetization.clear();
	
	vector<double>::iterator en_it;
	vector<double>::iterator mg_it;
	en_it = alat.energy.begin();
	mg_it = alat.magnetization.begin();
	int m=0;
	int i=0;
	bool Not_end_found=true;
	int skip_count=alat.energy.size()-1-count;
	while(en_it!=alat.energy.end() and mg_it!=alat.magnetization.end())
	{	
		if(Not_end_found)
		{
			m++;
			if(m>skip_count)
			{
				Not_end_found=false;
			}
			mg_it++;
			en_it++;
			continue;
		}
		i++;
		blattice.add_energy(*en_it);
		blattice.add_magnetization(*mg_it);
		mg_it++;
		en_it++;
	}
	blattice.set_sweep_count(i-1);
	return blattice;
}




