
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
//					if(m<5000)
//						continue;
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

long double get_partition_function(double desiered_J,double* desiered_H,vector<lattice> data_vector,long double * Z)
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
	dim=data_vector[0].get_diamension();
	
	vector<lattice>::iterator lat_it;
	int i=0;
	for(lat_it=data_vector.begin();lat_it!=data_vector.end();lat_it++)
	{
		J[i]=(*lat_it).get_J();
		H_ext = new double [dim];
		(*lat_it).get_ext_field(H_ext); 
		H.push_back(H_ext);
		count[i]=(*lat_it).energy.size();	
		i+=1;
	}
	
	long double t,eh1,eh2;
	vector<double>::iterator en_tot;
	vector<double *>::iterator mag_n_it;
	
	lat_it=data_vector.begin();
	vector<double *>::iterator H_ext_it;
	H_ext_it=H.begin();
	for(;lat_it!=data_vector.end() and H_ext_it!=H.end();lat_it++,H_ext_it++)
	{
		en_it=(*lat_it).energy.begin();
		mag_n_it=(*lat_it).magnetization_O_n.begin();
		for(  ;en_it!=(*lat_it).energy.end() and mgn_it!=(*lat_it).energy.end ;en_it++,mgn_it++)
		{
			t=0;
			for(int j=0;j<N;j++)
			{	
				eh2=0;eh1=0;
				for(int i=0;i<dim;i++)
				{
					eh1+=(*H_ext_it)[j]*(*mgn_it)[j];
					eh2+=(desiered_H[j]*(*mgn_it)[j];
				}
				t+=(count[j]*exp((J_d-J[j])*(*en_it+ eh1)/J[j] + (eh1-eh2) )/Z[j]);
			}
			
			Z_desired+=(1/t);
		}
	}
	
	delete [] J;
	delete [] count;
	for(H_ext_it=H.begin();H_ext_it!=H.end();H_ext_it++)
		delete [] H_ext_it;
	
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
		int indip_sample_period=(1.3*(tau+tau_err)+1);
		cout<<"REsampling at interval of "<<indip_sample_period<<"  \n";
		if( float(tau_err/tau)> tau_err_tolarance) 
			{
				cout<<"\n\n Tau Have err > Tol !! "<<tau_err*100/tau<<" %"<<"\n\n";
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

lattice resample_thermalized_xyModel(lattice alat,int tau,double tau_err,int thermalization_skip,float tau_err_tolarance)
{
	int indip_sample_period=(1.3*(tau+tau_err)+1);
	cout<<"Resampling at interval of "<<indip_sample_period<<"  \n";
	if( float(tau_err/tau)> tau_err_tolarance) 
		{
			cout<<"\n Tau Have err > Tol !! "<<tau_err*100/tau<<" %"<<"\n";
		}
	lattice_write parameters;
	parameters=alat.get_parameters_xyModel();
	lattice blattice;
	blattice.setup_xy_model(parameters._dimension,parameters._L,parameters._site_count,parameters._T,parameters._J,parameters._h4,
												parameters._h8,parameters._initalization,'X');
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
	blattice.magnetization_x.clear();
	blattice.magnetization_y.clear();
	
	vector<double>::iterator en_it;
	vector<double>::iterator mg_it;
	vector<double>::iterator mgx_it;
	vector<double>::iterator mgy_it;
	en_it = alat.energy.begin();
	mg_it = alat.magnetization.begin();
	mgx_it = alat.magnetization_x.begin();
	mgy_it = alat.magnetization_y.begin();
	int m=0;
	int i=0;
	bool Not_thermalized=true;
	while(en_it!=alat.energy.end() and mg_it!=alat.magnetization.end()
					and mgx_it!=alat.magnetization_x.end() and mgy_it!=alat.magnetization_y.end()	)
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
			mgx_it++;
			mgy_it++;
			en_it++;
			continue;
		}
		if(m==indip_sample_period)
		{
			i++;
			blattice.add_energy(*en_it);
			blattice.add_magnetization(*mg_it);
			blattice.magnetization_x.push_back(*mgx_it);
			blattice.magnetization_y.push_back(*mgy_it);
			m=0;
		}
		mg_it++;
		mgx_it++;
		mgy_it++;
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




