#include "isinglib.h"

#include <random>
#include <iostream>
#include <stack> 
#include <vector> 
#include <stdlib.h>
#include <set>
#include <string>
#include <fstream>

void display_site(site s)
{
	cout<<"spin = "<<s._spin<<" Neibs : "<<(s._neighbour[0])->_spin;
	for(int i=1;i<4;i++)
	{
	cout<<","<<(s._neighbour[i])->_spin;
	}
}


//default constructor
lattice::lattice()
{
	//uniform_real_distribution<double> distribution(0.0,1.0);
	_urd_seed= time(NULL);
	generator.seed(_urd_seed);
	update_history="";
}

//default destructor
void lattice::clear_lattice()
{
cout<<"\n clearing the whole lattice now !!";
	clear_data();
	for(int i=0;i<_site_count;i++)
		delete [] _state[i]._neighbour;
	delete [] _state;
	_site_count=0;
	_L=0;
	_state=NULL;
}

// constructor for lattice which initialiaze the system
lattice::lattice(int dim,int lattice_len,int N,double T_equilib,double J,int spin_count,char initialization,char model_type)
{
	_model_type =model_type;
	_state = new site[N];
	_J=J;
	_L=lattice_len;
	_site_count =N;
	_dimension=dim;
	_T=T_equilib;
	_beta=1.0/_T;
	_q=spin_count;
	//uniform_real_distribution<double> distribution(0.0,1.0);
	_initalization=initialization;
	initialize(_initalization);
	set_neighboures();
	if(_model_type=='I')
	{
		ising_standardize();
		add_new_energy(calc_energy());
	}	
	if(_model_type=='P')
	{
		add_new_energy(calc_potts_energy());
	}
	add_new_magnetization(calc_magnetization());
	set_LUT(_q,_J);
	_sweep_count=0;
	_urd_seed= time(NULL);
	generator.seed(_urd_seed);
	_wolf_add_probab=1-exp(-2*_beta*_J);
	_wolf_add_probab_potts=1-exp(-1*_beta*_J);
	update_history="";
}

// initializes the spin states
void lattice::initialize(char mode)
{

	if(mode=='Z')
	{
		_initalization_seed = -1;
		for(int i=0;i<_site_count;i++)
		{
			_state[i]._spin=0;
		}
		_sweep_count=0;
		return;
	}

	if(mode=='H')
	{
	    cout<<"\nHERE";
		_initalization_seed = time(NULL);
		for(int i=0;i<_site_count;i++)
		{
			_state[i]._spin=rand()%_q;
		}
		_sweep_count=0;
		return;
	}

}

// Function that initialize the neighbours of the site passed  .. pointers to nib at asite._neighbour
void lattice::set_neighboures()
{
	for (int i=0;i<_site_count;i++)
	{
	_state[i]._neighbour = new site* [4];
	_state[i]._neighbour[0] = _state + (i-1+_site_count)%_site_count;
	_state[i]._neighbour[1] = _state + (i+1)%_site_count;
	_state[i]._neighbour[2] = _state + (i+_L)%_site_count;
	_state[i]._neighbour[3] = _state + (i-_L+_site_count)%_site_count;
	}
}

//sets the temparature of the current state
void lattice::set_temparature(double temparature)
{
	_T=temparature;
	_beta=1.0/_T;
	_wolf_add_probab=1.0-exp(-2.0*_beta*_J);
	_wolf_add_probab_potts=1.0-exp(-1*_beta*_J);
}

void lattice::seed_generator(int seeder,int initialization_seed)
{
	if (initialization_seed!=-1)
		_initalization_seed=initialization_seed;
	_urd_seed=seeder;	
	generator.seed(seeder);
}

// Set the look uptable for fast calc
void lattice::set_LUT(int sp_count,int J)
{
	_LUT= new double [4];
	_LUT[0]=exp(-8);
	_LUT[1]=exp(-6);
	_LUT[2]=exp(-4);
	_LUT[3]=exp(-2);
}
// reading loook up table
double lattice::get_exp(double de)
{
	if(de<-7) return _LUT[0];
	if(de<-5) return _LUT[1];
	if(de<-3) return _LUT[2];
	return _LUT[3];
}



// calculates and return energy of current state
double lattice::calc_energy()
{
	double en=0;
	int nib_spin;
	for(int i=0;i<_site_count;i++)
		{
		nib_spin=_state[i]._neighbour[0]->_spin;
		nib_spin+=_state[i]._neighbour[1]->_spin;
		nib_spin+=_state[i]._neighbour[2]->_spin;
		nib_spin+=_state[i]._neighbour[3]->_spin;
		en +=_state[i]._spin*(nib_spin)*(-_J);
		}
	return en/2;
}

// calculates and return energy of current potts state
double lattice::calc_potts_energy()
{
	double en=0;
	int curr_spin,nib_spin;
	for(int i=0;i<_site_count;i++)
		{
			nib_spin=0;
			curr_spin=_state[i]._spin;
			if(curr_spin==_state[i]._neighbour[0]->_spin)
				nib_spin+=1;
			if(curr_spin==_state[i]._neighbour[1]->_spin)
				nib_spin+=1;
			if(curr_spin==_state[i]._neighbour[2]->_spin)
				nib_spin+=1;
			if(curr_spin==_state[i]._neighbour[3]->_spin)
				nib_spin+=1;
			en +=(nib_spin)*(-_J);
		}
	en+=(_L*_L)*(-_J)*2;
	return en/2;
}

// calculates and return magnetization of current state
double lattice::calc_magnetization()
{
	double mag=0;
	for(int i=0;i<_site_count;i++)
	{
		mag+=_state[i]._spin;
	}
	return mag;
}


//making the spins -1 and 1 rather than  0,1
void lattice::ising_standardize()
{
	for(int i=0;i<_site_count;i++)
		{
			if(_state[i]._spin==0)
				_state[i]._spin=-1;
		}
}

// clears the content of energy and magnetization keeping the last one intact;
void lattice::clear_data()
{
	double tem;
	if (!energy.empty())
		{	tem=energy.back();
			energy.clear();
			energy.push_back(tem);	}
	if (!magnetization.empty())
		{tem=magnetization.back();
		magnetization.clear();
		magnetization.push_back(tem);}
	update_history="";
}  

//returns the energy of the current state
double lattice::get_energy()
{
	return energy.back();
}

// adds energy of the new state using de
void lattice::add_new_energy(double de)
{
	if(energy.empty())
		energy.push_back(de);
	else
		energy.push_back(energy.back() + de);
}

// adds energy of the new state using dm
void lattice::add_new_magnetization(double dm)
{
	if(magnetization.empty())
		magnetization.push_back(dm);
	else
		magnetization.push_back(magnetization.back() + dm);
}

// adds energy of the new state
void lattice::add_energy(double e)
{
		energy.push_back(e);
}

// adds energy of the new state
void lattice::add_magnetization(double dm)
{
		magnetization.push_back(dm);
}

//returns the initialization char
char lattice::get_initialization()
{
	return _initalization;
}

//returns the sweep count
int lattice::get_sweep_count()
{
	return _sweep_count;
}

//sets the sweep count
void lattice::set_sweep_count(int swc)
{
	_sweep_count=swc;
}
//returns the site count
int lattice::get_site_count()
{
	return _site_count;
}
//returns the site array adress
site * lattice::get_state()
{
	return _state;
}
//returns the lattice length
int lattice::get_lattice_len()
{
	return _L;
}

//returns the magnetization of the current state
double lattice::get_magnetization()
{
	return magnetization.back();
}

//returns the temparature of the current state
double lattice::get_temparature()
{
	return _T;
}


// Displays the lattice parameters
void lattice::display_lattice()
{
	cout<<"\nLATTICE PARAMETERS\n_ _ _ _ _ _ _ _ _ _ _ _ _\n";
	cout<<"Lattice Diamension\t\t=\t"<<_dimension<<"\n";
	cout<<"Lattice Length\t\t\t=\t"<<_L<<"\n";
	cout<<"Total Number of sites\t\t=\t"<<_site_count<<"\n";
	cout<<"No. of possible spin states\t=\t"<<_q<<"\n";
	cout<<"Desired Temparature\t\t=\t"<<_T<<"\n";
	cout<<"Initialization\t\t\t=\t"<<_initalization<<"\n";
	cout<<"Sweeps stored\t\t\t=\t"<<get_sweep_count()<<"\n";
	cout<<"Energy\t\t\t\t=\t"<<get_energy()<<"\n";
	cout<<"Magnetization\t\t\t=\t"<<get_magnetization()<<"\n";
	cout<<"_ _ _ _ _ _ _ _ _ _ _ _\n";
}



// Displays the lattice parameters , lattice and neighbour  map
void lattice::display_lattice_detailed()
{
	cout<<"\n_ _ _ _ _ _ _ _ _ _ _ _ _\n\tLATTICE PARAMETERS\n";
	
	
	cout<<"Lattice Diamension\t\t=\t"<<_dimension<<"\n";
	cout<<"Lattice Length\t\t\t=\t"<<_L<<"\n";
	cout<<"Total Number of sites\t\t=\t"<<_site_count<<"\n";
	cout<<"No. of possible spin states\t=\t"<<_q<<"\n";
	cout<<"Desired Temparature\t\t=\t"<<_T<<"\n";
	cout<<"Initialization\t\t\t=\t"<<_initalization<<"\n";
	cout<<"Energy\t\t\t\t=\t"<<get_energy()<<"\n";
	cout<<"Magnetization\t\t\t=\t"<<get_magnetization()<<"\n";
	cout<<"Initialization\t\t\t=\t"<<_initalization<<"\n";
	cout<<"_beta\t\t\t\t=\t"<<_beta<<"\n";
	cout<<"_wolf_add_probab\t\t=\t"<<_wolf_add_probab<<"\n";
	cout<<"coupling J\t\t\t=\t"<<_J<<"\n";
	cout<<"sweeps\t\t\t\t=\t"<<_sweep_count<<"\n";
	cout<<"urd seed\t\t\t=\t"<<_urd_seed<<"\n";
	cout<<"initialization seed\t\t=\t"<<_initalization_seed<<"\n";
	cout<<"_ _ _ _ _ _ _ _ _ _ _ _\n";

	for(int i=0;i<_site_count;i++)
	{
		if(i%_L==0)
			cout<<"\n";
		cout<<_state[i]._spin<<",";
	}
	cout<<"\n";
	//to display neighbour details
	/*for(int i=0;i<_site_count;i++)
	{
		cout<<"\n"<<i<<"\t";
		display_site(_state[i]);
	}*/
}

void lattice::mc_steps(int count)
{
	site *mc_site;
	double de,dm,E,M;
	spin nib_spin;
	uniform_int_distribution<int> int_dist(0,_site_count-1);
	uniform_real_distribution<double> distribution(0.0,1.0);
	update_history+="# "+to_string(_sweep_count+1)+" -> Metropolis -> ";
	for (int k=0;k<count;k++)
	{	
	E=0;
	M=0;
	int x=0;
		for (int i=0;i<_site_count;i++)
		{
			x=int_dist(generator);
			mc_site=_state+x;
			nib_spin=mc_site->_neighbour[0]->_spin;
			//cout<<"\n["<<nib_spin<<",";
			nib_spin+=mc_site->_neighbour[1]->_spin;
		//	cout<<nib_spin<<",";
			nib_spin+=mc_site->_neighbour[2]->_spin;
			//cout<<nib_spin<<",";
			nib_spin+=mc_site->_neighbour[3]->_spin;
			//cout<<nib_spin<<"] ";
			de=mc_site->_spin*(nib_spin)*(2*_J);
			dm=mc_site->_spin*(-2);
			//cout<<"x = "<<x<<" nibs = "<<nib_spin<<" de = "<<de<<", dm = "<<dm<<" , isp : "<<mc_site->_spin;
			if(de<0)
			{
				mc_site->_spin*=-1;
			}
			else
			{
				if (distribution(generator)< exp(-1*_beta*de))
				{
					mc_site->_spin*=-1;
				}
				else
				{
					de=0;
					dm=0;
				}
			}
	//	cout<<" fsp : "<<mc_site->_spin<<", fde : "<<de<<" , fdm :"<<dm <<"\n";
		E+=de;
		M+=dm;
		}
		
	add_new_energy(E);
	add_new_magnetization(M);
	_sweep_count+=1;
	}
	update_history+=to_string(_sweep_count)+"\n";
}

void lattice::wolf_steps(int count)
{
	double de,dm;
	site *wlf_seed,*current_site;
	spin wlf_seed_spin;
	stack <site*> nib_stack;
	set <site*> flip_set;
	update_history+="# "+to_string(_sweep_count+1)+" -> Wolf -> ";
	uniform_real_distribution<double> distribution(0.0,1.0);
	for(int i=0;i<count;i++)
	{
		de=0.0;
		dm=0.0;
//		cout<<"\n_________ "<<i<<"\n";
		while(!nib_stack.empty())
		{
			nib_stack.pop();
		}
		flip_set.clear();

		wlf_seed=_state+rand()%_site_count;
		wlf_seed_spin=wlf_seed->_spin;
		nib_stack.push(wlf_seed);
		while(!nib_stack.empty())
		{
			current_site=nib_stack.top();
			nib_stack.pop();
			if(flip_set.count(current_site)!=0)
				{
					continue;
					cout<<"c";
				}
			flip_set.insert(current_site);
			for(int j=0;j<4;j++)
			{
				if(wlf_seed_spin!=current_site->_neighbour[j]->_spin)
					continue;
			
				if(distribution(generator)< _wolf_add_probab)
				{
//					cout<<"adding "<<current_site->_neighbour[j]<<"to stack\n";
					nib_stack.push(current_site->_neighbour[j]);
				}
			}
		}
		int border_spin=0;
		for (set<site*> :: iterator itr = flip_set.begin(); itr != flip_set.end(); itr++)
		{
			for(int j=0;j<4;j++)
				if(flip_set.count((*itr)->_neighbour[j])==0)
					border_spin+=(*itr)->_neighbour[j]->_spin;
			(*itr)->_spin *=-1;
		/*	cout<<"\n s ="<<s<<" flipping "<<(*itr)->_spin;
			(*itr)->_spin= -s;
			cout<<" TO "<<(*itr)->_spin;*/
		}
		de=wlf_seed_spin*(border_spin)*(2*_J);
		if(!flip_set.empty())
			dm=wlf_seed_spin*int(flip_set.size())*(-2);
		else
			dm=0;
		add_new_energy(de);
		add_new_magnetization(dm);
	_sweep_count++;
	}
	update_history+=to_string(_sweep_count)+"\n";
}

void lattice::mc_potts_steps(int count)
{
	site *mc_site;
	double de,dm,E,M;
	uniform_int_distribution<int> int_dist(0,_site_count-1);
	uniform_real_distribution<double> distribution(0.0,1.0);
	update_history+="# "+to_string(_sweep_count+1)+" ->PottsMetropolis -> ";
	for (int k=0;k<count;k++)
	{	
		E=0;
		M=0;
		int nib_spin=0,n_nib_spin=0;
		int nsp,osp;
		for (int i=0;i<_site_count;i++)
			{
			//	cout<<" X = "<<x<<" ";
				mc_site=_state+int_dist(generator);
				osp=mc_site->_spin;
				nsp=osp;
				while(osp==nsp)
					nsp=rand()%_q;;
			
				n_nib_spin=0;
			//	cout<<"osp,nsp : ("<<osp<<","<<nsp<<") : q : "<< _q;
				if(nsp==mc_site->_neighbour[0]->_spin)
					n_nib_spin+=1;
				if(nsp==mc_site->_neighbour[1]->_spin)
					n_nib_spin+=1;
				if(nsp==mc_site->_neighbour[2]->_spin)
					n_nib_spin+=1;
				if(nsp==mc_site->_neighbour[3]->_spin)
					n_nib_spin+=1;

				nib_spin=0;
				if(osp==mc_site->_neighbour[0]->_spin)
					nib_spin+=1;
				if(osp==mc_site->_neighbour[1]->_spin)
					nib_spin+=1;
				if(osp==mc_site->_neighbour[2]->_spin)
					nib_spin+=1;
				if(osp==mc_site->_neighbour[3]->_spin)
					nib_spin+=1;
			
				de=(n_nib_spin-nib_spin)*(-1*_J);
				dm=(nsp-osp);
			//	cout<<" de = "<<de;
				if(de<0)
				{
					mc_site->_spin=nsp;
				}
				else
				{
					if (distribution(generator)< exp(-1*_beta*de))
					{
						mc_site->_spin=nsp;
					}
					else
					{
			//			cout<<"  NOT CHANGING !!";
						de=0;
						dm=0;
					}
				}
		//	cout<<"fde : "<<de<<" FINAL STATE : "<<mc_site->_spin<<" Ecalc = "<<calc_potts_energy();
		//	cout<<"\n";
			E+=de;
			M+=dm;
	}
		
	add_new_energy(E);
	add_new_magnetization(M);
	_sweep_count+=1;
	}
	update_history+=to_string(_sweep_count)+"\n";
}

void lattice::wolf_potts_steps(int count)
{
	double de,dm;
	site *wlf_seed,*current_site;
	spin wlf_seed_spin;
	stack <site*> nib_stack;
	set <site*> flip_set;
	update_history+="# "+to_string(_sweep_count+1)+" -> PottWolf -> ";
	uniform_real_distribution<double> distribution(0.0,1.0);
	uniform_int_distribution<int> int_dist(0,_site_count-1);
	for(int i=0;i<count;i++)
	{
		de=0.0;
		dm=0.0;
//		cout<<"\n_________ "<<i<<"\n";
		if(distribution(generator)> _wolf_add_probab_potts)
		{
					add_new_energy(de);
					add_new_magnetization(dm);
					continue;
		}
		while(!nib_stack.empty())
		{
			nib_stack.pop();
		}
		flip_set.clear();

		wlf_seed=_state+int_dist(generator)%_site_count;
		wlf_seed_spin=wlf_seed->_spin;
		
		nib_stack.push(wlf_seed);
		
		while(!nib_stack.empty())
		{
			current_site=nib_stack.top();
			nib_stack.pop();
			if(flip_set.count(current_site)!=0)
				{
					continue;
				}
			flip_set.insert(current_site);
			for(int j=0;j<4;j++)
			{
				if(wlf_seed_spin!=current_site->_neighbour[j]->_spin)
					continue;
			
				if(distribution(generator)< _wolf_add_probab_potts)
				{
	//				cout<<"adding "<<current_site->_neighbour[j]->_spin<<"to stack\n";
					nib_stack.push(current_site->_neighbour[j]);
				}
			}
		}
		
		int border_spin_old=0;
		int border_spin_new=0;
		int nsp=wlf_seed_spin;
		
		while(nsp==wlf_seed_spin)
		{
			nsp=rand()%_q;
		}
		for (set<site*> :: iterator itr = flip_set.begin(); itr != flip_set.end(); itr++)
		{
			for(int j=0;j<4;j++)
				if(flip_set.count((*itr)->_neighbour[j])==0)
				{
				//	cout<<","<<(*itr)->_neighbour[j]->_spin;
					if((*itr)->_neighbour[j]->_spin==wlf_seed_spin)
						border_spin_old+=1;
					if((*itr)->_neighbour[j]->_spin==nsp)
						border_spin_new+=1;
				}
		//	cout<<","<<(*itr)->_spin<<"->";
			(*itr)->_spin=nsp;
			//cout<<(*itr)->_spin;			
		}
		de=(-1*_J)*(border_spin_new-border_spin_old);
		dm=int(flip_set.size())*(nsp-wlf_seed_spin);
	//	cout<<" osp = "<<wlf_seed_spin<<" nsp : "<<nsp<<"   de : "<<de<<"  dm "<<dm;
		add_new_energy(de);
		add_new_magnetization(dm);
	_sweep_count++;
//	cout<<"\n";
	}
	update_history+=to_string(_sweep_count)+"\n";
}




void lattice::writeToText(string fname)
{
	fstream ofile;
	ofile.open(fname.c_str(),ios::out);
	ofile<<"#T   , "<<_T<<"\n";
	ofile<<"#Q   , "<<_q<<"\n";
	ofile<<"#Initialization  ,"<<_initalization<<"\n";
	ofile<<"#Number of sites  , "<<_site_count<<"\n";
	ofile<<"#length of lattice  , "<<_L<<"\n";
	ofile<<"#coupling _beta  , "<<_J<<"\n";
	ofile<<"#coupling _wolf_add_probab  , "<<_wolf_add_probab<<"\n";
	ofile<<"#coupling _wolf_add_probab_potts  , "<<_wolf_add_probab_potts<<"\n";
	ofile<<"#coupling J  , "<<_J<<"\n";
	ofile<<"#sweeps  , "<<_sweep_count<<"\n";
	ofile<<"#urd seed  , "<<_urd_seed<<"\n";
	ofile<<"#initialization seed , "<<_initalization_seed<<"\n";
	ofile<<"#remarks , file made with the c++ code"<<"\n";
	//ofile<<"# UPDATE HISTORY \n"<<update_history;
	ofile<<"#sweep,E,M"<<"\n";
	
	vector<double>::iterator en_it,mg_it; 
	en_it = energy.begin();
	mg_it = magnetization.begin();
	int i=0;
	while(en_it!=energy.end() and mg_it!=magnetization.end())
	{
		ofile<<i<<","<<*en_it<<","<<*mg_it<<"\n";
		i++;
		mg_it++;
		en_it++;
	}
	ofile<<"#fstate\n";
	for(int i=0;i<_site_count;i++)
	{
		ofile<<_state[i]._spin;
		if((i+1)%_L==0)
			ofile<<"\n";
		else
			ofile<<",";
	}
	ofile.close();
}

lattice_write lattice::get_parameters()
{
  lattice_write params;
		params._T=_T;
		params._beta=_beta;
		params._J=_J;
		params._wolf_add_probab=_wolf_add_probab;
		params._wolf_add_probab_potts=_wolf_add_probab_potts;
		params._dimension=_dimension;
		params._site_count=_site_count;
		params._q=_q;
		params._L=_L;
		params._initalization=_initalization;
		params._sweep_count=_sweep_count;
		params._urd_seed=_urd_seed;
		params._initalization_seed=_initalization_seed;
		strcpy(params.update_history,update_history.c_str());
	//	cout<<"\n"<<params.update_history<<" -< \n";
		return params;
}










