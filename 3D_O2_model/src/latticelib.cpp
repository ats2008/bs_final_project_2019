#include "latticelib.h"

//void display_site(site s)
//{
//	cout<<"spin = "<<s._spin<<" Neibs : "<<(s._neighbour[0])->_spin;
//	for(int i=1;i<4;i++)
//	{
//	cout<<","<<(s._neighbour[i])->_spin;
//	}
//}

///         class :: lattice

//default constructor for lattice
lattice::lattice()
{
	_state =NULL;
	_urd_seed= time(NULL);
	_initalization_seed=_urd_seed+1;
	generator.seed(_urd_seed);
	update_history="";
}

//default destructor
void lattice::clear_lattice()
{
	if(_state!=NULL) delete [] _state;
}


//sets the temparature of the current state
void lattice::set_temparature(double temparature)
{
	_T=temparature;
	_beta=1.0/_T;
}

//sets the J of the current state
void lattice::set_J(double J)
{
	_J=J;
}

//sets the sweep count
void lattice::set_sweep_count(int swc)
{
	_sweep_count=swc;
}


void lattice::seed_generator(int seeder,int initialization_seed)
{
	if (initialization_seed!=-1)
		_initalization_seed=initialization_seed;
	_urd_seed=seeder;	
	generator.seed(seeder);
}

// clears the content of energy and magnetization keeping the last one intact;
void lattice::clear_data(bool full)
{
	double tem;
	if (!energy.empty())
		{	
			tem=energy.back();
			energy.clear();
			if(!full)
			energy.push_back(tem);	
		}
	if (!magnetization.empty())
		{
			tem=magnetization.back();
			magnetization.clear();
			if(!full)
			magnetization.push_back(tem);
		}
//	update_history="\n";
}  


// returns the diamension of the lattice
int lattice::get_diamension()
{
	return _dimension;
}
//returns the uniform distribution seed
int lattice::get_urd_seed()
{
	return _urd_seed;
}

//returns the initialization_seed
int lattice::get_initialization_seed()
{
	return _initalization_seed;
}

//returns the energy of the current state
double lattice::get_energy()
{
	return energy.back();
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
//	cout<<"here in get mag "<<"\n";
	return magnetization.back();
}

//returns the temparature of the current state
double lattice::get_temparature()
{
	return _T;
}

//returns the beta of the current state
double lattice::get_beta()
{
	return _beta;
}

//returns the J the current simulation
double lattice_O_n::get_J()
{
	return _J;
}

// adds energy of the new state using change de
void lattice::add_energy_change(double de)
{
	if(energy.empty())
	{
//		cout<<"ADDING NEW ENERGY TO EMPTY vector\n";
		energy.push_back(de);
	}
	else
	{
//		cout<<"adding delta e = "<< de<<"\n";
		energy.push_back(energy.back() + de);
	}
}


// Displays the basic lattice parameters
void lattice::display_lattice()
{
	cout<<"\nLATTICE PARAMETERS\n_ _ _ _ _ _ _ _ _ _ _ _ _\n";
	cout<<"Lattice Diamension\t\t=\t"<<_dimension<<"\n";
	cout<<"Lattice Length\t\t\t=\t"<<_L<<"\n";
	cout<<"Total Number of sites\t\t=\t"<<_site_count<<"\n";
	cout<<"Desired Temparature\t\t=\t"<<_T<<"\n";
	cout<<"Initialization\t\t\t=\t"<<_initalization<<"\n";
	cout<<"Energy\t\t\t\t=\t"<<get_energy()<<"\n";
	cout<<"Magnetization\t\t\t=\t"<<get_magnetization()<<"\n";
	cout<<"Initialization\t\t\t=\t"<<_initalization<<"\n";
	cout<<"_beta\t\t\t\t=\t"<<_beta<<"\n";
	cout<<"coupling J\t\t\t=\t"<<_J<<"\n";
	cout<<"sweeps\t\t\t\t=\t"<<_sweep_count<<"\n";
	cout<<"urd seed\t\t\t=\t"<<_urd_seed<<"\n";
	cout<<"initialization seed\t\t=\t"<<_initalization_seed<<"\n";
}



// Displays the lattice parameters , lattice and neighbour  map
void lattice::display_lattice_detailed()
{
	display_lattice();
}


void lattice::display_neighbours()
{
	cout<<"\n Neighbour List \n";
	vector <site *>::iterator nib_itr;
	for(int i=0;i<_site_count;i++)
	{
		cout<<"["<<_state[i]._id<<"]\t:\t{ ";
		nib_itr=_state[i]._neighbour.begin();
		cout<<(*nib_itr)->_id;
		nib_itr++;
		for(;nib_itr!=_state[i]._neighbour.end();nib_itr++)
			cout<<","<<(*nib_itr)->_id;
		cout<<" }\n";
	}
}

///*********************** BELOW FUNCTION NEEDS TO BE GENERALIZED ...(done)  ******************************************///

// Function that initialize the neighbours of the site passed  .. pointers to nib at asite._neighbour
void lattice::set_neighboures()
{
	int * shift_parts = new int [_dimension];
	shift_parts[0]=1;
//	cout<<"shift parts = "<<0<<"  :  "<<shift_parts[0]<<"\n";
	for(int i=1;i<_dimension;i++)
	{
		shift_parts[i]=shift_parts[i-1]*_L;
//		cout<<"shift parts = "<<i<<"  :  "<<shift_parts[i]<<"\n";
	}
//	cout<<"setting neighbours for dim = "<<_dimension<<" , L = "<<_L<<", site_count = "<<_site_count<<"\n";
	for (int i=0;i<_site_count;i++)
	{
		for(int j=0;j<_dimension;j++)
		{
			_state[i]._neighbour.push_back(_state + (i-shift_parts[j]+_site_count)%_site_count);
			_state[i]._neighbour.push_back(_state + (i+shift_parts[j])%_site_count);
//			cout<<" nebs set for "<<_state[i]._id<<" : "<<(_state+(i-shift_parts[j]+_site_count)%_site_count)->_id<<"\n";
//			cout<<" nebs set for "<<_state[i]._id<<" : "<<(_state+(i+shift_parts[j]+_site_count)%_site_count)->_id<<"\n";
		
		}
	}
	delete [] shift_parts;
}

///*********************** ABOVE FUNCTION NEEDS TO BE GENERALIZED  ******************************************///

//////////////////////////////////  class : O(n) Model ///////////////////////////////////


lattice_O_n::lattice_O_n()
{
	_H=NULL;
}
void lattice_O_n::initialize(char mode)
{
	site *itr=_state;
	int id=0;
	if(mode=='Z')
	{
		for(int i=0;i<_site_count;i++)
			{
				id++;
				(itr+i)->_id=id;
				(itr+i)->_components= new double [_n];
				(itr+i)->_components[0]=1.0;
				for(int j=1;j<_n;j++)
					(itr+i)->_components[j]=0.0;
			}
	}
	if(mode=='H')
	{
		uniform_real_distribution<double> distribution(0.0,1.0);
		generator.seed(_initalization_seed);
		distribution(generator);
		double temp;
		for(int i=0;i<_site_count;i++)
			{
			
				id++;
				(itr+i)->_id=id;
				temp=0;
				(itr+i)->_components= new double [_n];
				for(int j=0;j<_n;j++)
				{
					(itr+i)->_components[j]=distribution(generator)-0.5;
					temp+=(itr+i)->_components[j]*(itr+i)->_components[j];
				}
				temp=sqrt(temp);
				for(int j=0;j<_n;j++)
					(itr+i)->_components[j]/=temp;
			}
		generator.seed(_urd_seed);
	}
}


void lattice_O_n::setup_model(int dim,int n,int lattice_len,int N,int nib_count,
						double T_equilib,double J,double* H,char initialization,char model_type)
{
	_dimension=dim;
	_L=lattice_len;
	_site_count =N;
	_nib_count=nib_count;
	_n=n;
	_J=J;
	_T=T_equilib;
	_beta=1.0/_T;
	_initalization=initialization;
	_model_type = model_type;

	if(N>0)	_state = new site[N];
	if(n>0) _H =new double [n];
	for(int i=0;i<n;i++)
		_H[i]=H[i];
//	cout<<"H set"<<" \n ";
	initialize(_initalization);
//	cout<<"initialization set"<<" \n ";
	/*********************************************************/
	set_neighboures();
	/*********************************************************/
//	cout<<"Neibs set"<<" \n ";
	add_energy_change(calc_energy());
//	cout<<"Calc energy set as : "<<get_energy()<<" \n ";
//	cout<<"_n = "<<_n<<"\n";
	double * mag ;
	mag= new double [_n];
	calc_magnetization(mag);
	add_magnetization_change(mag);
	delete [] mag;
//	cout<<"Magnetization set as :\n ";
//	cout<<get_magnetization()<<" \n ";
	_sweep_count=0;
	_urd_seed= time(NULL);
	generator.seed(_urd_seed);
	update_history="\n";

}

//sets the ext-field H
void lattice_O_n::set_ext_field(double * H)
{
	for(int j=0;j<_n;j++)
		_H[j]=H[j];
}

void lattice_O_n::get_ext_field(double * H)
{
	for(int j=0;j<_n;j++)
		H[j]=_H[j];
}

int lattice_O_n::get_n()
{
	return  _n;
}

//returns the ext-field as string
string lattice_O_n::get_H_as_string()
{
	string hs;
	hs="[";
	hs+=to_string(round(_H[0]*100000)/100000).substr(0,7);
	for(int j=1;j<_n;j++)
	{
		hs+="|";
		hs+=to_string(round(_H[j]*1000)/1000).substr(0,7);
	}
	hs+="]";
	return hs;
}

// returns the _n
int lattice_O_n::get_spin_dim()
{
	return _n;
}
		

double lattice_O_n::calc_energy()
{
	double en=0;
	double* temp =new double [_n];
	for(int i=0;i<_site_count;i++)
		{
			for(int j=0;j<_n;j++)
						temp[j]=0;
			for(vector <site *>::iterator itr =_state[i]._neighbour.begin();
										  itr !=_state[i]._neighbour.end() ;itr++ )
				{
					for(int j=0;j<_n;j++)
						temp[j]+=(*itr)->_components[j];
				}
			for(int j=0;j<_n;j++)
				{
					en += -0.5*_J * _state[i]._components[j]*temp[j];
					en -= _H[j] * _state[i]._components[j];
				
				}
		}
	delete [] temp;
	return en;
}


void lattice_O_n::calc_magnetization(double* mag)
{
	for(int j=0;j<_n;j++)
				mag[j]=0;
//	cout<<"calc_magnetization with _site_count = "<<_site_count	<<"\n";
	
	for(int i=0;i<_site_count;i++)
		{
	    	for(int j=0;j<_n;j++)
				mag[j]+= _state[i]._components[j];
		}
//	cout<<"OUT\n";
}

void lattice_O_n::add_magnetization_change(double * dmag)
{
	
	double temp=0;
	double *newmag;
	
	if(magnetization_O_n.empty())
	{
		newmag=new double [_n];
		for(int j=0;j<_n;j++)
		{
			newmag[j]=dmag[j];
			temp+=newmag[j]*newmag[j];
//		cout<<" "<<dmag[j]<<" ["<<newmag[j]<<"]";
		}
//		cout<<"->dm  "<<"[newmag]"<<"\n";
		magnetization_O_n.push_back(newmag);
		magnetization.push_back(sqrt(temp));
	}
	else
	{
		double *oldmag;
		newmag=new double [_n];
		oldmag= magnetization_O_n.back();
		for(int j=0;j<_n;j++)
		{

			newmag[j]=oldmag[j]+dmag[j];
			temp+=newmag[j]*newmag[j];
//			cout<<" "<<dmag[j]<<" ["<<newmag[j]<<"]";
		}
//		cout<<"-> dm \n";
		magnetization_O_n.push_back(newmag);
		magnetization.push_back(sqrt(temp));
	}
}



//void lattice::writeLatticeToText(string fname)
//{
//	fstream ofile;
//	ofile.open(fname.c_str(),ios::out|ios::app);
//	ofile<<"#fstate\n";
//	for(int i=0;i<_site_count;i++)
//	{
//		ofile<<_state[i]._spin;
//		if((i+1)%_L==0)
//			ofile<<"\n";
//		else
//			ofile<<",";
//	}
//	ofile.close();
//}


void lattice_O_n::writeToText(string fname,bool lattice_skip)
{
	fstream ofile;
	ofile.open(fname.c_str(),ios::out);
	ofile<<"#dim,"<<_dimension<<"\n";
	ofile<<"#n,"<<_n<<"\n";
	ofile<<"#length of lattice, "<<_L<<"\n";
	ofile<<"#Number of sites, "<<_site_count<<"\n";
	ofile<<"#T, "<<_T<<"\n";
	ofile<<"#beta, "<<_beta<<"\n";
	ofile<<"#J, "<<_J<<"\n";
//	cout<<"here in txtwrite\n";
	for(int i=0;i<_n;i++)
		ofile<<"#H"<<i<<","<<_H[i]<<"\n";
	ofile<<"#Initialization,"<<_initalization<<"\n";
	ofile<<"#sweeps, "<<_sweep_count<<"\n";
	ofile<<"#urd_seed, "<<_urd_seed<<"\n";
	ofile<<"#initialization_seed, "<<_initalization_seed<<"\n";
	ofile<<"#remarks , file made with the c++ code : O(n) model"<<"\n";
	ofile<<"#UPDATE HISTORY , "<<update_history;
	ofile<<"#sweep,E,M";
	for(int i=0;i<_n;i++)
		ofile<<",M"<<i;
	ofile<<"\n";

	vector<double>::iterator en_it,mg_it;
	vector<double *>::iterator mgn_it; 
	en_it = energy.begin();
	mg_it = magnetization.begin();
	mgn_it = magnetization_O_n.begin();
	int i=0;
//	cout<<"starting the data write\n";
	while(en_it!=energy.end() and mg_it!=magnetization.end()  and mgn_it != magnetization_O_n.end())
	{
		ofile<<i<<","<<*en_it<<","<<*mg_it;
		for(int j=0;j<_n;j++)
			ofile<<","<<(*mgn_it)[j];
		ofile<<"\n";
		i++;
		en_it++;
		mg_it++;
		mgn_it++;
	}
//	cout<<"starting the lattice write\n";
	ofile<<"#fstate\n";
	if(!lattice_skip)
	for(int i=0;i<_site_count;i++)
	{
		ofile<<_state[i]._components[0];
		for(int j=1;j<_n;j++)
			ofile<<","<<_state[i]._components[j];
		ofile<<"\n";
	}
	ofile.close();
}

lattice_write_O_n lattice_O_n::get_parameters()
{
		lattice_write_O_n params;
		params._dimension=_dimension;
		params._L=_L;
		params._site_count=_site_count;
		params._nib_count=_nib_count;
		params._n=_n;
		params._T=_T;
		params._beta=1/_T;
		params._J=_J;
		params._model_type=_model_type;
		params._initalization=_initalization;
		params._sweep_count=_sweep_count;
		params._urd_seed=_urd_seed;
		params._initalization_seed=_initalization_seed;
		strcpy(params.update_history,update_history.c_str());
		for(int i=0;i<_n;i++)
		{
			params._H[i]=_H[i];
			if(i>=MAX_H_LEN) break;
		}
		return params;
}

void lattice_O_n::writeToBinary(string fname,bool lattice_skip)
{
	fstream ofile;
	ofile.open(fname.c_str(),ios::out | ios::binary);
	lattice_write_O_n parameters =get_parameters();
	if(lattice_skip)
		parameters._site_count=0;
	ofile.write((char *)(&parameters),sizeof(lattice_write_O_n));
	
	int site_count=parameters._site_count;
	for(int i=0;i<site_count;i++)
		for(int j=0;j<_n;j++)
		ofile.write((char *)(&_state[i]._components[j]),sizeof(double));
		
	vector<double>::iterator en_it,mg_it;
	vector<double * >::iterator mgn_it; 
	en_it = energy.begin();
	mg_it = magnetization.begin();
	mgn_it = magnetization_O_n.begin();
	double d=0;
	while(en_it!=energy.end() and mg_it!=magnetization.end() and mgn_it!=magnetization_O_n.end())
	{
		d=*en_it;
		ofile.write((char *)(&d),sizeof(double));
		d=*mg_it;
		ofile.write((char *)(&d),sizeof(d));
		for(int j=0;j<_n;j++)
			{
				d=(*mgn_it)[j];
				ofile.write((char *)(&d),sizeof(d));
			}
		en_it++;
		mg_it++;
		mgn_it++;
	}
//	cout<<"binary write done !\n";
	ofile.close();
}

void lattice_O_n::readFromBinary(string fname)
{
	fstream infile;
	infile.open(fname.c_str(),ios::in | ios::binary);
	if(!infile.is_open())
	{
		cout<<"\n ERROR HAPPEND !! FILE DOES NOT EXIST !! \n";
		cout<<"fname : "<<fname;
		return;
	}
	lattice_write_O_n parameters;
//	cout<<"reading the parameters\n";
	infile.read((char *)(&parameters),sizeof(lattice_write_O_n));
//	cout<<"parmeters read \n";
//	clear_lattice();
//	cout<<"cleared the lattice\n";
	setup_model(parameters._dimension,parameters._n,parameters._L,parameters._site_count,parameters._nib_count,parameters._T,parameters._J,
							parameters._H,'Z',parameters._model_type);
//	cout<<"setuped the lattice\n";
	update_history=parameters.update_history;
	seed_generator(parameters._urd_seed,parameters._initalization_seed);
	
	double d;
	for(int i=0;i<_site_count;i++)
		for(int j=0;j<_n;j++)
		{	
			infile.read((char *)(&d),sizeof(d));
			_state[i]._components[j]=d;
		}
	
	energy.clear();
	magnetization.clear();
	magnetization_O_n.clear();
	double * temp;

	infile.read((char *)(&d),sizeof(d));
	while(!infile.eof())
	{	
		energy.push_back(d);
		
		infile.read((char *)(&d),sizeof(d));
		magnetization.push_back(d);
		
		temp =new double [_n];
		for(int j=0;j<_n;j++)
		{	
			infile.read((char *)(&d),sizeof(d));
			temp[j]=d;
		}
		magnetization_O_n.push_back(temp);

		infile.read((char *)(&d),sizeof(d));
	}
	infile.close();
}

void lattice_O_n::clear_data(bool full)
{
	vector<double *>::iterator itr=magnetization_O_n.begin();
	double *temp,*ptr;
	if(!magnetization_O_n.empty())
	{
		ptr = magnetization_O_n.back();
		if(!full)
		{
			temp =new double [_n];
			for(int j=0;j<_n;j++)
			 {
				temp[j]=ptr[j];
			 }
	//		cout<<"going to delete the allocated n-mags\n";
		}
		while(itr != magnetization_O_n.end())
		{
			delete [] *itr;
			itr++;
		}
		magnetization_O_n.clear();
		if(!full)
			magnetization_O_n.push_back(temp);
	}
	
//	cout<<"exiting clear data of child going to parant\n";
	lattice::clear_data(full);
//	cout<<"exiting clear data \n";
}

void lattice_O_n::clear_lattice()
{
	clear_data(true);
	
	if (_H !=NULL)	delete [] _H;
	
	if(_state !=NULL)
	{
		for(int i=0;i<_site_count;i++)
		{
			delete [] _state[i]._components;
			_state[i]._components=NULL;
			_state[i]._neighbour.clear();
		}
	delete [] _state;
	}
	
	_state=NULL;
	_H=NULL;
	lattice::clear_lattice();
}

void lattice_O_n::display_lattice()
{
	lattice::display_lattice();
	cout<<"Magnetization [ n ]  \t\t=\t";
	if(!magnetization_O_n.empty())	
		{
			double  *temp =magnetization_O_n.back();
			cout<<"[ "<<temp[0];
			for(int i=1;i<_n;i++)
				cout<<","<<temp[i];
			cout<<" ]\n";
		}
	else 
		cout<<" EMPTY ! "<<"\n";
		
	cout<<"number spin componets\t\t=\t"<<_n<<"\n";
	cout<<"External Field       \t\t=\t[ "<<_H[0];
	for(int i=1;i<_n;i++)
		cout<<","<<_H[i];
	cout<<" ]\n";
}

void lattice_O_n::display_lattice_detailed()
{
	lattice::display_lattice_detailed();
	cout<<"Magnetization [ n ]  \t\t=\t";
	if(!magnetization_O_n.empty())	
		{
			double  *temp =magnetization_O_n.back();
			cout<<"[ "<<temp[0];
			for(int i=1;i<_n;i++)
				cout<<","<<temp[i];
			cout<<" ]\n";
		}
	else 
		cout<<" EMPTY ! "<<"\n";
		
	cout<<"number spin componets\t\t=\t"<<_n<<"\n";
	cout<<"External Field       \t\t=\t[ "<<_H[0];
	for(int i=1;i<_n;i++)
		cout<<","<<_H[i];
	cout<<" ]\n";
	cout<<"\t\t\t STATE\n";
	
	for(int i=0;i<_site_count;i++)
	{	
		cout<<"["<<_state[i]._id<<"]\t:\t{ ";
		cout<<_state[i]._components[0];
		for(int j=1;j<_n;j++)
		{
			cout<<","<<_state[i]._components[j];
		}
		cout<<"} \n";
	}
	
}

lattice_O_n lattice_O_n::get_bare_copy()
{
	lattice_O_n clat;
	clat.setup_model(_dimension,_n,_L,_site_count,_nib_count,_T,_J,_H,_initalization,_model_type);
	clat.seed_generator(_urd_seed,_initalization_seed);
	clat.clear_data();
	clat.set_neighboures();
	
	clat.energy.clear();
	clat.magnetization.clear();
	clat.magnetization_O_n.clear();
	clat.energy.push_back(energy.back());
	clat.magnetization.push_back(magnetization.back());
	double *ptr,*temp=new double [_n];
	ptr=magnetization_O_n.back();
	for(int j=0;j<_n;j++)
		temp[j]=ptr[j];
	clat.magnetization_O_n.push_back(temp);
	site *csite=clat.get_state();
	site *osite=_state;
	for(int i=0;i<_site_count;i++)
	{
		for(int j=0;j<_n;j++)
			csite->_components[j]=osite->_components[j];
		csite++;
		osite++;
	}
	clat.update_history=update_history;
	return clat;
}

void lattice_O_n::mc_steps(int count)
{

//	display_lattice_detailed();
	site *mc_site;
	double de,*M,*dm,E,temp_d;
	M=new double [_n];
	dm=new double [_n];
	double *new_spin_components,*nib_spin_components;
	new_spin_components=new double [_n];
	nib_spin_components=new double [_n];
	
	
	uniform_int_distribution<int> int_dist(0,_site_count-1);
	uniform_real_distribution<double> distribution(0.0,1.0);
	
	vector<site *>::iterator nib_itr; 
	
	for (int k=0;k<count;k++)
	{
		E=0;
		for(int j=0;j<_n;j++)
			M[j]=0.0;
		
		for (int i=0;i<_site_count;i++)
		{
			int yu=int_dist(generator);
			mc_site =_state+yu;
			for(int j=0;j<_n;j++)
				{
					dm[j]=0.0;
					nib_spin_components[j]=0;
				}
			nib_itr=(mc_site->_neighbour).begin();
//			cout<<"\nNEIB COUNT : "<<(mc_site->_neighbour).size()<<" _n = "<<_n<<"  _J = "<<_J<<" beta "<<_beta<<" yu = "<<yu<<"\n";
			for(; nib_itr!=(mc_site->_neighbour).end();nib_itr++)
			{
				for(int j=0;j<_n;j++)
					{
						nib_spin_components[j]+=(*nib_itr)->_components[j];
					}
			}
			
			temp_d=0;
			for(int j=0;j<_n;j++)
			{
				new_spin_components[j]=distribution(generator)-0.5;
				temp_d+=new_spin_components[j]*new_spin_components[j];
			}
			
			temp_d=sqrt(temp_d);
			de=0.0;
			for(int j=0;j<_n;j++)
			{
				new_spin_components[j]/=temp_d;
				dm[j]=new_spin_components[j] - mc_site->_components[j];
				de+=(-1.0*_J*dm[j]*nib_spin_components[j] - 1.0*_H[j]*dm[j]);
			}	
			if(de<0)
			{
				for(int j=0;j<_n;j++)
				{
					mc_site->_components[j]=new_spin_components[j];
				}
			}
			else
			{
				double rt=distribution(generator);
				if (rt< exp(-1*_beta*de))
				{
					for(int j=0;j<_n;j++)
					{
						mc_site->_components[j]=new_spin_components[j];
					}
				}
				else
				{
					de=0.0;
					for(int j=0;j<_n;j++)
						dm[j]=0.0;
				}
			}
			
			E+=de;
			for(int j=0;j<_n;j++)
						M[j]+=dm[j];
		}
		add_energy_change(E);
		add_magnetization_change(M);
		_sweep_count+=1;
	}
	
	delete []	M;
	delete []  dm;
	delete [] new_spin_components;
	delete [] nib_spin_components;
}

void lattice_O_n::wolf_steps(int count)
{
	double DE,de,*dm,*DM,temp_d;
	double *wlf_seed_spin_components,*temp_x;
	double component_sign,check,wlf_prob;
	uniform_real_distribution<double> distribution(0.0,1.0);
//	display_lattice();
	site *wlf_seed,*current_site;
	stack <site*> nib_stack;
	vector <site*>::iterator nib_itr;
	set <site*> flip_set;

	wlf_seed_spin_components= new double [_n];
	dm  = new double[_n];
	DM  = new double[_n];
	
//	cout<<"here at start of the wolf "<<"\n";

	for(int i=0;i<count;i++)
	{	
//		cout<<"@ "<<i<<"  wlf step and  E , M = "<<get_energy()<<"  "<<get_magnetization()<<"\n";
		
		while(!nib_stack.empty())
		{
			nib_stack.pop();
		}
		flip_set.clear();

		wlf_seed=_state+rand()%_site_count;
		nib_stack.push(wlf_seed);

		temp_d=0;
		for(int j=0;j<_n;j++)
		{
			wlf_seed_spin_components[j]=distribution(generator)-0.5;
			temp_d+=wlf_seed_spin_components[j]*wlf_seed_spin_components[j];
		}
		temp_d=sqrt(temp_d);
//		double p=0;
		for(int j=0;j<_n;j++)
		{
			wlf_seed_spin_components[j]/=temp_d;
//			p+=wlf_seed_spin_components[j]*wlf_seed_spin_components[j];
		}
//		cout<<"\t\t\t\t  p "<<p<<"\n";
		// making the cluster
		while(!nib_stack.empty())
		{
//			cout<<"at while(!nib_stack.empty())\n";
			current_site=nib_stack.top();
			nib_stack.pop();
			if(flip_set.count(current_site)!=0)
				{ 	
					continue;
				}
			flip_set.insert(current_site);
			
			component_sign=0;
			for(int j=0;j<_n;j++)
			{
				component_sign+=wlf_seed_spin_components[j]*current_site->_components[j];
			}
			
			
			nib_itr=current_site->_neighbour.begin();
			for(; nib_itr!=(current_site->_neighbour).end();nib_itr++)
			{
				check=0;
				for(int j=0;j<_n;j++)
					check+=(*nib_itr)->_components[j]*wlf_seed_spin_components[j];

				check*=component_sign;
				
				if(check<0)
					continue;
					
				wlf_prob = 1- exp(-2*_beta* _J * check);
				if(distribution(generator) < wlf_prob)
				{
					nib_stack.push(*nib_itr);
				}

			}
		}
		
		
		// Finding the flipped config variables
//		cout<<" Finding the flipped config variables"<<"\n";
		de=0.0;
		DE=0.0;
		for(int j=0;j<_n;j++)
			DM[j]=0;
		vector <double *> flipped_spins;
		for ( set<site*> :: iterator itr = flip_set.begin(); itr != flip_set.end(); itr++ )
		{	
//			cout<<" Finding the flipped config variables"<<"\n";
			temp_x = new double [_n];

			temp_d=0;
			for(int j=0;j<_n;j++)
				temp_d+=2*(*itr)->_components[j]*wlf_seed_spin_components[j] ;
			
			for(int j=0;j<_n;j++)
				temp_x[j]=(*itr)->_components[j]-temp_d*wlf_seed_spin_components[j];
			
			flipped_spins.push_back(temp_x);
//			cout<<" Pushed to flipped_spins"<<"\n";

			nib_itr=(*itr)->_neighbour.begin();
			for(; nib_itr!=((*itr)->_neighbour).end();nib_itr++)
			{
				if(flip_set.count(*nib_itr)==0)
				{
					for(int j=0;j<_n;j++)
						DE += -1.0*_J * ((temp_x[j]-(*itr)->_components[j]) *((*nib_itr)->_components[j]));
				}
			}
			for(int j=0;j<_n;j++)
			{
				dm[j]=(temp_x[j]-(*itr)->_components[j]) ;
				de+=-1.0*_H[j]*dm[j] ;
				DM[j]+=dm[j];
			}
		}
//		Deciding wether to flip the cluster
//		cout<<"Deciding wether to flip the cluster"<<"\n";

		vector <double *>::iterator d_itr;
		set<site*> :: iterator itr;
		if(de<=0)
			{
//				cout<<" FLIP LINE : ("<<wlf_seed_spin_components[0]<<","<<wlf_seed_spin_components[1]<<")  [ "<<wlf_seed_spin<<" ]\n";
				itr = flip_set.begin();
				d_itr =flipped_spins.begin();
				while (itr != flip_set.end() and d_itr!= flipped_spins.end())
				{
					for(int j=0;j<_n;j++)
						(*itr)->_components[j] = (*d_itr)[j];
					itr++;
					d_itr++;
				}
			}
		else
			{
				if (distribution(generator)< exp(-1*_beta*de))
				{
//					cout<<" FLIP LINE : ("<<wlf_seed_spin_components[0]<<","<<wlf_seed_spin_components[1]<<")  [ "<<wlf_seed_spin<<" ]\n";
					itr = flip_set.begin();
					d_itr =flipped_spins.begin();
					while (itr != flip_set.end() and d_itr!= flipped_spins.end())
					{
						for(int j=0;j<_n;j++)
							(*itr)->_components[j] = (*d_itr)[j];
						itr++;
						d_itr++;
					}
				}
				else
				{
					DE=0.0;
					de=0.0;
					for(int j=0;j<_n;j++)
						DM[j]=0.0;
				}
			}
//		cout<<"going to add energy and Magnetization"<<"\n";
		add_energy_change(de+DE);
		add_magnetization_change(DM);
//		cout<<"Added energy and Magnetization"<<"\n";
		d_itr =flipped_spins.begin();
		while(d_itr!=flipped_spins.end())
		{
				delete [] (*d_itr);
				d_itr++;
 		}
 		flipped_spins.clear();
	_sweep_count++;
	}
	
	delete [] wlf_seed_spin_components;
	delete [] DM;
	delete [] dm;
}
lattice_O_n lattice_O_n::resample_thermalized(int tau,double tau_err,int thermalization_skip,float tau_err_tolarance)
{
	int indip_sample_period=(1.3*(tau+tau_err)+1);
	cout<<"Resampling at interval of "<<indip_sample_period<<"  \n";
	if( float(tau_err/tau)> tau_err_tolarance) 
		{
			cout<<"\n Tau Have err > Tol !! "<<tau_err*100/tau<<" %"<<"\n";
		}
	lattice_O_n blattice=get_bare_copy();
	
	blattice.energy.clear();
	blattice.magnetization.clear();
	blattice.magnetization_O_n.clear();
	
	vector<double>::iterator en_it;
	vector<double>::iterator mg_it;
	vector<double *>::iterator mgn_it;
	double *new_d_vec=NULL;
	en_it = energy.begin();
	mg_it = magnetization.begin();
	mgn_it = magnetization_O_n.begin();
	int m=0;
	int i=0;
	bool Not_thermalized=true;
	while(en_it!=energy.end() and mg_it!=magnetization.end()
					and mgn_it!=magnetization_O_n.end())
	{	
		m++;
		if(Not_thermalized)
		{
			if(m>=thermalization_skip)
			{
				Not_thermalized=false;
				m=0;
			}
			en_it++;
			mg_it++;
			mgn_it++;
			continue;
		}
		if(m==indip_sample_period)
		{
			i++;
			blattice.energy.push_back(*en_it);
			blattice.magnetization.push_back(*mg_it);
			new_d_vec=new double [_n];
			for(int j=0;j<_n;j++)
				new_d_vec[j]=(*mgn_it)[j];
			blattice.magnetization_O_n.push_back(new_d_vec);
			m=0;
		}
		en_it++;
		mg_it++;
		mgn_it++;
	}
	blattice.set_sweep_count(i-1);
	return blattice;
}

lattice_O_n lattice_O_n::resample_thermalized_averages(int tau,double tau_err,int thermalization_skip,float tau_err_tolarance)
{
	int indip_sample_period=(1.3*(tau+tau_err)+1);
	cout<<"Resampling at integration window of "<<indip_sample_period<<"  \n";
	if( float(tau_err/tau)> tau_err_tolarance) 
		{
			cout<<"\n Tau Have err > Tol !! "<<tau_err*100/tau<<" %"<<"\n";
		}
	lattice_O_n blattice=get_bare_copy();
	
	blattice.energy.clear();
	blattice.magnetization.clear();
	blattice.magnetization_O_n.clear();
	
	vector<double>::iterator en_it;
	vector<double>::iterator mg_it;
	vector<double *>::iterator mgn_it;
	double *new_d_vec=NULL;
	en_it = energy.begin();
	mg_it = magnetization.begin();
	mgn_it = magnetization_O_n.begin();
	int m=0;
	int i=0;
	bool Not_thermalized=true;
	
	double temp_en,temp_mag,*temp_magn;
	temp_magn=new double [_n];
	while(en_it!=energy.end() and mg_it!=magnetization.end()
					and mgn_it!=magnetization_O_n.end())
	{	
		m++;
		if(Not_thermalized)
		{
			if(m>=thermalization_skip)
			{
				Not_thermalized=false;
				m=0;
			}
			en_it++;
			mg_it++;
			mgn_it++;
			continue;
		}
		if(m==indip_sample_period)
		{
			i++;
			blattice.energy.push_back(temp_en/indip_sample_period);
			blattice.magnetization.push_back(temp_mag/indip_sample_period);
			new_d_vec=new double [_n];
			for(int j=0;j<_n;j++)
			{
				new_d_vec[j]=temp_magn[j]/indip_sample_period;
				temp_magn[j]=0;
			}
			blattice.magnetization_O_n.push_back(new_d_vec);
			temp_en=0;
			temp_mag=0;
			m=0;
		}
		
		for(int j=0;j<_n;j++)
		{
			temp_magn[j]+=(*mgn_it)[j];
		}
		temp_en+=(*en_it);
		temp_mag+=(*mg_it);
		en_it++;
		mg_it++;
		mgn_it++;
	}
	if(m!=0)
	{
		blattice.energy.push_back(temp_en/m);
		blattice.magnetization.push_back(temp_mag/m);
		new_d_vec=new double [_n];
		for(int j=0;j<_n;j++)
		{
			new_d_vec[j]=temp_magn[j]/m;
			temp_magn[j]=0;
		}
		blattice.magnetization_O_n.push_back(new_d_vec);
	}
	blattice.set_sweep_count(i-1);
	delete [] temp_magn;
	return blattice;
}




lattice_O_n lattice_O_n::resample_data(int indip_sample_period,int off_set)
{
	cout<<"Resampling at intervals of  "<<indip_sample_period<<"  \n";
	lattice_O_n blattice=get_bare_copy();
	
	blattice.energy.clear();
	blattice.magnetization.clear();
	blattice.magnetization_O_n.clear();
	
	vector<double>::iterator en_it;
	vector<double>::iterator mg_it;
	vector<double *>::iterator mgn_it;
	double *new_d_vec=NULL;
	en_it = energy.begin();
	mg_it = magnetization.begin();
	mgn_it = magnetization_O_n.begin();
	int m=0;
	int i=0;
	
	double temp_en,temp_mag,*temp_magn;
	temp_magn=new double [_n];
	while(en_it!=energy.end() and mg_it!=magnetization.end()
					and mgn_it!=magnetization_O_n.end())
	{	
		m++;
		if(m==indip_sample_period)
		{
			i++;
			blattice.energy.push_back(temp_en/indip_sample_period);
			blattice.magnetization.push_back(temp_mag/indip_sample_period);
			new_d_vec=new double [_n];
			for(int j=0;j<_n;j++)
			{
				new_d_vec[j]=temp_magn[j]/indip_sample_period;
				temp_magn[j]=0;
			}
			blattice.magnetization_O_n.push_back(new_d_vec);
			temp_en=0;
			temp_mag=0;
			m=0;
		}
		
		for(int j=0;j<_n;j++)
		{
			temp_magn[j]+=(*mgn_it)[j];
		}
		temp_en+=(*en_it);
		temp_mag+=(*mg_it);
		en_it++;
		mg_it++;
		mgn_it++;
	}
	if(m!=0)
	{
		blattice.energy.push_back(temp_en/m);
		blattice.magnetization.push_back(temp_mag/m);
		new_d_vec=new double [_n];
		for(int j=0;j<_n;j++)
		{
			new_d_vec[j]=temp_magn[j]/m;
			temp_magn[j]=0;
		}
		blattice.magnetization_O_n.push_back(new_d_vec);
	}
	blattice.set_sweep_count(i-1);
	delete [] temp_magn;
	return blattice;
}






void lattice_O_n::bootstarap_mgn(double *rslt,int n,int bstrap_count,int bstrap_idd_count)
{
	
	uniform_int_distribution<int> int_dist(0,energy.size()-1);
	double * bstrap_array= new double [bstrap_idd_count];
	double * Avgbstrap_rslts= new double [bstrap_count];
	double * Varbstrap_rslts= new double [bstrap_count];
	int site_count=int(pow(_L,_dimension));
	
	for (int i=0;i<bstrap_count;i++)
		{
			for(int j=0;j<bstrap_idd_count;j++)
			{
					bstrap_array[j]=magnetization_O_n[int_dist(generator)][n];
			}
			Avgbstrap_rslts[i]=average(bstrap_array,bstrap_idd_count);
			Varbstrap_rslts[i]=variance(bstrap_array,bstrap_idd_count);
		}
	rslt[0]=average(Avgbstrap_rslts,bstrap_count)/site_count;
	rslt[1]=sqrt(variance(Avgbstrap_rslts,bstrap_count))/site_count;
	rslt[2]=average(Varbstrap_rslts,bstrap_count)/site_count;
	rslt[3]=sqrt(variance(Varbstrap_rslts,bstrap_count))/site_count;
	delete [] bstrap_array;
	delete [] Avgbstrap_rslts;
	delete [] Varbstrap_rslts;
}

void lattice::bootstarap_basic(double *rslt,int bstrap_count,int bstrap_idd_count)
{
	
	uniform_int_distribution<int> int_dist(0,energy.size()-1);
	
	double * Ebstrap_array= new double [bstrap_idd_count];
	double * Ebstrap_rslts= new double [bstrap_count];
	
	double * Mbstrap_array= new double [bstrap_idd_count];
	double * Mbstrap_rslts= new double [bstrap_count];

	double * SPbstrap_rslts= new double [bstrap_count];
	double * SUbstrap_rslts= new double [bstrap_count];
	
	
	for (int i=0;i<bstrap_count;i++)
		{
			for(int j=0;j<bstrap_idd_count;j++)
			{
					Ebstrap_array[j]=energy[int_dist(generator)];
					Mbstrap_array[j]=magnetization[int_dist(generator)];
			}
			
			Ebstrap_rslts[i]=average(Ebstrap_array,bstrap_idd_count);
			SPbstrap_rslts[i]=variance(Ebstrap_array,bstrap_idd_count);
			
			Mbstrap_rslts[i]=average(Mbstrap_array,bstrap_idd_count);
			SUbstrap_rslts[i]=variance(Mbstrap_array,bstrap_idd_count);
		}
	int site_count=int(pow(_L,_dimension));
	rslt[0]=average(Ebstrap_rslts,bstrap_count)/site_count;
	rslt[1]=sqrt(variance(Ebstrap_rslts,bstrap_count)/site_count);
	rslt[2]=average(Mbstrap_rslts,bstrap_count)/site_count;
	rslt[3]=sqrt(variance(Mbstrap_rslts,bstrap_count))/site_count;
	
	rslt[4]=average(SPbstrap_rslts,bstrap_count)*_beta*_beta/site_count;
	rslt[5]=sqrt(variance(SPbstrap_rslts,bstrap_count))*_beta*_beta/site_count;
	rslt[6]=average(SUbstrap_rslts,bstrap_count)*_beta/site_count;
	rslt[7]=sqrt(variance(SUbstrap_rslts,bstrap_count))*_beta/site_count;
	
	delete [] Ebstrap_array;
	delete [] Ebstrap_rslts;
	delete [] Mbstrap_array;
	delete [] Mbstrap_rslts;
	delete [] SPbstrap_rslts;
	delete [] SUbstrap_rslts;
}

void lattice_O_n::bootstarap_bindersC(double *rslt,int n,int bstrap_count,int bstrap_idd_count)
{
	uniform_int_distribution<int> int_dist(0,energy.size()-1);
	double * bstrap_array2= new double [bstrap_idd_count];
	double * bstrap_array4= new double [bstrap_idd_count];
	double * Avgbstrap_rslts= new double [bstrap_count];
	double s4,s2;
	
	for (int i=0;i<bstrap_count;i++)
		{
			for(int j=0;j<bstrap_idd_count;j++)
			{
					int r=int_dist(generator);
					bstrap_array2[j]=pow(magnetization_O_n[r][n],2);
					bstrap_array4[j]=pow(magnetization_O_n[r][n],4);
			}
			s2=average(bstrap_array2,bstrap_idd_count);
			s4=average(bstrap_array4,bstrap_idd_count);
			Avgbstrap_rslts[i]=1-s4/(3*s2*s2);
		}
	rslt[0]=average(Avgbstrap_rslts,bstrap_count);
	rslt[1]=sqrt(variance(Avgbstrap_rslts,bstrap_count));
	delete [] bstrap_array2;
	delete [] bstrap_array4;
	delete [] Avgbstrap_rslts;
	
}


lattice_O_n lattice_O_n::get_bootstrap_copy(int seed_val)
{
	lattice_O_n  alat;
	uniform_int_distribution<int> int_dist(0,energy.size()-1);
	int curr_state,data_size=0;
	double *mag_n;
	generator.seed(seed_val);
	alat=get_bare_copy();
	alat.clear_data(true);
	data_size=energy.size();
	for(int i=0;i<data_size;i++)
	{
		curr_state=int_dist(generator);
//		cout<<","<<curr_state;
		alat.energy.push_back(energy[curr_state]);
		alat.magnetization.push_back(magnetization[curr_state]);
		mag_n = new double[_n];
		for(int j=0;j<_n;j++)
			mag_n[j]=magnetization_O_n[curr_state][j];
		alat.magnetization_O_n.push_back(mag_n);
	}
	return alat;
}
