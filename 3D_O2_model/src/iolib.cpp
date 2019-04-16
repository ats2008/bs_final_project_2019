#include "iolib.h"

void WriteLatticeToBinary(string fname,lattice & alattice)
{
	fstream ofile;
	ofile.open(fname.c_str(),ios::out | ios::binary);
	lattice_write parameters;
	parameters=alattice.get_parameters();

	ofile.write((char *)(&parameters),sizeof(lattice_write));
	int site_count=alattice.get_site_count();
	site * state;
	state=alattice.get_state();
	for(int i=0;i<site_count;i++)
		ofile.write((char *)(&state[i]._spin),sizeof(spin));

	vector<double>::iterator en_it,mg_it; 
	en_it = alattice.energy.begin();
	mg_it = alattice.magnetization.begin();
	double d;
	while(en_it!=alattice.energy.end() and mg_it!=alattice.magnetization.end())
	{
		d=*en_it;
		ofile.write((char *)(&d),sizeof(d));
		//cout<<"writing th E,M = "<<d<<",";
		d=*mg_it;
		ofile.write((char *)(&d),sizeof(d));
		//cout<<d<<" \n";
		mg_it++;
		en_it++;
	}
	ofile.close();
}

lattice ReadLattticeFromBinary(string fname)
{
	fstream infile;
	infile.open(fname.c_str(),ios::in | ios::binary);
	if(!infile.is_open())
	{
		cout<<"\n ERROR HAPPEND !! FILE DOES NOT EXIST !! \n";
		cout<<"fname : "<<fname;
		exit(0);
	}
	lattice_write parameters;
	
	infile.read((char *)(&parameters),sizeof(lattice_write));

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
	site * state;
	state=blattice.get_state();
	for(int i=0;i<parameters._site_count;i++)
		{
			infile.read((char *)(&state[i]._spin),sizeof(spin));
		}
	blattice.set_neighboures();

	int i=0;
	double d=0;
	blattice.energy.clear();
	blattice.magnetization.clear();
	infile.read((char *)(&d),sizeof(d));
	while(!infile.eof())
	{	
		i+=1;
		blattice.add_energy(d);
	//	cout<<"got "<<i<<" th - > E = "<<d;
		infile.read((char *)(&d),sizeof(d));
		blattice.add_magnetization(fabs(d));
	//	cout<<"   M = "<<d<<"  tellg = "<<infile.tellg()<<"\n"; 
		infile.read((char *)(&d),sizeof(d));
	}
//	cout<<"\n E TOP : "<<blattice.energy.back();
//	cout<<"\n M TOP : "<<blattice.magnetization.back();
	blattice.set_sweep_count(i-1);
	infile.close();
	return blattice;
}


///////////////////////////////////  XY Model ////////////////////////////////

void WriteLatticeToBinary_xyModel(string fname,lattice alattice,bool lattice_skip)
{
//	cout<<"\n here \n";
	fstream ofile;
	ofile.open(fname.c_str(),ios::out | ios::binary);
//	cout<<"\n fname = "<<fname<<"\n";
	lattice_write parameters;
	parameters=alattice.get_parameters_xyModel();
	if(lattice_skip) parameters._site_count=0;
	ofile.write((char *)(&parameters),sizeof(parameters));

	int site_count;
//	cout<<"wrote struct \n";
	site_count=parameters._site_count;
	site * state;
	state=alattice.get_state();
	for(int i=0;i<site_count;i++)
		ofile.write((char *)(&state[i]._spin),sizeof(spin));

	vector<double>::iterator en_it,mg_it,mgx_it,mgy_it; 
	en_it = alattice.energy.begin();
	mg_it = alattice.magnetization.begin();
	mgx_it = alattice.magnetization_x.begin();
	mgy_it = alattice.magnetization_y.begin();
	double d=0;
//	cout<<"Number of elts : "<<alattice.energy.size()<<"\n";
	while(en_it!=alattice.energy.end() and mg_it!=alattice.magnetization.end() and 
								 mgx_it!=alattice.magnetization_x.end() and mgy_it!=alattice.magnetization_y.end())
	{
		d=*en_it;
		ofile.write((char *)(&d),sizeof(double));
//		ofile<<d;
		d=*mg_it;
		ofile.write((char *)(&d),sizeof(d));
		d=*mgx_it;
		ofile.write((char *)(&d),sizeof(d));
		d=*mgy_it;
		ofile.write((char *)(&d),sizeof(d));
		mg_it++;
		mgx_it++;
		mgy_it++;
		en_it++;
	}
	ofile.close();
}

lattice ReadLattticeFromBinary_xyModel(string fname)
{
	fstream infile;
	infile.open(fname.c_str(),ios::in | ios::binary);
	if(!infile.is_open())
	{
		cout<<"\n ERROR HAPPEND !! FILE DOES NOT EXIST !! \n";
		cout<<"fname : "<<fname;
		exit(0);
	}
	lattice_write parameters;
	
	infile.read((char *)(&parameters),sizeof(lattice_write));
	lattice blattice;
	//alat.setup_xy_model(dim,L,N,T,J,0.0,0.0,initialization)
	blattice.setup_xy_model(parameters._dimension,parameters._L,parameters._site_count,parameters._T,parameters._J,parameters._h4,
													parameters._h8,parameters._initalization,'X');
	blattice.update_history=parameters.update_history;
	blattice.seed_generator(parameters._urd_seed,parameters._initalization_seed);
	site * state;
	state=blattice.get_state();
	for(int i=0;i<parameters._site_count;i++)
		{
			infile.read((char *)(&state[i]._spin),sizeof(spin));
		}
	blattice.set_neighboures();

	int i=0;
	double d=0;
	blattice.energy.clear();
	blattice.magnetization.clear();
	blattice.magnetization_x.clear();
	blattice.magnetization_y.clear();
	infile.read((char *)(&d),sizeof(d));
	while(!infile.eof())
	{	
		i+=1;
		blattice.energy.push_back(d);
		
		infile.read((char *)(&d),sizeof(d));
		blattice.magnetization.push_back(d);
		
		infile.read((char *)(&d),sizeof(d));
		blattice.magnetization_x.push_back(d);

		infile.read((char *)(&d),sizeof(d));
		blattice.magnetization_y.push_back(d);
		
		infile.read((char *)(&d),sizeof(d));
	}
	blattice.set_sweep_count(i-1);
	infile.close();
	return blattice;
}








