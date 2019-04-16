#include <random>
#include <iostream>
#include <stack> 
#include <vector> 
#include <stdlib.h>
#include <set>
#include <string>
#include <cstring>
#include <fstream>

#ifndef LATTICE_H
#define LATTICE_H


#define RANDOM_ENGINE default_random_engine


// ISING SPECIFIC LOOK UPTABLE DEFINED (J=1 case) .. willneed change if alterd
// Defining the lattice class 

typedef  int spin;
using namespace std;

class site{
	public:
		spin _spin;
		site **_neighbour;
		site()
		{
			_neighbour=NULL;
		}
		
};
void display_site(site s);

struct lattice_write
{
		double _T;
		double _beta;
		double _J;
		double _wolf_add_probab;
		double _wolf_add_probab_potts;
		int _dimension;
		int _site_count;
		int _q;
		int _L;
		char _initalization;
		int _sweep_count;
		int _urd_seed;
		int _initalization_seed;
		char update_history[2096];
};


class lattice
{
	public:
		lattice();
		void clear_lattice();
		lattice(int dim,int lattice_len,int N,double T_equilib,double J=1.0,int spin_count=2,char initialization='Z',char model_type='I');
		vector<double> energy;
		vector<double> magnetization;
		string update_history;
		
		double get_energy();
		double get_magnetization();
		double get_temparature();
		char get_initialization();
		int get_sweep_count();
		int get_site_count();
		int get_lattice_len();
		site * get_state();
		
		void set_neighboures();
		void set_temparature(double Temparature);
		void clear_data();
		void set_sweep_count(int swc);
		void seed_generator(int seeder, int initialization_seed=-1);
		void set_LUT(int sp_count=2,int J=1);
		
		void add_new_energy(double de);
		void add_new_magnetization(double dm);
		void add_energy(double e);
		void add_magnetization(double m);
		double calc_energy();
		double calc_potts_energy();
		double calc_magnetization();

		void display_lattice();
		lattice_write get_parameters();
		void display_lattice_detailed();
		void writeToText(string fname);
		
		void mc_steps(int count=1);
		void wolf_steps(int count=1);
		void mc_potts_steps(int count=1);
		void wolf_potts_steps(int count=1);

	private:

		site* _state;
		double *_LUT;
		double _T;
		double _beta;
		double _J;
		double _wolf_add_probab;
		double _wolf_add_probab_potts;
		int _dimension;
		int _site_count;
		int _q;
		int _L;
		int _sweep_count;
		char _initalization;
		int _urd_seed;
		int _initalization_seed;
		char _model_type;
		RANDOM_ENGINE generator;
		double get_exp(double de);
		void initialize(char mode='Z');
		void ising_standardize();

};


#endif
