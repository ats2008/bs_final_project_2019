#include <random>
#include <iostream>
#include <stack> 
#include <vector> 
#include <stdlib.h>
#include <set>
#include <math.h>
#include <string>
#include <cstring>

#include <fstream>
#include <site.h>

#include "stattools.h"

#ifndef LATTICELIB_H
#define LATTICELIB_H

#define THERMALIZATION_SKIP 0
#define TAU_ERR_TOLARANCE 0.10

//#define RANDOM_ENGINE default_random_engine
#define RANDOM_ENGINE mt19937_64
#define MAX_H_LEN 10

// ISING SPECIFIC LOOKUP TABLE DEFINED (J=1 case) .. willneed change if alterd
// Defining the lattice class 

void display_site(site s);

struct lattice_write_O_n
{
		int _dimension;
		int _L;
		int _site_count;
		int _nib_count;
		char _model_type;
		double _T;
		double _beta;
		double _J;
		char _initalization;
		int _sweep_count;
		int _urd_seed;
		int _initalization_seed;
		char update_history[256];
		int _n;
		double _H[MAX_H_LEN];
};

class lattice
{
	public:
		lattice();
		vector<double> energy;
		vector<double> magnetization;
		string update_history;
	
		double get_energy();
		double get_magnetization();
		double get_temparature();
		double get_beta();
		int get_diamension();
		char get_initialization();
		int get_sweep_count();
		int get_site_count();
		int get_lattice_len();
		int get_initialization_seed();
		int get_urd_seed();
		site * get_state();
		void set_neighboures();
		void set_temparature(double Temparature);
		void set_sweep_count(int swc);
		void set_J(double J);
		void seed_generator(int seeder, int initialization_seed=-1);
		void add_energy_change(double de);
		void bootstarap_basic(double *rslt,int bstrap_count,int bstrap_idd_count);
		
		void clear_data(bool full =false);
		void clear_lattice();
		void display_lattice();
		void display_lattice_detailed();
		void display_neighbours();
		
		site* _state;
		
	protected:
		int _dimension;
		int _L;
		double _T;
		double _beta;
		double _J;
		int _nib_count;
		int _site_count;
		int _sweep_count;
		char _initalization;
		int _urd_seed;
		int _initalization_seed;
		char _model_type;
		RANDOM_ENGINE generator;
};


class lattice_O_n : public lattice
{
	public:
		lattice_O_n();
		vector<double * > magnetization_O_n;
		

		void setup_model(int dim,int n,int lattice_len,int N,int nib_count,
						double T_equilib,double J,double* H,char initialization='Z',char model_type='O');
		
		void set_ext_field(double * H);
		void get_ext_field(double * H);
		int get_n();
		string get_H_as_string();
		double get_J();
		
		int get_spin_dim();
				
		void add_magnetization_change(double * dm);
		lattice_write_O_n get_parameters();
		lattice_O_n get_bare_copy();

		double calc_energy();
		void calc_magnetization(double *magnetization_array);

		
		void clear_data(bool full=false);
		void clear_lattice();
		void display_lattice();
		void display_lattice_detailed();
		
		

		void writeToText(string fname,bool lattice_skip=false);
		void writeToBinary(string fname,bool lattice_skip=false);
		void readFromBinary(string fname);
		
		void mc_steps(int count=1);
		void wolf_steps(int count=1);
		
		lattice_O_n resample_thermalized(int tau,double tau_err,int thermalization_skip= THERMALIZATION_SKIP ,
											float tau_err_tolarance=TAU_ERR_TOLARANCE );
		lattice_O_n resample_thermalized_averages(int tau,double tau_err,int thermalization_skip= THERMALIZATION_SKIP ,
											float tau_err_tolarance=TAU_ERR_TOLARANCE );
		lattice_O_n resample_data(int indip_sample_period=1,int off_set=0);
		
		void bootstarap_mgn(double *rslt,int n,int bstrap_count,int bstrap_idd_count);
		void bootstarap_bindersC(double *rslt,int n,int bstrap_count,int bstrap_idd_count);

		lattice_O_n get_bootstrap_copy(int seed_val);
		
	protected:
		int _n;
		double *_H;
		void initialize(char mode='Z');
};


void RG_block_xy_model(lattice orig_lat,lattice &step_lat,int L_int=2) ;


#endif
