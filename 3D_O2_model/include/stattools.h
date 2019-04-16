 #ifndef _STATOOLS_
 #define _STATOOLS_
 
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <sys/stat.h>
#include "latticelib.h"


#include <gsl/gsl_statistics.h>

#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace boost::multiprecision;
#define Z_MAX_ITR 60
#define Z_TOLARENCE 1e-5

double average(double *array,int len);
double variance(double *array,int len);
cpp_dec_float_100 average(cpp_dec_float_100 *array,int len);
cpp_dec_float_100 variance(cpp_dec_float_100 *array,int len);


long double* resolve_Z(vector<lattice_O_n> data_vector,long double *Z0=NULL,int MAX_ITR=Z_MAX_ITR,double TOLARENCE=Z_TOLARENCE);
double get_expectation_value(double desiered_J,double* desiered_H,double (*func)(double E,double *M),
								vector<lattice_O_n> data_vector,long double * Z=NULL);
long double get_partition_function(double desiered_J,double* desiered_H,vector<lattice_O_n> data_vector,long double * Z=NULL);


cpp_dec_float_100 * resolve_Z_large(vector<lattice_O_n> data_vector,cpp_dec_float_100 *Z0=NULL,string Z_file="",int MAX_ITR=Z_MAX_ITR,double TOLARENCE=Z_TOLARENCE);
cpp_dec_float_100 get_partition_function_large(double desiered_J,double* desiered_H,vector<lattice_O_n> data_vector,long double * Z=NULL);
cpp_dec_float_100 get_expectation_value_large(double desiered_J,double* desiered_H,double (*func)(double E,double *M),
								vector<lattice_O_n> data_vector,cpp_dec_float_100* Z=NULL);

void save_Z(cpp_dec_float_100 *Z,int N,string fname);
cpp_dec_float_100 * read_Z(int *N,string fname);

#endif
