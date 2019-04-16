 #ifndef _STATOOLS_
 #define _STATOOLS_
 
 #define FIXED_POINT_MAX_ITR  200
 #define FIXED_POINT_MAX_TOLARENCE 0.000001
 #define THERMALIZATION_SKIP 3000
 #define TAU_ERR_TOLARANCE 0.10
 #define DEFAULT_END_COUNT 10000
 
 #include "isinglib.h"
 #include <math.h>
 #include <iostream>
 #include <vector>
 using namespace std;
 
 
long double _logadd(long double l1,long double l2);
long double* resolve_Z(vector<lattice> data_vector,long double *Z0=NULL,int MAX_ITR=FIXED_POINT_MAX_ITR,double tol=FIXED_POINT_MAX_TOLARENCE);
long double get_partition_function(double desiered_T,vector<lattice> data_vector,long double * Z = NULL);
double get_expectation_value(double desiered_T,double (*func)(double E,double M),vector<lattice> data_vector,long double * Z=NULL);
lattice resample_thermalized(lattice alat,int tau,double tau_err,int thermalization_skip=THERMALIZATION_SKIP,float tau_err_tolarance=TAU_ERR_TOLARANCE);
lattice resample_end(lattice alat,int count = DEFAULT_END_COUNT);
 #endif
