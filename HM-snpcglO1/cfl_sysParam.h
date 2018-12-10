#ifndef CFL_SYSPARAM_H_ 
#define CFL_SYSPARAM_H_

#include "cfl.h"

double cfl_hdf5_system_parameter_snpcglO1(hid_t group_id, const char* dset_name, cfl_init I, int cfg_list_entry, int dimension, double hbar_scale, double as_scale, double mass);
double cfl_hdf5_system_parameter_snpcglO2(hid_t group_id, const char* dset_name, cfl_init I, int cfg_list_entry, int dimension, double hbar_scale, double as_scale, double mass);

void cfl_hdf5_system_parameter(const char* filename, cfl_init I, cfl_potential* V_pre, cfl_potential* V_post);

#endif
