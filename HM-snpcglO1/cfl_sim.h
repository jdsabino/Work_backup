#ifndef CFL_SIM_H_
#define CFL_SIM_H_

#include "cfl.h"
#include "cfl_sysParam.h"
#include "usr_update.h"
#include "usr_multiSim.h"
#include "usr_chData.h"

cfl_ran_complex _USR_complex_WN                 = &cfl_ptools_complex_gaussian_WN_ziggurat;
snpcglO1_update_global _USR_f_global            = &snpcglO1_conf_update_D;
snpcglO1_update_local _USR_f_local              = &snpcglO1_conf_update_local_D;

//----------multiSim (e.g. for parameter sweep)
_Bool _USR_domultiSim                           = true;
//________________________________________________________________________________________________//
int _USR_multiSim_loops                         = 105;
snpcglO1_multisim _USR_multiSim                 = &snpcglO1_multisim_D;
//________________________________________________________________________________________________//




//----------initialise potentials
_Bool _USR_potential_read_file                  = false;
char* _USR_potential_file_pre                   = "data_in/potential_pre.h5";
char* _USR_potential_file_post                  = "data_in/potential_post.h5";
//________________________________________________________________________________________________//
_Bool _USR_potential_harmonic[2]                = {true , true};

_Bool _USR_potential_box[2]                     = {true , true};
cfl_potential_box_walltype _USR_boxwall_type[2] = {CFL_BOX_WALL_ERF , CFL_BOX_WALL_ERF};

_Bool _USR_potential_ioffe[2]                   = {false , false};

_Bool _USR_disorder_truncWN[2]                  = {false , false};

cfl_potential_type _USR_potential_real[2]       = {CFL_POT_BOX , CFL_POT_BOX};
cfl_potential_type _USR_potential_imag[2]       = {CFL_POT_DEFAULT , CFL_POT_DEFAULT};
cfl_potential_type _USR_disorder[2]             = {CFL_POT_DEFAULT , CFL_POT_DEFAULT};

_Bool _USR_potential_memoryLock[2]              = {false, false};
double _USR_potential_memoryLock_init_time[2]   = {0. , 0.};

_Bool _USR_save_potential[2]                    = {false , false};
//________________________________________________________________________________________________//



/*----------read in T=0 data 
 * single core 
 * one file with one complex dataset in shape (c,t,x,y,z) is read in
 */
_Bool _USR_read_data_1                          = false;
//________________________________________________________________________________________________//
char* _USR_read_data_1_file                     = "data_in/init_data_1.h5";
char* _USR_read_data_1_dset                     = "/dset_in";
//________________________________________________________________________________________________//



//----------imaginary time evolution
_Bool _USR_doITE                                = true;
//________________________________________________________________________________________________//
_Bool _USR_save_kgrid_ITE                       = true;
_Bool _USR_add_const_ITE                        = true;
_Bool _USR_save_steps_ITE                       = false;                                          //only use outside multithreading!
char* _USR_save_filename_ITE                    = "ITE.h5";                                       //only use outside multithreading!
conf_snpcglO1_mask _USR_mask_ITE                = SNPCGLO1_IGPE;
conf_snpcglO1_update _USR_update_ITE            = SNPCGLO1_STATIC;
snpcglO1_solver_ts _USR_solver_ITE              = &snpcglO1_propagate_tsfp_O1;
snpcglO1_solver_kprop_ts _USR_kprop_ITE         = &snpcglO1_kprop_init_tsfp_O1;
snpcglO1_solver_rprop_ts _USR_rprop_ITE         = &snpcglO1_rprop_ts_fourier_norm_constant;
//________________________________________________________________________________________________//



//----------chemical potential
_Bool _USR_doCHEMPOT                            = true;
//________________________________________________________________________________________________//
conf_snpcglO1_mask _USR_mask_CHEMPOT            = SNPCGLO1_GPE;
conf_snpcglO1_update _USR_update_CHEMPOT        = SNPCGLO1_STATIC;
_Bool _USR_setCHEMPOT                           = false;
double _USR_chempot                             = 0.;                                             //in [Hz]
//________________________________________________________________________________________________//



/*----------read in real time data 
 * parallel (multi-core)
 * each core has its own hdf5 file XXX_(#core number).h5
 * each file has (runs per core) datasets in the form (#run of core)_XXX
 * T=0 data from imaginary time evolution is loaded by default
 */
_Bool _USR_read_data_2                          = false;
//________________________________________________________________________________________________//
char* _USR_read_data_2_file                     = "coreDATA";
char* _USR_read_data_2_dset                     = "SGPE_rgrid";
//________________________________________________________________________________________________//



//----------stochastic time evolution
_Bool _USR_doSGPE                               = false;
//________________________________________________________________________________________________//
_Bool _USR_save_kgrid_SGPE                      = true;
_Bool _USR_init_zero_SGPE                       = false;
_Bool _USR_save_steps_SGPE                      = false;                                          //only use outside multithreading!
char* _USR_save_filename_SGPE                   = "SGPE";                                         //only use outside multithreading!
conf_snpcglO1_mask _USR_mask_SGPE               = SNPCGLO1_SGPE;
conf_snpcglO1_update _USR_update_SGPE           = SNPCGLO1_STATIC;
snpcglO1_solver_ts _USR_solver_SGPE             = &snpcglO1_propagate_tsfp_O1;
snpcglO1_solver_kprop_ts _USR_kprop_SGPE        = &snpcglO1_kprop_init_tsfp_O1;
snpcglO1_solver_rprop_ts _USR_rprop_SGPE        = &snpcglO1_rprop_ts_fourier_norm_constant;
//________________________________________________________________________________________________//



//----------edit data before real time evolution
_Bool _USR_edit_data                            = true;
//________________________________________________________________________________________________//
snpcglO1_chData _USR_chdata                     = &snpcglO1_changedata_D;
//________________________________________________________________________________________________//



//----------set real time evolution solver
_Bool _USR_doRT                                 = true;
//________________________________________________________________________________________________//
_Bool _USR_save_kgrid_RT                        = true;
conf_snpcglO1_mask _USR_mask_RT                 = SNPCGLO1_GPE;
conf_snpcglO1_update _USR_update_RT             = SNPCGLO1_UPDATE;
snpcglO1_solver_ts _USR_solver_RT               = &snpcglO1_propagate_tsfp_O2;
snpcglO1_solver_kprop_ts _USR_kprop_RT          = &snpcglO1_kprop_init_tsfp_O2;
snpcglO1_solver_rprop_ts _USR_rprop_RT          = &snpcglO1_rprop_ts_fourier_norm_constant;
//________________________________________________________________________________________________//

#endif




