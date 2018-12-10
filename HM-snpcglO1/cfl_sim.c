/* snpcglO1 multi purpose program */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
//random numbers
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//openmp
#include <omp.h>
//libconfig
#include <libconfig.h>


#include "cfl_sim.h"

int main(int argc, char *argv[]) {
  
  //----------start time measurement
  time_t clock_start=clock();
  time_t start, end;
  double diff;
  time(&start);
  
  //----------initialise random number seed
  gsl_rng * RandGSL;
  RandGSL = gsl_rng_alloc(gsl_rng_ranlxd2);
  gsl_rng_set(RandGSL, (unsigned long int)time(NULL));
  
  
  //----------set function pointer for parameter config struct
  cfl_ran_complex complex_WN = _USR_complex_WN;
  snpcglO1_update_global f_global = _USR_f_global;
  snpcglO1_update_local f_local = _USR_f_local;
  
  
  //----------initialise config file
  cfl_init I = cfl_cfg_init("cfl_INIT.cfg");
  cfl_cfg_prepare_data_dir(&I, true, true, true);
  
  //----------runs and openMP settings
  omp_lock_t omp_lock;
  omp_init_lock(&omp_lock);
  int openmp_threads = cfl_cfg_get_int(I.cfg, "OMP_THREADS");
  int N_SIM_ALL = cfl_cfg_get_int(I.cfg, "RUNS_PER_THREAD");
  int max_threads = omp_get_max_threads();
  const int run_threads = (openmp_threads<=max_threads) ? openmp_threads : max_threads;
  
  
  //----------set imaginary time evolution solver
  snpcglO1_solver_ts			solver_ITE	= _USR_solver_ITE;
  snpcglO1_solver_kprop_ts		kprop_ITE 	= _USR_kprop_ITE;
  snpcglO1_solver_rprop_ts 		rprop_ITE 	= _USR_rprop_ITE;
  //----------set stochastic time evolution solver
  snpcglO1_solver_ts			solver_SGPE	= _USR_solver_SGPE;
  snpcglO1_solver_kprop_ts		kprop_SGPE 	= _USR_kprop_SGPE;
  snpcglO1_solver_rprop_ts 		rprop_SGPE 	= _USR_rprop_SGPE;
  //----------set real time evolution solver
  snpcglO1_solver_ts			solver_RT	= _USR_solver_RT;
  snpcglO1_solver_kprop_ts		kprop_RT 	= _USR_kprop_RT;
  snpcglO1_solver_rprop_ts 		rprop_RT 	= _USR_rprop_RT;
  
  
  //----------initialize FFTW
  cfl_fftw fftw_xtc = cfl_set_zfftw_3d(I.grid_x,CFL_FFTW_PATIENT);
  
  //----------initialize pre- and post-quench potentials
  cfl_potential V_pre, V_post;
  
  if(_USR_potential_read_file) {
      V_pre = cfl_potential_read_HDF5(_USR_potential_file_pre);
      V_post = cfl_potential_read_HDF5(_USR_potential_file_post);
  }
  else {
    V_pre = cfl_potential_init(I.grid_x, CFL_POT_DEFAULT, CFL_POT_DEFAULT, CFL_POT_DEFAULT, 0);
    V_post = cfl_potential_init(I.grid_x, CFL_POT_DEFAULT, CFL_POT_DEFAULT, CFL_POT_DEFAULT, 0);
    if(_USR_potential_harmonic[0]) {
	cfl_potential_harmonic Vh_pre = cfl_potential_harmonic_init_cfg(I.grid_x,I.cfg,0);
	cfl_potential_add(&V_pre, 1, CFL_POT_HARMONIC, Vh_pre);
    }
    if(_USR_potential_box[0]) {
	cfl_potential_box Vb_pre = cfl_potential_box_init_cfg(I.grid_x,I.cfg, 0, _USR_boxwall_type[0]);
	cfl_potential_add(&V_pre, 1, CFL_POT_BOX, Vb_pre);
    }
    if(_USR_potential_harmonic[1]) {
	cfl_potential_harmonic Vh_post = cfl_potential_harmonic_init_cfg(I.grid_x,I.cfg,1);
	cfl_potential_add(&V_post, 1, CFL_POT_HARMONIC, Vh_post);
    }
    if(_USR_potential_box[1]) {
	cfl_potential_box Vb_post = cfl_potential_box_init_cfg(I.grid_x,I.cfg, 1, _USR_boxwall_type[1]);
	cfl_potential_add(&V_post, 1, CFL_POT_BOX, Vb_post);
    }
    if(_USR_potential_ioffe[0] || _USR_potential_ioffe[1]) {
	cfl_potential_ioffe Vioffe = cfl_potential_ioffe_init_cfg(I.grid_x,I.cfg, 0);
	if(_USR_potential_ioffe[0])
	    cfl_potential_add(&V_pre, 1, CFL_POT_IOFFE, Vioffe);
	if(_USR_potential_ioffe[1])
	    cfl_potential_add(&V_post, 1, CFL_POT_IOFFE, Vioffe);
    }
    if(_USR_disorder_truncWN[0] || _USR_disorder_truncWN[1]) {
	cfl_disorder_truncWN DtruncWN = cfl_disorder_truncWN_init_cfg(I.grid_x, RandGSL, I.cfg, 0);
	if(_USR_disorder_truncWN[0])
	    cfl_potential_add(&V_pre, 1, CFL_DIS_TRUNCWN, DtruncWN);
	if(_USR_disorder_truncWN[1])
	    cfl_potential_add(&V_post, 1, CFL_DIS_TRUNCWN, DtruncWN);
    }
    
    cfl_potential_switch(_USR_potential_real[0], _USR_potential_imag[0], _USR_disorder[0], &V_pre);
    if(_USR_potential_memoryLock[0])
	cfl_potential_set_mem(_USR_potential_memoryLock_init_time[0], &V_pre, true);
    
    cfl_potential_switch(_USR_potential_real[1], _USR_potential_imag[1], _USR_disorder[1], &V_post);
    if(_USR_potential_memoryLock[1])
	cfl_potential_set_mem(_USR_potential_memoryLock_init_time[1], &V_post, true);
   
    //----------write out potential
    if(_USR_save_potential[0])
	cfl_potential_write_HDF5(&V_pre, "potential_pre.h5");
    if(_USR_save_potential[1])
	cfl_potential_write_HDF5(&V_post, "potential_post.h5");
  }
  
  
  //----------ser pre- and post-quench parameter config structs
  conf_snpcglO1 conf_Mpre = conf_snpcglO1_init_cfg(V_pre,I,0);
  conf_snpcglO1_init_ran_complex(&conf_Mpre, complex_WN);
  conf_snpcglO1_init_update_conf(&conf_Mpre, f_global, f_local);
  conf_snpcglO1 conf_Mpost = conf_snpcglO1_init_cfg(V_post,I,1);
  conf_snpcglO1_init_ran_complex(&conf_Mpost, complex_WN);
  conf_snpcglO1_init_update_conf(&conf_Mpost, f_global, f_local);
  
  
  //----------print some simulation info
  printf("OpenMP info: Maximum number of cores available: %i\n",max_threads);
  printf("             Number of threads for run-loop   : %i\n",run_threads);
  printf("Starting simulation with:   Timestep:    %.4f ms\n",(I.Tgridding*conf_Mpost.step*1000));
  printf("                            Observables: %.4f ms\n",(I.Tgridding*conf_Mpost.step*conf_Mpost.T_evolve*1000));
  printf("                            Total time:  %.4f ms\n",(I.Tgridding*conf_Mpost.step*conf_Mpost.T_evolve*conf_Mpost.N_evolve*1000));
  
  //----------array for SGPE/RT_data output
  cfl_grid_parameter grid_out_SGPE = cfl_set_parameter_gtcr(I.grid_x,1,I.grid_x.component,(run_threads*N_SIM_ALL));
  cfl_complex_grid* data_out_SGPE = NULL;
  cfl_complex_grid* dataK_out_SGPE = NULL;
  if(_USR_doSGPE) {
      data_out_SGPE = cfl_CG_allocate(grid_out_SGPE);
      if(_USR_save_kgrid_SGPE)
	  dataK_out_SGPE = cfl_CG_allocate(grid_out_SGPE);
  }
  cfl_grid_parameter grid_out = cfl_set_parameter_gtcr(I.grid_x,(conf_Mpost.N_evolve+1),I.grid_x.component,(run_threads*N_SIM_ALL));
  cfl_complex_grid* data_out_RT = NULL;
  cfl_complex_grid* dataK_out_RT = NULL; 
  if(_USR_doRT) {
      data_out_RT = cfl_CG_allocate(grid_out);
      if(_USR_save_kgrid_RT)
	  dataK_out_RT = cfl_CG_allocate(grid_out);
  }
      
  int multiSim_loops = (_USR_domultiSim) ? _USR_multiSim_loops : 1;
  for(int mSim=0; mSim<multiSim_loops; mSim++) {
      
      char filename_sim[128];
      if(multiSim_loops==1)
	  sprintf(filename_sim, "simDATA.h5");
      else {
	  printf("---------------------------\n");
	  printf("multiSim: %i\n",mSim);
	  sprintf(filename_sim, "%imSim_simDATA.h5",mSim);
      }
      
      //----------set data_out arrays to zero
      if(data_out_SGPE!=NULL) { cfl_CGv_zero(data_out_SGPE); }
      if(dataK_out_SGPE!=NULL) { cfl_CGv_zero(dataK_out_SGPE); }
      if(data_out_RT!=NULL) { cfl_CGv_zero(data_out_RT); }
      if(dataK_out_RT!=NULL) { cfl_CGv_zero(dataK_out_RT); }
      
      //----------write physical dimensions for grid
      cfl_hdf5_system_parameter(filename_sim, I, &V_pre, &V_post);
      
      //----------change parameter in multiSim
      (*_USR_multiSim)(&I, &conf_Mpre, &conf_Mpost, &V_pre, &V_post, mSim);
      
      //----------imaginary time evolution
      cfl_complex_grid* rgrid = cfl_CG_allocate(I.grid_x);
      cfl_complex_grid* kgrid = cfl_CG_allocate(I.grid_x);
      
      //----------read in T=0 single core initial data
      if(_USR_read_data_1) {
	  cfl_hdf5_read_complex(rgrid, _USR_read_data_1_file, _USR_read_data_1_dset);
	  cfl_zfftw_execute(&fftw_xtc, rgrid, kgrid, CFl_FFTW_FORWARD);
	  cfl_CGv_zdscal(kgrid,fftw_xtc.norm);
      }
      
      if(_USR_doITE) {
	  conf_snpcglO1 conf_ITE = conf_snpcglO1_set_mask(&conf_Mpre,_USR_mask_ITE);
	  (*kprop_ITE)(rgrid->param, &conf_ITE);  
	  _Bool ITE_conv = snpcglO1_ts_fourier_ITE(solver_ITE, rprop_ITE, RandGSL, I.T_evolve_ITE, 
						   rgrid, kgrid, &fftw_xtc, &conf_ITE, I.convergence_ITE, 
					    I.max_evolve_ITE, _USR_add_const_ITE, _USR_save_steps_ITE, _USR_save_filename_ITE,_USR_update_ITE);
	  if(!ITE_conv) {
	      fprintf(stderr,"ERROR: ITE is not converged!\n");
	      cfl_CG_free(kgrid);
	      cfl_CG_free(rgrid);
	      conf_snpcglO1_free(&conf_Mpost);
	      conf_snpcglO1_free(&conf_Mpre);
	      cfl_free_fftw(&fftw_xtc,true);
	      cfl_potential_free(&V_post);
	      cfl_potential_free(&V_pre);
	      gsl_rng_free(RandGSL);
	      cfl_cfg_free(&I, false);
	      exit(EXIT_FAILURE);
	  }
	  conf_snpcglO1_free(&conf_ITE);
      }
      
      //----------Calculate chemical potential and set in conf_Mpre
      if(_USR_doCHEMPOT && !_USR_setCHEMPOT) {
	  double cp = snpcglO1_chempot_tsfp_set(rgrid, kgrid, &conf_Mpre,_USR_mask_CHEMPOT,_USR_update_CHEMPOT);
	  printf("chemical potential:  %e [kHz]\n",cp/(1000.*2.*M_PI*I.Tgridding));
      }
      if(_USR_setCHEMPOT && !_USR_doCHEMPOT) {
	  double chem_pot = _USR_chempot;
	  if(gsl_complex_abs2(conf_Mpre.cp)==0.)
	      conf_Mpre.cp = gsl_complex_rect(chem_pot,(2.*M_PI*I.Tgridding*chem_pot));
	  else
	      gsl_complex_mul_real(conf_Mpre.cp,chem_pot);
      }
      
      cfl_hdf5_write_complex(rgrid, filename_sim, "INIT_rgrid", -1);
      if(_USR_save_kgrid_ITE)
	  cfl_hdf5_write_complex(kgrid, filename_sim, "INIT_kgrid", -1);
      cfl_CG_free(kgrid);
      cfl_CG_free(rgrid);
      
      
      if(_USR_doSGPE || _USR_doRT) {
	  #pragma omp parallel for schedule(guided) num_threads(run_threads)
	  for(int run=0; run<run_threads; run++) {
	      
	      omp_set_lock(&omp_lock);
	      gsl_rng* RandGSL_thread;
	      RandGSL_thread = gsl_rng_alloc(gsl_rng_ranlxd2);
	      int rand_offset = gsl_ran_flat(RandGSL, 0, 1000000);
	      gsl_rng_set(RandGSL_thread, (unsigned long int)(time(NULL)+rand_offset+run));
	      
	      conf_snpcglO1 conf_SGPE = conf_snpcglO1_set_mask(&conf_Mpre,_USR_mask_SGPE);
	      (*kprop_SGPE)(I.grid_x, &conf_SGPE);
	      
	      conf_snpcglO1 conf_RT = conf_snpcglO1_set_mask(&conf_Mpost,_USR_mask_RT);
	      (*kprop_RT)(I.grid_x, &conf_RT);
	      
	      cfl_complex_grid* rgrid_init = cfl_CG_allocate(I.grid_x);
	      cfl_hdf5_read_complex(rgrid_init, filename_sim, "INIT_rgrid");
	      cfl_complex_grid* rgrid = cfl_CG_allocate(I.grid_x);
	      cfl_complex_grid* kgrid = cfl_CG_allocate(I.grid_x);
	      
	      char filename_coreIN[128];
	      if(multiSim_loops==1)
		  sprintf(filename_coreIN, "data_in/%s_%i.h5",_USR_read_data_2_file,run);
	      else
		  sprintf(filename_coreIN, "data_in/%imSim_%s_%i.h5",mSim,_USR_read_data_2_file,run);
	      
	      omp_unset_lock(&omp_lock);
	      
	      for(int run_t=0; run_t<N_SIM_ALL; run_t++) {
		  
		  if(run==0) { printf("%i ",run_t); fflush(stdout); }
		  
		  conf_snpcglO1_reset_time(&conf_SGPE);
		  conf_snpcglO1_reset_time(&conf_RT);
		  cfl_CGv_zero(rgrid);
		  cfl_CGv_zero(kgrid);
		  
		  //----------T=0 ITE data
		  cfl_CGv_copy(rgrid,rgrid_init);
		  
		  //----------read in data for SGPE/RT evolution (multiple cores)
		  if(_USR_read_data_2) {
		      omp_set_lock(&omp_lock);
		      char dset[128];
		      sprintf(dset,"%i_%s",run_t,_USR_read_data_2_dset);
		      cfl_hdf5_read_complex(rgrid, filename_coreIN, dset);
		      omp_unset_lock(&omp_lock);
		      cfl_zfftw_execute(&fftw_xtc, rgrid, kgrid, CFl_FFTW_FORWARD);
		      cfl_CGv_zdscal(kgrid,fftw_xtc.norm);
		  }
		  
		  //----------SGPE growth evolution for finite temperature initial state with conf_Mpre
		  if(_USR_doSGPE) {
		      char filename_SGPE[128];
		      if(_USR_save_steps_SGPE)
			  sprintf(filename_SGPE, "%s_%i_%i.h5",_USR_save_filename_SGPE,run,run_t);
		      
		      _Bool SGPE_conv = snpcglO1_ts_fourier_SGPE(solver_SGPE, rprop_SGPE, RandGSL_thread, I.T_evolve_SGPE, 
								 rgrid, kgrid, &fftw_xtc, &conf_SGPE, I.convergence_SGPE, 
						   I.conv_evolve_SGPE, I.max_evolve_SGPE, _USR_init_zero_SGPE, 
						   _USR_save_steps_SGPE, filename_SGPE, _USR_update_SGPE);
		      if(!SGPE_conv)
			  printf("core: %i , run: %i SGPE convergence failed!\n",run,run_t);
		      
		      cfl_MRCG_xc2MR(data_out_SGPE,rgrid,0,(run*N_SIM_ALL+run_t));
		      if(_USR_save_kgrid_SGPE)
			  cfl_MRCG_xc2MR(dataK_out_SGPE,kgrid,0,(run*N_SIM_ALL+run_t));
		  }
		  
		  //----------edit data before real time evolution
		  if(_USR_edit_data) {
		      (*_USR_chdata)(rgrid,kgrid,&conf_Mpost);
		      cfl_zfftw_execute(&fftw_xtc, rgrid, kgrid, CFl_FFTW_FORWARD);
		      cfl_CGv_zdscal(kgrid,fftw_xtc.norm);
		  }
		  
		  //----------real time evolution with conf_Mpost
		  if(_USR_doRT) {
		      for(int t=0; t<conf_RT.N_evolve; t++) {
			  
			  cfl_MRCG_xc2MR(data_out_RT,rgrid,t,(run*N_SIM_ALL+run_t));
			  if(_USR_save_kgrid_RT)
			      cfl_MRCG_xc2MR(dataK_out_RT,kgrid,t,(run*N_SIM_ALL+run_t));
			  
			  (*solver_RT)(rprop_RT, RandGSL_thread, conf_RT.T_evolve, rgrid, kgrid, &fftw_xtc, &conf_RT, _USR_update_RT);
			  
		      }
		      
		      cfl_MRCG_xc2MR(data_out_RT,rgrid,conf_RT.N_evolve,(run*N_SIM_ALL+run_t));
		      if(_USR_save_kgrid_RT)
			  cfl_MRCG_xc2MR(dataK_out_RT,kgrid,conf_RT.N_evolve,(run*N_SIM_ALL+run_t));
		      
		  }
		  
	      }
	      cfl_CG_free(kgrid);
	      cfl_CG_free(rgrid);
	      cfl_CG_free(rgrid_init);
	      conf_snpcglO1_free(&conf_RT);
	      conf_snpcglO1_free(&conf_SGPE);
	      gsl_rng_free(RandGSL_thread);
	  }
	  
      }
      printf("\n");
      
      if(_USR_doSGPE) {
	  cfl_hdf5_write_complex(data_out_SGPE, filename_sim, "SGPE_rgrid", -1);
	  if(_USR_save_kgrid_SGPE)
	      cfl_hdf5_write_complex(dataK_out_SGPE, filename_sim, "SGPE_kgrid", -1);
      }
      if(_USR_doRT) {
	  cfl_hdf5_write_complex(data_out_RT, filename_sim, "RT_rgrid", -1);
	  if(_USR_save_kgrid_RT)
	      cfl_hdf5_write_complex(dataK_out_RT, filename_sim, "RT_kgrid", -1);
      }
      
  }
  
  //---------end time measurement
  time(&end);
  diff=difftime(end, start);
  double time_c=(double)(clock()- clock_start)/CLOCKS_PER_SEC;
  
  printf("\n");
  printf("Simulation finished in %.2f seconds (from clocks per sec) and %.2f seconds (from time difference)\n",time_c,diff);
  printf("Clearing memory...\n\n");
  
  if(data_out_SGPE!=NULL) { cfl_CG_free(data_out_SGPE); }
  if(dataK_out_SGPE!=NULL) { cfl_CG_free(dataK_out_SGPE); }
  if(data_out_RT!=NULL) { cfl_CG_free(data_out_RT); }
  if(dataK_out_RT!=NULL) { cfl_CG_free(dataK_out_RT); }
  conf_snpcglO1_free(&conf_Mpost);
  conf_snpcglO1_free(&conf_Mpre);
  cfl_free_fftw(&fftw_xtc,true);
  cfl_potential_free(&V_post);
  cfl_potential_free(&V_pre);  
  gsl_rng_free(RandGSL);
  cfl_cfg_free(&I, false);
  
  return 0;
}
 

