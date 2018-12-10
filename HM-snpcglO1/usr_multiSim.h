#ifndef USR_MULTISIM_H_
#define USR_MULTISIM_H_

#include "cfl.h"
#define F 20.
#define STEP 2
#define TIME 25.
typedef void (*snpcglO1_multisim)(cfl_init*, conf_snpcglO1*, conf_snpcglO1*, cfl_potential*, cfl_potential*, int);

static inline void snpcglO1_multisim_D(cfl_init* I, conf_snpcglO1* conf_pre, conf_snpcglO1* conf_post, cfl_potential* V_pre, cfl_potential* V_post, int sim_index) {


  conf_pre->init_param_update[1] = F + STEP*sim_index;
  conf_post->init_param_update[1] = F + STEP*sim_index;

  double t_shake = 0.;
  double period = 1/conf_pre->init_param_update[1]*1000;

  printf("period: %f \n", period);
  while ( t_shake <= TIME){

    t_shake += period;

  }

  conf_post->init_param_update[2] = t_shake;
  printf("Shaking for: %f ms\n", t_shake);

    

 

    
}

#endif
