#ifndef USR_UPDATE_H_ 
#define USR_UPDATE_H_

#include "cfl.h"
#include <math.h>

#define PI 3.14159265359

//-----Global update function
inline void snpcglO1_conf_update_D(conf_snpcglO1* conf, cfl_grid_parameter* param, double time) { 
    /* [time]=ms
     * 
     * generic start involves first setting the reference value the first time the function is called
     * for now the maximum number of reference values is 20
     * use the array to store the initial/reference value of any parameter (double!) you change
     * 
     * if(conf->init_param_update_set[0]==false) {
     * 		conf->init_param_update[0] = ( your parameter to change, e.g. conf->gG.dat[0] );
     *		conf->init_param_update_set[0] = 1;
     * }
     */
	
      if(conf->init_param_update_set[0]==false) {
      		conf->init_param_update[0] = conf->V.settings_box.box_height.dat[0];
		
     		conf->init_param_update_set[0] = 1;
      }

      double v_amp = conf->init_param_update[0]*(1-2.5/3.)*0.5;
      double f = conf->init_param_update[1];
      double t_shk = 100.;//conf->init_param_update[2];

      if (time <= t_shk)
	conf->V.settings_box.box_height.dat[0] = conf->init_param_update[0] + (1 - cos(2*PI*f*time/1000))*v_amp;
      else
	conf->V.settings_box.box_height.dat[0] = conf->init_param_update[0] + (1 - cos(2*PI*f*t_shk/1000))*v_amp;
   
}

//-----Local update function
inline void snpcglO1_conf_update_local_D(conf_snpcglO1* conf, cfl_grid_parameter* param, double time, size_t x, size_t y, size_t z) { 
    /* [time]=ms, [x_um]=[y_um]=[z_um]=um
     * double x_um = (((int)(x))-((int)(param->xD/2)))*param->Sgridding;
     * double y_um = (((int)(y))-((int)(param->yD/2)))*param->Sgridding;
     * double z_um = (((int)(z))-((int)(param->zD/2)))*param->Sgridding;
     * 
     * * if(conf->init_param_update_set[0]==false) {
     * 		conf->init_param_update[0] = ( your parameter to change, e.g. conf->gG.dat[0] );
     *		conf->init_param_update_set[0] = 1;
     * }
     */
    
    
    
}


#endif 

