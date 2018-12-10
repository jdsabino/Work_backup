#include "cfl_sysParam.h"

double cfl_hdf5_system_parameter_snpcglO1(hid_t group_id, const char* dset_name, cfl_init I, int cfg_list_entry, int dimension, double hbar_scale, double as_scale, double mass) {
    
    int Rank = 1;
    hsize_t dims[Rank];
    dims[0] = 6;
    hid_t DSPACE = H5Screate_simple(Rank, dims, NULL);
	
    config_setting_t *setting_pre;
    setting_pre = config_lookup(&I.cfg, "snpcgl_propagator");
    config_setting_t *prop_pre = config_setting_get_elem(setting_pre, cfg_list_entry);
    double N_pre = I.particles;
    double w_perp_pre = (2.*M_PI*cfl_cfg_getE_float(prop_pre, "V_PERP",-1));
    double a_perp_pre = sqrt(hbar_scale/(mass*w_perp_pre));
    double T_pre = cfl_cfg_getE_float(prop_pre, "TEMPERATURE",-1);
    double J_pre = 0.;
    double g_pre;
    switch(dimension) {
	case 1:
	    g_pre = 2.*hbar_scale*w_perp_pre*as_scale;
	    break;
	case 2:
	    if(I.zD>1)
		g_pre = (2.*sqrt(2.*M_PI)*((hbar_scale*hbar_scale)/mass)*(as_scale/a_perp_pre));
	    else
		g_pre = (4.*M_PI*((hbar_scale*hbar_scale)/mass)*(as_scale/I.Sgridding));
	    break;
	case 3:
	    g_pre = (4.*M_PI*((hbar_scale*hbar_scale)/mass)*as_scale);
	    break;
	default:
	    fprintf(stderr,"ERROR #hdf5_system_parameter: Physical units not defined for D>3.\n");
	    exit(EXIT_FAILURE);
	    break;
    }
    double sys_param_pre[6] = {N_pre,w_perp_pre,a_perp_pre,T_pre,J_pre,g_pre};
    
    hid_t DSET = H5Dcreate2(group_id, dset_name, H5T_NATIVE_DOUBLE, DSPACE, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(DSET, H5T_NATIVE_DOUBLE, DSPACE, DSPACE, H5P_DEFAULT, &sys_param_pre);
    H5Dclose(DSET);
    H5Sclose(DSPACE);
    
    double step = cfl_cfg_getE_float(prop_pre, "TIMESTEP",-1);
    return step;
}

double cfl_hdf5_system_parameter_snpcglO2(hid_t group_id, const char* dset_name, cfl_init I, int cfg_list_entry, int dimension, double hbar_scale, double as_scale, double mass) {
    
    int Rank = 2;
    hsize_t dims[Rank];
    dims[0] = 2; dims[1] = 6;
    hid_t DSPACE = H5Screate_simple(Rank, dims, NULL);
    
    config_setting_t *setting_post;
    setting_post = config_lookup(&I.cfg, "snpcgl_o2_propagator");
    config_setting_t *prop_post = config_setting_get_elem(setting_post, cfg_list_entry);
    double N_post[2], w_perp_post[2], a_perp_post[2], T_post[2], g_post[2];
    for(int i=0; i<2; i++) {
	N_post[i] = I.particles*cfl_cfg_getE_float(prop_post, "PARTICLE_IMBALANCE",i); 
	w_perp_post[i] = (2.*M_PI*cfl_cfg_getE_float(prop_post, "V_PERP",i));
	a_perp_post[i] = sqrt(hbar_scale/(mass*w_perp_post[i]));
	T_post[i] = cfl_cfg_getE_float(prop_post, "TEMPERATURE",i);
	switch(dimension) {
	    case 1:
		g_post[i] = 2.*hbar_scale*w_perp_post[i]*as_scale;
		break;
	    case 2:
		if(I.zD>1)
		    g_post[i] = (2.*sqrt(2.*M_PI)*((hbar_scale*hbar_scale)/mass)*(as_scale/a_perp_post[i]));
		else
		    g_post[i] = (4.*M_PI*((hbar_scale*hbar_scale)/mass)*(as_scale/I.Sgridding));
		break;
	    case 3:
		g_post[i] = (4.*M_PI*((hbar_scale*hbar_scale)/mass)*as_scale);
		break;
	    default:
		fprintf(stderr,"ERROR #hdf5_system_parameter: Physical units not defined for D>3.\n");
		exit(EXIT_FAILURE);
		break;
	}
    }
    double J_post = (2.*M_PI*cfl_cfg_getE_float(prop_post, "JJ_RE",-1));
    double sys_param_post[12] = {N_post[0],w_perp_post[0],a_perp_post[0],T_post[0],J_post,g_post[0],N_post[1],w_perp_post[1],a_perp_post[1],T_post[1],J_post,g_post[1]};
    
    hid_t DSET = H5Dcreate2(group_id, dset_name, H5T_NATIVE_DOUBLE, DSPACE, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(DSET, H5T_NATIVE_DOUBLE, DSPACE, DSPACE, H5P_DEFAULT, &sys_param_post);
    H5Dclose(DSET);
    H5Sclose(DSPACE);
    
    double step = cfl_cfg_getE_float(prop_post, "TIMESTEP",-1);
    return step;
}
    

void cfl_hdf5_system_parameter(const char* filename, cfl_init I, cfl_potential* V_pre, cfl_potential* V_post) {
    
    hid_t file_id;
    herr_t status;
    if( access( filename, F_OK ) != -1 )
	file_id = H5Fopen(filename,H5F_ACC_RDWR,H5P_DEFAULT);
    else
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    status = H5Eset_auto(H5E_DEFAULT,NULL, NULL);
    status = H5Gget_objinfo(file_id, "/parameter", 0, NULL);
    if(status==0) {
	fprintf(stderr,"ERROR #1 cfl_hdf5_write_simInfo: Parameter group already exists.\n");
	exit(EXIT_FAILURE);
    }
    else {
	hid_t DSET;
	hid_t group_param = H5Gcreate2(file_id, "/parameter", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
	int Rank_Gdims = 1;
	hsize_t dims_Gdims[Rank_Gdims];
	dims_Gdims[0] = 6;
	hid_t DSPACE_Gdims = H5Screate_simple(Rank_Gdims, dims_Gdims, NULL);
	
	double step = 1.;
	int dimension = 0;
	if(I.xD>1) { dimension += 1; }
	if(I.yD>1) { dimension += 1; }
	if(I.zD>1) { dimension += 1; }
	
//-----physical constants/units
	double length_scale = I.length_scale;
	double temp_scale = I.temperature_scale;
	double hbar_scale = I.hbar_scale;
	double kb_scale = I.k_B_scale*I.hbar_scale;
	double as_scale = I.a_s_scale;
	double mass = I.mass;
	double phys_units[6] = {length_scale,temp_scale,hbar_scale,kb_scale,as_scale,mass};
	
	DSET = H5Dcreate2(group_param, "physical_units", H5T_NATIVE_DOUBLE, DSPACE_Gdims, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(DSET, H5T_NATIVE_DOUBLE, DSPACE_Gdims, DSPACE_Gdims, H5P_DEFAULT, &phys_units);
	H5Dclose(DSET);
	//------------------------------------------------------------------------------------------
	
//-----system parameters (pre-quench, snpcglO1)
	step = cfl_hdf5_system_parameter_snpcglO1(group_param, "system_parameter_pre", I, 0, dimension, hbar_scale, as_scale, mass);
	//------------------------------------------------------------------------------------------
	
//-----system parameters (post-quench, snpcglO1)
	step = cfl_hdf5_system_parameter_snpcglO1(group_param, "system_parameter_post", I, 1, dimension, hbar_scale, as_scale, mass);
	//------------------------------------------------------------------------------------------
	
//-----grid parameter (only for post-quench interesting)
	double dx = (I.xD>1) ? (I.Sgridding*I.dx_scale) : 1.;
	double dy = (I.yD>1) ? (I.Sgridding*I.dy_scale) : 1.;
	double dz = I.Sgridding;
	double d3r = dx*dy*dz;
	double dt_num = I.Tgridding;
	double dt = (I.Tgridding*step*((int)(I.obs_time/(1000.*I.Tgridding*step)))*1000.);
	double grid_dims[6] = {dx,dy,dz,d3r,dt,dt_num};
	
	DSET = H5Dcreate2(group_param, "grid_dimensions", H5T_NATIVE_DOUBLE, DSPACE_Gdims, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(DSET, H5T_NATIVE_DOUBLE, DSPACE_Gdims, DSPACE_Gdims, H5P_DEFAULT, &grid_dims);
	H5Dclose(DSET);
	//------------------------------------------------------------------------------------------
	
//-----Potential parameters (pre-quench)
	hid_t group_V_pre = H5Gcreate2(group_param, "V_pre", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	cfl_potential_write_sysParam(V_pre, group_V_pre);
	cfl_hdf5_write_complex(V_pre->V, filename, "/parameter/V_pre/V", -1);
	H5Gclose (group_V_pre);
	//------------------------------------------------------------------------------------------
	
//-----Potential parameters (post-quench)
	hid_t group_V_post = H5Gcreate2(group_param, "V_post", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	cfl_potential_write_sysParam(V_post, group_V_post);
	cfl_hdf5_write_complex(V_post->V, filename, "/parameter/V_post/V", -1);
	H5Gclose (group_V_post);
	//------------------------------------------------------------------------------------------
	
	H5Gclose (group_param);
	H5Sclose(DSPACE_Gdims);
    }
    
    H5Fclose(file_id);
    
}
