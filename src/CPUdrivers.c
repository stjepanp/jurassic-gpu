#include "jr_common.h" // inline definitions of common functions for CPU and GPU code
#include "scatter_lineofsight.h"

	// ################ CPU driver routines - keep consistent with GPUdrivers.cu ##############

    __host__
	void radiance_to_brightness_CPU(ctl_t const *ctl, obs_t *obs) {
#pragma omp parallel for
		for(int ir = 0; ir < obs->nr; ir++) { // loop over rays
			for(int id = 0; id < ctl->nd; id++) { // loop over detectors
				// convert in-place
				obs->rad[ir][id] = brightness_core(obs->rad[ir][id], ctl->nu[id]);
			} // id
		} // ir
	} // radiance_to_brightness

	__host__
	void surface_terms_CPU(trans_table_t const *tbl, obs_t *obs, double const tsurf[], int const nd) {
#pragma omp  parallel for
		for(int ir = 0; ir < obs->nr; ir++) { // loop over rays
			for(int id = 0; id < nd; id++) { // loop over detectors
				add_surface_core(obs, tbl, tsurf[ir], ir, id);
			} // id
		} // ir
	} // surface_terms_CPU

	/*__host__
	double multi_continua_CPU(char const fourbit, ctl_t const *ctl, 
        pos_t const *los, int const ig_co2, int const ig_h2o, int const id) {
        continua_core_non_temp(fourbit, ctl, los, ig_co2, ig_h2o, id);
	} // multi_continua_CPU */

  __host__ __device__ __ext_inline__
  double continua_core_fourbit(char const fourbit, ctl_t const *ctl, pos_t const *los, 
                               int const ig_co2, int const ig_h2o, int const id) {
      char const CO2 = 8 & fourbit; // 1 << 3 
      char const H2O = 4 & fourbit; // 1 << 2
      char const N2 = 2 & fourbit;  // 1 << 1
      char const O2 = 1 & fourbit;  // 1 << 0
      double const p = los->p;
      double const t = los->t;
      double const ds = los->ds;
      double beta_ds = los->k[ctl->window[id]]*ds;  // extinction
      // make sure that ig_co2 and ig_h2o are both >= 0
      if(CO2) beta_ds += continua_ctmco2(ctl->nu[id], p, t, los->u[ig_co2]);  // co2 continuum
      if(H2O) beta_ds += continua_ctmh2o(ctl->nu[id], p, t, los->q[ig_h2o], los->u[ig_h2o]);  // h2o continuum
      if(N2)  beta_ds += continua_ctmn2(ctl->nu[id], p, t)*ds;  // n2 continuum
      if(O2)  beta_ds += continua_ctmo2(ctl->nu[id], p, t)*ds;  // o2 continuum
      return     beta_ds;
  } // continua_core_fourbit

	__host__
	void apply_kernels_CPU(trans_table_t const *tbl, ctl_t const *ctl, obs_t *obs, 
        pos_t const (*restrict los)[NLOS], int const np[], 
        int const ig_co2, int const ig_h2o, char const fourbit,
        double const (*restrict aero_beta)[ND]) { // aero_beta is added

#pragma omp parallel for
		for(int ir = 0; ir < obs->nr; ir++) { // loop over independent rays
			double tau_path[ND][NG]; // private for each ray
			for(int id = 0; id < ND; id++) { // loop over detectors
				obs->rad[ir][id] = 0.0;
				obs->tau[ir][id] = 1.0;
				for(int ig = 0; ig < NG; ig++) { // loop over gases
					tau_path[id][ig] = 1.0;
				} // ig
			} //  id

			for(int ip = 0; ip < np[ir]; ++ip) { // loop over line-of-sight points
				for(int id = 0; id < ctl->nd; id++) { // loop over detector channels

					// compute extinction coefficient
					double const beta_ds = continua_core_fourbit(fourbit, ctl, &(los[ir][ip]), ig_co2, ig_h2o, id);
					
          //Added:
          double const aero_ds = los[ir][ip].aerofac * aero_beta[los[ir][ip].aeroi][id] * los[ir][ip].ds;
          
          // compute transmission with the EGA method
					double const tau_gas = apply_ega_core(tbl, &(los[ir][ip]), tau_path[id], ctl->ng, id);
					// compute the source term
					double const planck = src_planck_core(tbl, los[ir][ip].t, id);
					// perform integration
					new_obs_core(obs, ir, id, beta_ds + aero_ds, planck, tau_gas);

				} // id --> could be vectorized over detector channels
#ifdef GPUDEBUG
				assert(los[ir][ip].ip == ip); // sanity check
				assert(los[ir][ip].ir == ir); // sanity check
#endif
			} // ip --> non-parallelisable due to loup carried dependency
		} // ir --> OpenMP parallel over rays

	} // apply_kernels_CPU

	__host__
	void raytrace_rays_CPU(ctl_t const *ctl, atm_t const *atm, obs_t *obs, 
                           pos_t los[NR][NLOS], double tsurf[], int np[]) {
#pragma omp parallel for
		for(int ir = 0; ir < obs->nr; ir++) { // loop over rays
			np[ir] = traceray(ctl, atm, obs, ir, los[ir], &tsurf[ir]);
		} // ir
	} // raytrace_rays_CPU

	__host__
	void hydrostatic1d_CPU(ctl_t const *ctl, atm_t *atm, int const nr, int const ig_h2o) {
			if(ctl->hydz < 0) return; // Check reference height
			for(int ir = 0; ir < nr; ir++) { // Apply hydrostatic equation to individual profiles
				hydrostatic_1d_h2o(ctl, atm, 0, atm->np, ig_h2o);
			} // ir
	} // hydrostatic1d_CPU
  
  __host__
  void get_los_and_convert_los_t_to_pos_t_CPU(pos_t (*pos)[NLOS], int *np, double *tsurf, ctl_t const *ctl, 
                                              atm_t const *atm, obs_t *obs, aero_t *aero) {
#pragma omp parallel for
    for(int ir = 0; ir < obs->nr; ir++) {
      los_t *one_los;
      one_los = (los_t*) malloc(sizeof(los_t));
      // raytracer copied from jurassic-scatter
      jur_raytrace(ctl, atm, obs, aero, one_los, ir);
      np[ir] = one_los->np;
      tsurf[ir] = one_los->tsurf;
      for(int ip = 0; ip < one_los->np; ip++) { 
        convert_los_to_pos_core(&pos[ir][ip], one_los, ip);
      }
      free(one_los);
    }
  }

	// ################ end of CPU driver routines ##############

	// The full forward model on the CPU ////////////////////////////////////////////
	__host__
	void formod_CPU(ctl_t const *ctl, atm_t *atm, obs_t *obs,
                  aero_t const *aero) {
    printf("DEBUG #%d formod_CPU was called!\n", ctl->MPIglobrank);
  
    if (ctl->checkmode) {
      printf("# %s: checkmode = %d, no actual computation is performed!\n", __func__, ctl->checkmode);
      return; // do nothing here
    } // checkmode

		assert(obs);

		char mask[NR][ND];
		save_mask(mask, obs, ctl);

		trans_table_t const *tbl = get_tbl(ctl);
		double *t_surf = (double*)malloc((obs->nr)*sizeof(double));
		int *np = (int*)malloc((obs->nr)*sizeof(int));
		pos_t (*los)[NLOS] = (pos_t (*)[NLOS])malloc((obs->nr)*(NLOS)*sizeof(pos_t));

		// gas absorption continua configuration
		static int ig_co2 = -999, ig_h2o = -999;
		if((ctl->ctm_h2o) && (-999 == ig_h2o)) ig_h2o = jur_find_emitter(ctl, "H2O");
		if((ctl->ctm_co2) && (-999 == ig_co2)) ig_co2 = jur_find_emitter(ctl, "CO2");
		// binary switches for the four gases
		char const fourbit = (char) (
                  ( (1 == ctl->ctm_co2) && (ig_co2 >= 0) )*8   // CO2
				+ ( (1 == ctl->ctm_h2o) && (ig_h2o >= 0) )*4   // H2O
				+   (1 == ctl->ctm_n2)                    *2   // N2
				+   (1 == ctl->ctm_o2)                    *1); // O2

 
    hydrostatic1d_CPU(ctl, atm, obs->nr, ig_h2o); // in this call atm might get modified
    // if formod function was NOT called from jurassic-scatter project
    if(ctl->leaf_nr == -1) {
      raytrace_rays_CPU(ctl, atm, obs, los, t_surf, np);
    } else {  
      get_los_and_convert_los_t_to_pos_t_CPU(los, np, t_surf, ctl, atm, obs, aero);
    }

    // "beta_a" -> 'a', "beta_e" -> 'e'
    char const beta_type = ctl->sca_ext[5];

    apply_kernels_CPU(tbl, ctl, obs, los, np, ig_co2, ig_h2o, fourbit,
                      (('a' == beta_type) ? aero->beta_a : aero->beta_e));
    surface_terms_CPU(tbl, obs, t_surf, ctl->nd);

		free(los);
		free(np);
		free(t_surf);

    if(ctl->write_bbt && ctl->leaf_nr == -1)  // convert radiance to brightness (in-place)
		  radiance_to_brightness_CPU(ctl, obs);

		apply_mask(mask, obs, ctl);
	} // formod_CPU

	__host__ void formod_GPU(ctl_t const *ctl, atm_t *atm, obs_t *obs,
                            aero_t const *aero, int n)
#ifdef hasGPU
    ; // declaration only, will be provided by GPUdrivers.o at link time 
#else
    { // definition here
        static int warnGPU = 1;
        if (ctl->useGPU > 0) { // USEGPU > 0 means use-GPU-always
            fprintf(stdout, "CUDA not found during compilation, cannot run on GPUs. ABORT\n\n");
            fprintf(stdout, "USEGPU = 1 (use-GPU-always) found in controls\n"
                            "USEGPU = -1 (use-GPU-if-possible) could help\n\n");
            fflush(stdout);
            fprintf(stderr, "CUDA not found during compilation, cannot run on GPUs. ABORT\n\n");
            exit(EXIT_FAILURE);
        } else {               // USEGPU < 0 means use-GPU-if-possible
            assert(ctl->useGPU < 0); 
            // automatic decision: fallback solution is CPU
            if (warnGPU) { // this warning appears only once per process
                printf("CUDA not found during compilation, continue on CPUs instead!\n");
                warnGPU = 0; // switch this warning off
            } // warnGPU
            //TODO: Add parallelism!
            for(int i = 0; i < n; i++) {
              //aero and los are optionl
              formod_CPU(ctl, atm, &obs[i], aero);
            }
        } //
    } // formod_GPU
#endif

	__host__
	void jur_formod(ctl_t const *ctl, atm_t *atm, obs_t *obs, aero_t const *aero, int n) {
        if (ctl->checkmode) {
            static int nr_last_time = -999;
            if (obs->nr != nr_last_time) {
                printf("# %s: %d max %d rays , %d of max %d gases, %d of max %d channels\n", 
                        __func__, obs->nr, NR, ctl->ng, NG, ctl->nd, ND);
                // the number of rays can be zero if we skipped to read obs.tab, checkmode=1
                nr_last_time = obs->nr;
            } // only report if nr changed
        } // checkmode
        if (ctl->useGPU) { //has to be changed
            formod_GPU(ctl, atm, obs, aero, n);
        } else { // USEGPU = 0 means use-GPU-never
            for(int i = 0; i < n; i++) {
              //aero and los are optionl
              formod_CPU(ctl, atm, &obs[i], aero);
            }
        } // useGPU
    } // formod

//we could use the same trick as above but it's not necessary
__host__
trans_table_t* get_tbl_on_GPU(ctl_t const *ctl)
;

__host__
void jur_table_initialization(ctl_t *ctl) {
  double tic = omp_get_wtime(); 
  trans_table_t *tbl;
  if(ctl->useGPU) {
    printf("DEBUG #%d call initilaze GPU..\n", ctl->MPIglobrank);
    tbl = get_tbl_on_GPU(ctl); 
  }
  else {
    printf("DEBUG #%d call initilaze CPU..\n", ctl->MPIglobrank);
    tbl = get_tbl(ctl);
  }
  double toc = omp_get_wtime();
  printf("TIMER #%d jurassic-gpu table initialization time: %lf\n", ctl->MPIglobrank, toc - tic);
}
