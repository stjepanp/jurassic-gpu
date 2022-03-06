  template<int CO2, int H2O, int N2, int O2>
  __host__ __device__ __ext_inline__
  double continua_core_temp(ctl_t const *ctl, pos_t const *los, int const ig_co2, int const ig_h2o, int const id) {
      double const p = los->p;
      double const t = los->t;
      double const ds = los->ds;
      double beta_ds = los->k[ctl->window[id]]*ds;													// extinction
      // make sure that ig_co2 and ig_h2o are both >= 0
      if ((CO2)) beta_ds += continua_ctmco2(ctl->nu[id], p, t, los->u[ig_co2]);						// co2 continuum
      if ((H2O)) beta_ds += continua_ctmh2o(ctl->nu[id], p, t, los->q[ig_h2o], los->u[ig_h2o]);		// h2o continuum
      if ((N2))  beta_ds += continua_ctmn2(ctl->nu[id], p, t)*ds;									// n2 continuum
      if ((O2))  beta_ds += continua_ctmo2(ctl->nu[id], p, t)*ds;									// o2 continuum
      return     beta_ds;
  } // continua_core_bbbb where each b is either 0 or 1
  
  template<int CO2, int H2O, int N2, int O2>
    void __global__ // GPU-kernel
    fusion_kernel_GPU(trans_table_t const *tbl, ctl_t const *ctl,
        obs_t *obs, pos_t const (*restrict los)[NLOS],
        int const np[], int const ig_co2, int const ig_h2o,
        double const (*restrict aero_beta)[ND]) {
			double tau_path[NG];
			for(int ir = blockIdx.x; ir < obs->nr; ir += gridDim.x) { // grid stride loop over blocks = rays
				for(int id = threadIdx.x; id < ND; id += blockDim.x) { // grid stride loop over threads = detectors
					obs->rad[ir][id] = 0.0;
					obs->tau[ir][id] = 1.0;
				} // id
				for(int ig = 0; ig < NG; ++ig) {
					tau_path[ig] = 1.0;
				} // ig
				for(int ip = 0; ip < np[ir]; ++ip) {
					for(int id = threadIdx.x; id < ctl->nd; id += blockDim.x) { // grid stride loop over threads = detectors
						double const beta_ds = continua_core_temp <CO2, H2O, N2, O2>
                                   (ctl, &(los[ir][ip]), ig_co2, ig_h2o, id);           // function args
						
            //Added:
            double const aero_ds = los[ir][ip].aerofac * aero_beta[los[ir][ip].aeroi][id] * los[ir][ip].ds;
            
            double const tau_gas = apply_ega_core(tbl, &(los[ir][ip]), tau_path, ctl->ng, id);
						double const planck = src_planck_core(tbl, los[ir][ip].t, id);
						new_obs_core(obs, ir, id, beta_ds + aero_ds, planck, tau_gas);
					} // id --> parallel over detectors=threads
				} // ip --> non-parallelisable
			} // ir --> parallel over rays==blocks
    } // fusion_kernel_GPU
