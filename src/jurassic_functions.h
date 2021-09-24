/*
	 This file is part of JURASSIC.

	 JURASSIC is free software: you can redistribute it and/or modify
	 it under the terms of the GNU General Public License as published by
	 the Free Software Foundation, either version 3 of the License, or
	 (at your option) any later version.

	 JURASSIC is distributed in the hope that it will be useful,
	 but WITHOUT ANY WARRANTY; without even the implied warranty of
	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	 GNU General Public License for more details.

	 You should have received a copy of the GNU General Public License
	 along with JURASSIC. If not, see <http://www.gnu.org/licenses/>.

	 Copright (C) 2003-2015 Forschungszentrum Juelich GmbH
	 */

/*!
	\file
	JURASSIC library declarations.

	\mainpage

	The JUelich RApid Spectral SImulation Code (JURASSIC) is a fast radiative
	transfer model for the mid-infrared spectral region.

	This reference manual provides information on the algorithms
	and data structures used in the code. Further information can be found at:
http://www.fz-juelich.de/ias/jsc/jurassic
*/
#pragma once

#include "jurassic.h"

/*! Compose state vector or parameter vector. */
size_t jur_atm2x(
		ctl_t const *ctl,
		atm_t const *atm,
		gsl_vector * x,
		int *iqa,
		int *ipa);

/*! Add elements to state vector. */
void jur_atm2x_help(
		atm_t const *atm,
		double zmin,
		double zmax,
		double const *value,
		int val_iqa,
		gsl_vector * x,
		int *iqa,
		int *ipa,
		size_t * n);

/*! Compute brightness temperature. */
double jur_brightness(
		double const rad,
		double const nu) __deprecated__;

/*! Convert Cartesian coordinates to geolocation. */
// void cart2geo(
// 		double const *x,
// 		double *z,
// 		double *lon,
// 		double *lat) __deprecated__;

/*! Interpolate climatological data. */
void jur_climatology(
		ctl_t const * ctl,
		atm_t * atm_mean);

/*! Compute carbon dioxide continuum (optical depth). */
double jur_ctmco2(
		double nu,
		double p,
		double t,
		double u) __deprecated__;

/*! Compute water vapor continuum (optical depth). */
double jur_ctmh2o(
		double nu,
		double p,
		double t,
		double q,
		double u) __deprecated__;

/*! Compute nitrogen continuum (absorption coefficient). */
double jur_ctmn2(
		double nu,
		double p,
		double t) __deprecated__;

/*! Compute oxygen continuum (absorption coefficient). */
double jur_ctmo2(
		double nu,
		double p,
		double t) __deprecated__;

/*! Copy and initialize atmospheric data. */
void jur_copy_atm(
		ctl_t const *ctl,
		atm_t * atm_dest,
		atm_t const *atm_src,
		int const init);

/*! Copy and initialize observation data. */
void jur_copy_obs(
		ctl_t const *ctl,
		obs_t * obs_dest,
		obs_t const *obs_src,
		int const init);

/*! Find index of an emitter. */
int jur_find_emitter(
		ctl_t const * const ctl,
		char const * const emitter);

/*! Determine ray paths and compute radiative transfer. */
void jur_formod(ctl_t const *ctl, 
    atm_t *atm, 
    obs_t *obs, 
    aero_t const *aero, 
    int n);

/*! Compute radiative transfer for a pencil beam. */
void fr_formod_pencil(
		ctl_t const *ctl,
		atm_t * atm,
		obs_t * obs,
		int const ir);

/*! Apply RFM for radiative transfer calculations. */
void jur_formod_rfm(
		ctl_t const *ctl,
		atm_t * atm,
		obs_t * obs);

/*! Compute Planck source function. */
void jur_formod_srcfunc(
		ctl_t const *ctl,
		trans_table_t const *tbl,
		double const t,
		double *src);

/*! Convert geolocation to Cartesian coordinates. */
// void geo2cart(
// 		double z,
// 		double lon,
// 		double lat,
// 		double *x) __deprecated__;

// /*! Determine gravity of Earth. */
// double gravity(
// 		double z,
// 		double lat) __deprecated__;

/*! Set hydrostatic equilibrium. */
void jur_hydrostatic(
		ctl_t const *ctl,
		atm_t * atm) __deprecated__;

/*! Set hydrostatic equilibrium for individual profile. */
void jur_hydrostatic_1d(
		ctl_t const *ctl,
		atm_t * atm,
		int const ip0,
		int const ip1) __deprecated__;

/*! Determine name of state vector quantity for given index. */
void jur_idx2name(
		ctl_t * ctl,
		int idx,
		char *quantity);

/*! Initialize look-up tables. */
void jur_init_tbl(
		ctl_t const * ctl,
		trans_table_t * tbl);

/*! Interpolate complete atmospheric data set. */
void jur_intpol_atm(
		ctl_t * ctl,
		atm_t * atm_dest,
		atm_t * atm_src);

/*! Interpolate atmospheric data for given geolocation. */
void jur_intpol_atm_geo(
		ctl_t const *ctl,
		atm_t *atm,
		double const z0,
		double const lon0,
		double const lat0,
		double *p,
		double *t,
		double *q,
		double *k);

/*! Interpolate 1D atmospheric data (vertical profile). */
void jur_intpol_atm_1d(
		ctl_t const *ctl,
		atm_t const *atm,
		int const idx0,
		int const n,
		double const z0,
		double *p,
		double *t,
		double *q,
		double *k);

/*! Interpolate 2D atmospheric data (satellite track). */
void jur_intpol_atm_2d(
		ctl_t const *ctl,
		atm_t *atm,
		double const z0,
		double const lon0,
		double const lat0,
		double *p,
		double *t,
		double *q,
		double *k);

/*! Interpolate 3D atmospheric data (Lagrangian grid). */
void jur_intpol_atm_3d(
		ctl_t const *ctl,
		atm_t *atm,
		double const z0,
		double const lon0,
		double const lat0,
		double *p,
		double *t,
		double *q,
		double *k);

/*! Interpolate emissivity from look-up tables. */
double jur_intpol_tbl_eps(
		trans_table_t const *tbl,
		int const ig,
		int const id,
		int const ip,
		int const it,
		double const u);

/*! Interpolate column density from look-up tables. */
double jur_intpol_tbl_u(
		trans_table_t const *tbl,
		int const ig,
		int const id,
		int const ip,
		int const it,
		double const eps);

/*! Convert seconds to date. */
void jur_jsec2time(
		double const jsec,
		int *year,
		int *mon,
		int *day,
		int *hour,
		int *min,
		int *sec,
		double *remain);

/*! Compute Jacobians. */
void jur_kernel(
		ctl_t const *ctl,
		atm_t *atm,
		obs_t *obs,
		gsl_matrix * k) __deprecated__;

/*! Compose measurement vector. */
size_t jur_obs2y(
		ctl_t const *ctl,
		obs_t const *obs,
		gsl_vector * y,
		int *ida,
		int *ira);

/*! Compute Planck function. */
double jur_planck(
		double const t,
		double const nu);

void jur_altitude_range(
		atm_t const *atm,
		double *zmin,
		double *zmax) __deprecated__;

/*! Change segment lengths according to trapezoid rule */
void jur_trapezoid_rule(
		int const np,
		double ds[]) __deprecated__;

/*! Read atmospheric data. */
void jur_read_atm(
		const char *dirname,
		const char *filename,
		ctl_t * ctl,
		atm_t * atm);

/*! Read forward model control parameters. */
void jur_read_ctl(
		int argc,
		char *argv[],
		ctl_t * ctl);

/*! Read matrix. */
void jur_read_matrix(
		const char *dirname,
		const char *filename,
		gsl_matrix * matrix);

/*! Read observation data. */
void jur_read_obs(
		const char *dirname,
		const char *filename,
		ctl_t * ctl,
		obs_t * obs);

/*! Read observation data in RFM format. */
double jur_read_obs_rfm(
		const char *basename,
		double z,
		double *nu,
		double *f,
		int n);

/*! Read RFM spectrum. */
void jur_read_rfm_spec(
		const char *filename,
		double *nu,
		double *rad,
		int *npts);

/*! Read shape function. */
int jur_read_shape(
		const char *filename,
		double *x,
		double *y,
        int const checkmode);

// /*! Compute refractivity (return value is n - 1). */
// double refractivity(
// 		double p,
// 		double t) __deprecated__;

/*! Search control parameter file for variable entry. */
double jur_scan_ctl(
		int argc,
		char *argv[],
		const char *varname,
		int arridx,
		const char *defvalue,
		char *value);

/*! Convert date to seconds. */
void jur_time2jsec(
		int year,
		int mon,
		int day,
		int hour,
		int min,
		int sec,
		double remain,
		double *jsec);

/*! Measure wall-clock time. */
double jur_timer(
		const char *name,
		const char *file,
		const char *func,
		int line,
		int mode);

/*! Write atmospheric data. */
void jur_write_atm(
		const char *dirname,
		const char *filename,
		ctl_t * ctl,
		atm_t * atm);

/*! Write atmospheric data in RFM format. */
void jur_write_atm_rfm(
		const char *filename,
		ctl_t const *ctl,
		atm_t const *atm);

/*! Write matrix. */
void jur_write_matrix(
		const char *dirname,
		const char *filename,
		ctl_t * ctl,
		gsl_matrix * matrix,
		atm_t * atm,
		obs_t * obs,
		const char *rowspace,
		const char *colspace,
		const char *sort);

/*! Write observation data. */
void jur_write_obs(
		const char *dirname,
		const char *filename,
		ctl_t * ctl,
		obs_t * obs);

/*! Decompose parameter vector or state vector. */
void jur_x2atm(
		ctl_t const *ctl,
		gsl_vector *x,
		atm_t * atm);

/*! Extract elements from state vector. */
void jur_x2atm_help(
		atm_t * atm,
		double zmin,
		double zmax,
		double *value,
		gsl_vector * x,
		size_t * n);

/*! Decompose measurement vector. */
void jur_y2obs(ctl_t * ctl, gsl_vector * y, obs_t * obs);

/* Helpers */
FILE* jur_mkFile(const char*, const char*, const char*);
void jur_shell(const char*, const char*);
// double c01(double);

/* stringify the value of a macro, two expansion levels needed */
#define xstr(a) str(a)
#define str(b) #b

// Added:
void jur_table_initialization(ctl_t *ctl);  
