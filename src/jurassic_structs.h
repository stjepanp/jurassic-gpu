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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>

/* ------------------------------------------------------------
	 Macros...
	 ------------------------------------------------------------ */
/*! Allocate memory. */
#define ALLOC(ptr, type, n)				\
	if((ptr=(type*)malloc((size_t)(n)*sizeof(type)))==NULL)	\
ERRMSG("Out of memory!");

/*! Compute Cartesian distance between two vectors. */
#define DIST(a, b) sqrt(DIST2(a, b))

/*! Compute squared distance between two vectors. */
#define DIST2(a, b)							\
	((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))

/*! Compute dot product of two vectors. */
#define DOTP(a, b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

/*! Print error message and quit program. */
#define ERRMSG(msg) {							\
	printf("\nError (%s, %s, l%d): %s\n\n",				\
			__FILE__, __func__, __LINE__, msg);				\
	exit(EXIT_FAILURE);							\
}

/*! Compute exponential interpolation. */
#define EXP(x0, y0, x1, y1, x)					\
	(((y0)>0 && (y1)>0)						\
	 ? ((y0)*exp(log((y1)/(y0))/((x1)-(x0))*((x)-(x0))))					\
	 : LIN(x0, y0, x1, y1, x))

/*! Compute linear interpolation. */
#define LIN(x0, y0, x1, y1, x) ((y0) + ((x) - (x0))*((y1) - (y0))/((x1) - (x0)))

/*! Compute norm of a vector. */
#define NORM(a) sqrt(DOTP(a, a))

/*! Print macro for debugging. */
#define PRINT(format, var)						\
	printf("Print (%s, %s, l%d): %s= " format "\n",				\
			__FILE__, __func__, __LINE__, #var, var);

/*! Start or stop a timer. */
#define TIMER(name, mode) jur_timer(name, __FILE__, __func__, __LINE__, mode)

/*! Read string tokens. */
#define TOK(line, tok, format, var, saveptr) {			\
	if(((tok)=strtok_r((line), " \t",saveptr))) {			\
		if(sscanf(tok, format, &(var))!=1) continue;	\
	} else ERRMSG("Error while reading!");		\
}

#define __deprecated__ __attribute__((deprecated))

/* stringify the value of a macro, two expansion levels needed */
#define xstr(a) str(a)
#define str(b) #b

// double eip(const double, const double, const double, const double, const double) __deprecated__;
// double lip(const double, const double, const double, const double, const double) __deprecated__;

/* ------------------------------------------------------------
	 Constants...
	 ------------------------------------------------------------ */

/*! First spectroscopic constant (c_1 = 2 h c^2) [W/(m^2 sr cm^-4)]. */
#define C1 1.19104259e-8

/*! Second spectroscopic constant (c_2 = h c / k) [K/cm^-1]. */
#define C2 1.43877506

/*! Standard gravity [m/s^2]. */
#define G0 9.80665

/*! Standard pressure [hPa]. */
#define P0 1013.25

/*! Standard temperature [K]. */
#define T0 273.15

/*! Mean radius of Earth [km]. */
#define RE 6367.421

/*! Mass of Earth [kg]. */
#define ME 5.976e24



/* ------------------------------------------------------------
	 Dimensions...
	 ------------------------------------------------------------ */

/*! Maximum number of radiance channels. */
#ifndef ND
  #define ND 100 
#endif

/*! Maximum number of emitters. */
#ifndef NG
  #define NG 30
#endif  

/*! Maximum number of atmospheric data points. */
#define NP 9600

/*! Maximum number of ray paths. */
#define NR 1088

/*! Maximum number of spectral windows. */
#define NW 1

/*! Maximum length of ASCII data lines. */
#define LEN 5000

/*! Maximum size of measurement vector. */
#define M (NR*ND)

/*! Maximum size of state vector. */
#define N (NQ*NP)

/*! Maximum number of quantities. */
#define NQ (2+NG+NW)

/*! Maximum number of LOS points. */
#define NLOS 400

/*! Maximum number of shape function grid points. */
#define NSHAPE 2048

/*! Number of ray paths used for FOV calculations. */
#define NFOV 5

/*! Maximum number of pressure levels in emissivity tables. */
#define TBLNP 40

/*! Maximum number of temperatures in emissivity tables. */
#define TBLNT 30

/*! Maximum number of column densities in emissivity tables. */
#define TBLNU 304

/*! Maximum number of source function temperature levels. */
#define TBLNS 1201

/*! Maximum number of RFM spectral grid points. */
#define RFMNPTS 10000000

/*! Maximum length of RFM data lines. */
#define RFMLINE 100000

/* ------------------------------------------------------------
	 Quantity indices...
	 ------------------------------------------------------------ */

/*! Index for pressure. */
#define IDXP 0

/*! Index for temperature. */
#define IDXT 1

/*! Indices for volume mixing ratios. */
#define IDXQ(ig) (2+ig)

/*! Indices for extinction. */
#define IDXK(iw) (2+ctl->ng+iw)

/* ------------------------------------------------------------
	 Structs...
	 ------------------------------------------------------------ */

typedef struct { /// Atmospheric data. //////////////////////////////////////////
	double time[NP];			 /// Time (seconds since 2000-01-01T00:00Z).
	double z[NP];					 /// Altitude [km].
	double lon[NP];				 /// Longitude [deg].
	double lat[NP];				 /// Latitude [deg].
	double p[NP];					 /// Pressure [hPa].
	double t[NP];					 /// Temperature [K].
	double q[NG][NP];			 /// Volume mixing ratio.
	double k[NW][NP];			 /// Extinction [1/km].
	int np;                /// Number of data points.
	int init;							 /// Init flag for interpolation (0=no, 1=yes).
} atm_t; // ////////////////////////////////////////////////////////////////////

/*! Forward model control parameters. */
typedef struct {

	/*! Number of emitters. */
	int ng;

	/*! Name of each emitter. */
	char emitter[NG][LEN];

	/*! Number of radiance channels. */
	int nd;

	/*! Number of spectral windows. */
	int nw;

	/*! Centroid wavenumber of each channel [cm^-1]. */
	double nu[ND];

	/*! Window index of each channel. */
	int window[ND];

	/*! Basename for table files and filter function files. */
	char tblbase[LEN];

	/*! Reference height for hydrostatic pressure profile (-999 to skip) [km]. */
	double hydz;

	/*! Compute CO2 continuum (0=no, 1=yes). */
	int ctm_co2;

	/*! Compute H2O continuum (0=no, 1=yes). */
	int ctm_h2o;

	/*! Compute N2 continuum (0=no, 1=yes). */
	int ctm_n2;

	/*! Compute O2 continuum (0=no, 1=yes). */
	int ctm_o2;

	/*! Interpolation method (1=profile, 2=satellite track, 3=Lagrangian grid). */
	int ip;

	/*! Influence length for vertical interpolation [km]. */
	double cz;

	/*! Influence length for horizontal interpolation [km]. */
	double cx;

	/*! Take into account refractivity (0=no, 1=yes). */
	int refrac;

	/*! Maximum step length for raytracing [km]. */
	double rayds;

	/*! Vertical step length for raytracing [km]. */
	double raydz;

	/*! Field-of-view data file. */
	char fov[LEN];

	/*! Minimum altitude for pressure retrieval [km]. */
	double retp_zmin;

	/*! Maximum altitude for pressure retrieval [km]. */
	double retp_zmax;

	/*! Minimum altitude for temperature retrieval [km]. */
	double rett_zmin;

	/*! Maximum altitude for temperature retrieval [km]. */
	double rett_zmax;

	/*! Minimum altitude for volume mixing ratio retrieval [km]. */
	double retq_zmin[NG];

	/*! Maximum altitude for volume mixing ratio retrieval [km]. */
	double retq_zmax[NG];

	/*! Minimum altitude for extinction retrieval [km]. */
	double retk_zmin[NW];

	/*! Maximum altitude for extinction retrieval [km]. */
	double retk_zmax[NW];

	/*! Use brightness temperature instead of radiance (0=no, 1=yes). */
	int write_bbt;

	/*! Write matrix file (0=no, 1=yes). */
	int write_matrix;

	/*! Forward model (1=CGA, 2=EGA, 3=RFM). */
	int formod;

	/*! Path to RFM binary. */
	char rfmbin[LEN];

	/*! HITRAN file for RFM. */
	char rfmhit[LEN];

	/*! Emitter cross-section files for RFM. */
	char rfmxsc[NG][LEN];

	/*! Use GPU-accelerated formod implementation (0=no, 1=yes) */
	int useGPU;
    
    /*! do not perform input, computation, nor output, just make sure files are there */
	int checkmode;

	/*! MPI rank information */
	int MPIglobrank;  // global rank
	int MPIlocalrank; // node-local Rank

	/* binary IO */
    int read_binary;
    int write_binary;

    /*! Shared memory controler for GPU kernels */
    int gpu_nbytes_shared_memory;

} ctl_t;


/*! Point on the Line-of-sight data without storing */
typedef struct {
	double z;		/*! Altitude [km]. */
	double lon;	/*! Longitude [deg]. */
	double lat;	/*! Latitude [deg]. */
	double p;		/*! Pressure [hPa]. */
	double t;		/*! Temperature [K]. */
	double q[NG];	/*! Volume mixing ratio. */
	double k[NW];	/*! Extinction [1/km]. */
	double ds;	/*! Segment length [km]. */
	double u[NG];	/*! Column density [molecules/cm^2]. */
#ifdef CURTIS_GODSON
	double cgp[NG];	/*! Curtis-Godson pressure [hPa]. */
	double cgt[NG];	/*! Curtis-Godson temperature [K]. */
	double cgu[NG];	/*! Curtis-Godson column density [molecules/cm^2]. */
#endif
#ifdef GPUDEBUG
	int ip, ir;  // debug helpers
#endif
} pos_t;

typedef struct { /// Observation geometry and radiance data. ///////////////////
	double time[NR];		/// Time (seconds since 2000-01-01T00:00Z). 
	double obsz[NR];		/// Observer altitude [km]. 
	double obslon[NR];	/// Observer longitude [deg]. 
	double obslat[NR];	/// Observer latitude [deg]. 
	double vpz[NR];			/// View point altitude [km]. 
	double vplon[NR];		/// View point longitude [deg]. 
	double vplat[NR];		/// View point latitude [deg]. 
	double tpz[NR];			/// Tangent point altitude [km]. 
	double tplon[NR];		/// Tangent point longitude [deg]. 
	double tplat[NR];		/// Tangent point latitude [deg]. 
	double tau[NR][ND]; /// Transmittance of ray path.		// transposed
	double rad[NR][ND]; /// Radiance [W/(m^2 sr cm^-1)].	// transposed
	int nr;							/// Number of ray paths.
} obs_t; // ////////////////////////////////////////////////////////////////////

typedef float real_tblND_t;

/*! Emissivity look-up tables. */
typedef struct {

	/*! Number of pressure levels. */
	int32_t np[NG][ND];

	/*! Number of temperatures. */
	int32_t nt[NG][TBLNP][ND];

	/*! Number of column densities. */
	int32_t nu[NG][TBLNP][TBLNT][ND];

	/*! Pressure [hPa]. */
	double p[NG][TBLNP][ND];

	/*! Temperature [K]. */
	double t[NG][TBLNP][TBLNT][ND];

	/*! Column density [molecules/cm^2]. */
	real_tblND_t u[NG][TBLNP][TBLNT][TBLNU][ND];
    
	/*! Emissivity. */
	real_tblND_t eps[NG][TBLNP][TBLNT][TBLNU][ND];

	/*! Source function radiance [W/(m^2 sr cm^-1)]. */
	double sr[TBLNS][ND];

	/*! Source function temperature [K]. */
	double st[TBLNS];

#ifdef  FAST_INVERSE_OF_U
	/*! u0inv[g][p][t][d] * u[g][p][t][0][d] == 1 must hold! */ // FAST_INVERSE_OF_U
	double u0inv[NG][TBLNP][TBLNT][ND];                         // FAST_INVERSE_OF_U
    /*! We assume a logarithmic increment by 2^(1/6) */         // FAST_INVERSE_OF_U
#endif

} tbl_t;
