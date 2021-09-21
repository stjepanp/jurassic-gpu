#ifndef JR_SCATTER_GPU_H
#define JR_SCATTER_GPU_H

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_statistics.h>

#ifdef MPI
#include <mpi.h>
#endif


/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/* Allocate memory. */
#define ALLOC(ptr, type, n)				\
  if((ptr=malloc((size_t)(n)*sizeof(type)))==NULL)	\
    ERRMSG("Out of memory!");

/* Compute angle between two vectors. */
#define ANGLE(a, b)						\
  acos(GSL_MIN(GSL_MAX(DOTP(a, b)/NORM(a)/NORM(b), -1), 1))

/* Compute Cartesian distance between two vectors. */
#define DIST(a, b) sqrt(DIST2(a, b))

/* Compute squared distance between two vectors. */
#define DIST2(a, b)							\
  ((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))

/* Compute dot product of two vectors. */
#define DOTP(a, b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

/* Print error message and quit program. */
#define ERRMSG(msg) {                                                   \
    printf("\nError (%s, %s, l%d): %s\n\n",                             \
           __FILE__, __func__, __LINE__, msg);                      \
    exit(EXIT_FAILURE);                                                 \
  }

/* Compute exponential interpolation. */
#define EXP(x0, y0, x1, y1, x)					\
  (((y0)>0 && (y1)>0)						\
   ? ((y0)*exp(log((y1)/(y0))/((x1)-(x0))*((x)-(x0))))          \
   : LIN(x0, y0, x1, y1, x))

/* Read binary data. */
#define FREAD(ptr, type, nmemb, stream) {                               \
    if(fread(ptr, sizeof(type), (size_t)nmemb, stream)!=(size_t)nmemb)  \
      ERRMSG("Error while reading!");                                   \
  }

/* Write binary data. */
#define FWRITE(ptr, type, nmemb, stream) {				\
    if(fwrite(ptr, sizeof(type), (size_t)nmemb, stream)!=(size_t)nmemb)	\
      ERRMSG("Error while writing!");					\
  }

/* Compute linear interpolation. */
#define LIN(x0, y0, x1, y1, x)			\
  ((y0)+((y1)-(y0))/((x1)-(x0))*((x)-(x0)))

/* Write log message. */
#define LOGMSG(lev, cmd) {if(lev<=VERBLEV) cmd;}

/* Execute netCDF library command and check result. */
#define NC(cmd) {				     \
    if((cmd)!=NC_NOERR)				     \
      ERRMSG(nc_strerror(cmd));			     \
  }

/* Compute norm of a vector. */
#define NORM(a) sqrt(DOTP(a, a))

/* Print macro for debugging. */
#define PRINT(format, var)                                              \
  printf("Print (%s, %s, l%d): %s= "format"\n",                         \
         __FILE__, __func__, __LINE__, #var, var);

//Added:

/*! Start or stop a timer. */
#define TIMER(name, mode) jur_timer(name, __FILE__, __func__, __LINE__, mode)

#define __deprecated__ __attribute__((deprecated))

/* stringify the value of a macro, two expansion levels needed */
#define xstr(a) str(a)
#define str(b) #b

/*! Read string tokens. */
// in scatter version there were 4 parameters
#define TOK_FIVE_ARGS(line, tok, format, var, saveptr) {			\
	if(((tok)=strtok_r((line), " \t",saveptr))) {			\
		if(sscanf(tok, format, &(var))!=1) continue;	\
	} else ERRMSG("Error while reading!");		\
}

/* Read string tokens. */
#define TOK(line, tok, format, var) {			\
    if(((tok)=strtok((line), " \t"))) {			\
      if(sscanf(tok, format, &(var))!=1) continue;	\
    } else ERRMSG("Error while reading!");		\
  }

/* ------------------------------------------------------------
   Redefinable constants...
   ------------------------------------------------------------ */

/* Verbosity level. */
#ifndef VERBLEV
#define VERBLEV 2
#endif

/* ------------------------------------------------------------
   Constants...
   ------------------------------------------------------------ */

/* First spectroscopic constant (c_1 = 2 h c^2) [W/(m^2 sr cm^-4)]. */
#define C1 1.19104259e-8

/* Second spectroscopic constant (c_2 = h c / k) [K/cm^-1]. */
#define C2 1.43877506

/* Standard pressure [hPa]. */
#define P0 1013.25

/* Mean radius of Earth [km]. */
#define RE 6367.421

/* Mass of Earth [kg]. */
#define ME 5.976e24

/* Temperature of the Sun [K]. */
#define TSUN 5780.

//Added:
/* Standard gravity [m/s^2]. */
#define G0 9.80665

/* Standard temperature [K]. */
#define T0 273.15


/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum length of ASCII data lines. */
#define LEN 5000 //Ok!

/* Maximum size of measurement vector. */
#define MMAX (NRMAX*NDMAX) //TODO: M (NR*ND)
#define M MMAX

/* Maximum size of state vector. */
#define NMAX (NQMAX*NPMAX) //TODO: N (NQ*NP)
#define N NMAX

/* Maximum number of quantities. */
#define NQMAX (2+NGMAX+NWMAX) //TODO: NQ (2+NG+NW)
#define NQ NQMAX

/* Maximum number of spectral windows. */
#define NWMAX 5 //TODO: NW 1 
#define NW NWMAX

/* Maximum number of radiance channels. */
#define NDMAX 130 //TODO: ifndef ND 100
#ifndef ND
  #define ND NDMAX 
#endif

/* Maximum number of emitters. */
#define NGMAX 25 //TODO: ifndef NG 30
#ifndef NG
  #define NG NGMAX
#endif  

/* Maximum number of LOS points. */
#define NLOS 10000 //TODO: NLOS 400

/* Maximum number of atmospheric data points. */
#define NPMAX 1000 //TODO: NP 9600
#define NP NPMAX

/* Maximum number of ray paths. */
#define NRMAX 1000 //TODO: NR 1088
#define NR NRMAX

/* Maximum number of shape function grid points. */
#define NSHAPE 10000 //TODO: NSHAPE 2048 | TODO: ?

/* Number of ray paths used for FOV calculations. */
#define NFOV 50 //TODO: NFOV 5

/* Maximum number of pressure levels in emissivity tables. */
#define TBLNPMAX 50 //TODO: TBLNP 40
#define TBLNP TBLNPMAX

/* Maximum number of source function temperature levels. */
#define TBLNSMAX 1201 //TODO: TBLNS 1201
#define TBLNS TBLNSMAX

/* Maximum number of temperatures in emissivity tables. */
#define TBLNTMAX 30 //TODO: TBLNT 30
#define TBLNT TBLNTMAX

/* Maximum number of column densities in emissivity tables. */
#define TBLNUMAX 430 //TODO: TBLNU 304
#define TBLNU TBLNUMAX

/* Maximum number of scattering models. */
#define SCAMOD 30

/* Maximum number of aerosol/cloud layers */
#define NLMAX 10

/* Number of scattering angles (from 0 to 180 deg). */
#define NTHETA 181

/* Number of points for Gauss-Hermite integration. */
#define NRAD 170

/* Maximum number of refractive indices. */
#define REFMAX 5000

//Added:

/* Maximum number of RFM spectral grid points. */
#define RFMNPTS 10000000

/* Maximum length of RFM data lines. */
#define RFMLINE 100000

/* ------------------------------------------------------------
   Quantity indices...
   ------------------------------------------------------------ */

/* Pressure. */
#define IDXP 0

/* Temperature. */
#define IDXT 1

/* Volume mixing ratios. */
#define IDXQ(ig) (2+ig)

/* Extinction. */
#define IDXK(iw) (2+ctl->ng+iw)

/* Particle concentration. */
#define IDXNN  (2+ctl->ng+ctl->nw+1)

/* Particle size (radius). */
#define IDXRR  (2+ctl->ng+ctl->nw+2)

/* Particle size distribution width. */
#define IDXSS  (2+ctl->ng+ctl->nw+3)


/* ------------------------------------------------------------
	 Structs...
	 ------------------------------------------------------------ */
typedef struct {
    void* input;
    void* result;
    int ir;
} queue_item_t;

typedef struct {
    queue_item_t* items;
    int capacity;
    int begin;
    int end;
} queue_t;

typedef struct { /// Atmospheric data. /////////////////////////////////////////
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
} atm_t; ///////////////////////////////////////////////////////////////////////

typedef struct { /// Forward model control parameters. /////////////////////////
	int ng;                 /// Number of emitters.
	char emitter[NG][LEN];  /// Name of each emitter.
  int nd;                 /// Number of radiance channels.
	int nw;                 /// Number of spectral windows.
	double nu[ND];          /// Centroid wavenumber of each channel [cm^-1].
	int window[ND];         /// Window index of each channel. 
	char tblbase[LEN];      /// Basename for table files and filter function files.
	double hydz;            /// Reference height for hydrostatic pressure profile (-999 to skip) [km].
	int ctm_co2;            /// Compute CO2 continuum (0=no, 1=yes).
	int ctm_h2o;            /// Compute H2O continuum (0=no, 1=yes).
	int ctm_n2;             /// Compute N2 continuum (0=no, 1=yes).
	int ctm_o2;             /// Compute O2 continuum (0=no, 1=yes).
	int ip;                 /// Interpolation method (1=profile, 2=satellite track, 3=Lagrangian grid).
	double cz;              /// Influence length for vertical interpolation [km].
	double cx;              /// Influence length for horizontal interpolation [km].
	int refrac;             /// Take into account refractivity (0=no, 1=yes).
	double rayds;           /// Maximum step length for raytracing [km].
	double raydz;           /// Vertical step length for raytracing [km].
	char fov[LEN];          /// Field-of-view data file.
	double retp_zmin;       /// Minimum altitude for pressure retrieval [km].
	double retp_zmax;       /// Maximum altitude for pressure retrieval [km].
	double rett_zmin;       /// Minimum altitude for temperature retrieval [km].
	double rett_zmax;       /// Maximum altitude for temperature retrieval [km].
	double retq_zmin[NG];   /// Minimum altitude for volume mixing ratio retrieval [km].
	double retq_zmax[NG];   /// Maximum altitude for volume mixing ratio retrieval [km].
	double retk_zmin[NW];   /// Minimum altitude for extinction retrieval [km].
	double retk_zmax[NW];   /// Maximum altitude for extinction retrieval [km].
	int write_bbt;          /// Use brightness temperature instead of radiance (0=no, 1=yes).
	int write_matrix;       /// Write matrix file (0=no, 1=yes).
	int formod;             /// Forward model (1=CGA, 2=EGA, 3=RFM).
	char rfmbin[LEN];       /// Path to RFM binary.
	char rfmhit[LEN];       /// HITRAN file for RFM.
	char rfmxsc[NG][LEN];   /// Emitter cross-section files for RFM.
	int useGPU;             /// Use GPU-accelerated formod implementation (0=no, 1=yes)
	int checkmode;          /// do not perform input, computation, nor output, just make sure files are there 
	int MPIglobrank;        /// MPI global rank
	int MPIlocalrank;       /// MPI node-local Rank
  int read_binary;        /// binary IO
  int write_binary;       /// binary IO
  int gpu_nbytes_shared_memory; /// Shared memory controler for GPU kernels
  
  /// ---------------- for scattering ------------------
  int sca_n;              /// Number of scattering models.
  int sca_mult;           /// Number of recursions for multiple scattering.
                          /// (0=no scattering, 1=single scattering, 2<=multiple scattering)
  char sca_ext[LEN];      /// Extinction coefficient type if sca_mult=0
  double transs;          /// Sampling step for transition layers [km].
  double retnn_zmin;      /// Minimum altitude for particle [km].
  double retnn_zmax;      /// Maximum altitude for particle retrieval [km].
  int retnn;              /// Retrieval of particle concentration (0=no, 1=yes)
  int retrr;              /// Retrieval of particle size (0=no, 1=yes)
  int retss;              /// Retrieval of particle size distribution width (0=no, 1=yes)
  int leaf_nr;            /// Number of leaf rays, for example number of secondary rays if sca_mult=1
  int queue_state;        /// We have multiple queues, but all of them are in same state
} ctl_t; ///////////////////////////////////////////////////////////////////////

typedef struct {    /// Point on the Line-of-sight data without storing //////////
	double z;		      /// Altitude [km].
	double lon;	      /// Longitude [deg].
	double lat;	      /// Latitude [deg].
	double p;		      /// Pressure [hPa].
	double t;		      /// Temperature [K].
	double q[NG];	    /// Volume mixing ratio.
	double k[NW];	    /// Extinction [1/km].
  int aeroi;        /// Aerosol/cloud layer index
  double aerofac;   /// Aerosol/cloud layer scaling factor for transition layer
  double ds;	      /// Segment length [km].
	double u[NG];	    /// Column density [molecules/cm^2].
#ifdef CURTIS_GODSON
	double cgp[NG];	  /// Curtis-Godson pressure [hPa].
	double cgt[NG];	  /// Curtis-Godson temperature [K].
	double cgu[NG];	  /// Curtis-Godson column density [molecules/cm^2].
#endif
#ifdef GPUDEBUG
	int ip, ir;       /// debug helpers
#endif
} pos_t; //////////////////////////////////////////////////////////////////////

typedef struct { /// Observation geometry and radiance data. //////////////////
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
} obs_t; ///////////////////////////////////////////////////////////////////////

typedef float real_tblND_t;

typedef struct {  /// Transposed emissivity look-up tables. - GPU version  /////
	int32_t np[NG][ND];                             /// Number of pressure levels.
	int32_t nt[NG][TBLNP][ND];                      /// Number of temperatures.
	int32_t nu[NG][TBLNP][TBLNT][ND];               /// Number of column densities.
	double p[NG][TBLNP][ND];                        /// Pressure [hPa].
	double t[NG][TBLNP][TBLNT][ND];                 /// Temperature [K].
	real_tblND_t u[NG][TBLNP][TBLNT][TBLNU][ND];    /// Column density [molecules/cm^2].
	real_tblND_t eps[NG][TBLNP][TBLNT][TBLNU][ND];  /// Emissivity.
	double sr[TBLNS][ND];                           /// Source function radiance [W/(m^2 sr cm^-1)].
	double st[TBLNS];                               /// Source function temperature [K].
#ifdef  FAST_INVERSE_OF_U
	/// u0inv[g][p][t][d] * u[g][p][t][0][d] == 1 must hold!  /// FAST_INVERSE_OF_U
	double u0inv[NG][TBLNP][TBLNT][ND];                       /// FAST_INVERSE_OF_U
  /// We assume a logarithmic increment by 2^(1/6)          /// FAST_INVERSE_OF_U
#endif
} trans_table_t; ///////////////////////////////////////////////////////////////////////

/// --------------------------------- for scattering -----------------------------------
typedef struct { /// Emissivity look-up tables. ////////////////////////////////////////
  int np[NGMAX][NDMAX];                                 /// Number of pressure levels.
  int nt[NGMAX][NDMAX][TBLNPMAX];                       /// Number of temperatures.
  int nu[NGMAX][NDMAX][TBLNPMAX][TBLNTMAX];             /// Number of column densities.
  double p[NGMAX][NDMAX][TBLNPMAX];                     /// Pressure [hPa].
  double t[NGMAX][NDMAX][TBLNPMAX][TBLNTMAX];           /// Temperature [K].
  float u[NGMAX][NDMAX][TBLNPMAX][TBLNTMAX][TBLNUMAX];  /// Column density [molecules/cm^2].
  float eps[NGMAX][NDMAX][TBLNPMAX][TBLNTMAX][TBLNUMAX];/// Emissivity.
  double st[TBLNSMAX];                                  /// Source function temperature [K].
  double sr[NDMAX][TBLNSMAX];                           /// Source function radiance [W/(m^2 sr cm^-1)].
} tbl_t; ///////////////////////////////////////////////////////////////////////////////

typedef struct { /// Line-of-sight data. ///////////////////////////////////////////////
  int np;                 /// Number of LOS points.
  double z[NLOS];         /// Altitude [km].
  double lon[NLOS];       /// Longitude [deg].
  double lat[NLOS];       /// Latitude [deg].
  double p[NLOS];         /// Pressure [hPa].
  double t[NLOS];         /// Temperature [K].
  double q[NLOS][NGMAX];  /// Volume mixing ratio.
  double k[NLOS][NWMAX];  /// Extinction [1/km].
  int aeroi [NLOS];       /// Aerosol/cloud layer index
  double aerofac[NLOS];   /// Aerosol/cloud layer scaling factor for transition layer
  double tsurf;           /// Surface temperature [K].
  double ds[NLOS];        /// Segment length [km].
  double u[NLOS][NGMAX];  /// Column density [molecules/cm^2].
} los_t; ///////////////////////////////////////////////////////////////////////////////

typedef struct { /// Aerosol and Cloud properties. /////////////////////////////////////
  /// Aerosol and cloud input parameters
  int nm;                         /// Number of aerosol/cloud models
  double top_mod[SCAMOD];         /// Model top altitude [km]
  double bottom_mod[SCAMOD];      /// Model bottom altitude [km]
  double trans_mod[SCAMOD];       /// Model transition layer thickness [km]
  char type[SCAMOD][LEN];         /// Optical properties source
  char filepath[SCAMOD][LEN];     /// Refractive index file or optical properties file
  double nn[SCAMOD];              /// Number concentration [cm-3] or extinction coefficient [km-1]
  double rr[SCAMOD];              /// Median radius of log-normal size distribution [mum]
  double ss[SCAMOD];              /// Width of log-normal size distribution
  /// Aerosol and cloud optical properties for radiative transfer
  int nl;                         /// Number of aerosol/cloud layers
  int nmod[NLMAX];                /// Number of modes per layer
  double top[NLMAX];              /// Layer top altitude [km]
  double bottom[NLMAX];           /// Layer bottom altitude [km]
  double trans[NLMAX];            /// Transition layer thickness [km]
  double beta_e[NLMAX][NDMAX];    /// Extinction coefficient [1/km]
  double beta_s[NLMAX][NDMAX];    /// Scattering coefficient [1/km]
  double beta_a[NLMAX][NDMAX];    /// Absorption coefficient [1/km]
  double p[NLMAX][NDMAX][NTHETA]; /// Phase function for each layer, angle and wave number 
 } aero_t; /////////////////////////////////////////////////////////////////////////////

typedef struct { /// Retrieval control parameters. /////////////////////////////////////
  char dir[LEN];            /// Working directory.
  int kernel_recomp;        /// Recomputation of kernel matrix (number of iterations).
  int conv_itmax;           /// Maximum number of iterations.
  double conv_dmin;         /// Minimum normalized step size in state space.
  double resmax;            /// Threshold for radiance residuals [%] (-999 to skip filtering).
  int err_ana;              /// Carry out error analysis (0=no, 1=yes).
  double err_formod[NDMAX]; /// Forward model error [%].
  double err_noise[NDMAX];  /// Noise error [W/(m^2 sr cm^-1)].
  double err_press;         /// Pressure error [%].
  double err_press_cz;      /// Vertical correlation length for pressure error [km].
  double err_press_ch;      /// Horizontal correlation length for pressure error [km].
  double err_temp;          /// Temperature error [K].
  double err_temp_cz;       /// Vertical correlation length for temperature error [km].
  double err_temp_ch;       /// Horizontal correlation length for temperature error [km].
  double err_nn;            /// Particle concentration error [cm-3].
  double err_rr;            /// Particle radius error [m-6].
  double err_ss;            /// Particle size distribution width error.
  double err_q[NGMAX];      /// Volume mixing ratio error [%].
  double err_q_cz[NGMAX];   /// Vertical correlation length for volume mixing ratio error [km].
  double err_q_ch[NGMAX];   /// Horizontal correlation length for volume mixing ratio error [km].
  double err_k[NWMAX];      /// Extinction error [1/km].
  double err_k_cz[NWMAX];   /// Vertical correlation length for extinction error [km]. 
  double err_k_ch[NWMAX];   /// Horizontal correlation length for extinction error [km].
} ret_t; 

#endif
