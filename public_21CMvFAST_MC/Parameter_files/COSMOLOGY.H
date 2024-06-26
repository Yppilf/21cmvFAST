#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include <math.h>
#include "Variables.h"


/****** New in v1.1. ******/
/*** WDM parameters ***/
#define P_CUTOFF (int) 0 // supress the power spectrum? 0= CDM; 1=WDM
#define M_WDM (float) 2 // mass of WDM particle in keV.  this is ignored if P_CUTOFF is set to zero
#define g_x (float) 1.5 /* degrees of freedom of WDM particles; 1.5 for fermions */


// Some cosmological parameters are now obtained from a textfile which is taken as input. Required
// for sampling the cosmological parameters in 21CMMC.

/* COSMOLOGICAL PARAMETERS */
//#define SIGMA8 (0.815) // rms mass fluctuation on 8/hlittle Mpc scales at z=0
//#define hlittle (0.678) // little hubble h
//#define OMm  (0.308) // at z=0
//#define OMl  (float) (1.0-OMm) // at z=0
//#define OMb  (float) ((0.02226/hlittle)/hlittle) // at z=0
#define OMn  (0.0)
#define OMk  (0.0)
#define OMr  (8.6e-5)
#define OMtot (1.0)
#define Y_He (0.245)
//#define POWER_INDEX (0.968) /* power law index for the power spectrum */
#define wl   (-1.0) /* dark energy equation of state parameter (wl = -1 for vacuum ) */

/*
   Note that the best fit b and c ST params for these 3D realisations have a redshift,
   and a DELTA_R_FACTOR (see ANAL_PARAMS.H) dependence, as (will be) shown in Mesinger+.
   For converged mass functions at z~5-10, set DELTA_R_FACTOR=1.1 and SHETH_b~0.15 SHETH_c~0.05 (work in progress)

   For most purposes, a larger step size is quite sufficient and provides an excelent match
   to N-body and smoother mass functions, though the b and c parameters should be changed
   to make up for some "stepping-over" massive collapsed halos (see Mesinger, Perna, Haiman (2005)
   and Mesinger et al., in preparation)

   For example, at z~7-10, one can set DELTA_R_FACTOR=1.3 and SHETH_b=0.15 SHETH_c=0.25, to increase the speed of the halo finder.
*/
#define SHETH_b (0.15) /*  1-D realisation best fit from Barkana et al. 2001: SHETH_b 0.34 */
#define SHETH_c (0.05) /* 1-D realisation best fit from Barkana et al. 2001: SHETH_c 0.81 */

#define Zreion_HeII (double) 3 /* redshift of helium reionization, currently only used for tau_e */



#define OUTPUT_FOLDER "/scratch/s4950836/output_3"
//JBM:output folder name. If you want output in programs (as default) write "."

// TS: Input folder name for boxes, incase files need to get moved.
// If you want the default, write "../Boxes"
#define BOXES_INPUT_FOLDER "/scratch/s4950836/Boxes"

// TS: Where to store the Walker files, as they can take up quite a bit of space
// For default, leave an empty string or "./"
#define WALKER_FOLDER "/scratch/s4950836/Walker/"

// TS: When having issues with development, I added this function to drive_21cmMC_streamlined to not overwrite files 
// when already created in a previous run. 
// Note this function isn't fully implemented yet. If set to 1, default functionality will still apply.
// Note this does not mean the output isn't calculated(Which would be a useful functionality), only not rewritten to file
#define OVERWRITE_OUTPUT (0)


/* STRUCTURE */
#define Deltac (1.68) /* at z=0, density excess at virialization */
#define FILTER (0) /* smoothing: 0=tophat, 1=gaussian */
#define POWER_SPECTRUM (5) /* EH=0 BBKS=1  EFSTATHIOU=2  PEEBLES=3  WHITE=4 */
//JBM: POWER_SPECTRUM=5 reads from CLASS transfer functions//
#define USE_RELATIVE_VELOCITIES (1)
//JBM:whether to include DM-b relative velocities or not. 0 no 1 yes. //we take all relative velocities at z=zi.
#define ZERO_REL_VELOCITY (0)
//JBM: if you want USE_RELATIVE_VELOCITIES but fixing v=0 (so just molecular cooling).
#define AVG_REL_VELOCITY (0)
//JBM: if you want USE_RELATIVE_VELOCITIES but fixing v=vavg (almost same bckgrd as v case but no fluctuations).
#define FEEDBACK_LEVEL (3)
//JBM:level of Ly-Werner feedback. 0 (none), 1 (low) 2 (medium) 3 (high), 4 (very high), 10 (atomic only).
#define VBCMAX 90
#define VBCMIN 0
//JBM: max and min velocity in km/s, we cap above it (to have simple interp. tables)
#define KTOP 2000.0
#define KBOT 1e-5
//JBM:max and min k in Mpc-1 for which we integrate. Fixed to avoid k>k_CLASS.
#define ZLYMANWERNER 10.0
//JBM:redshift at which our F_LW data ends, and we assume feedback saturates to atomic (T_vir input).
#define VCB_RMS      (29.7) 
//JBM:value of rms vcb [<v^2>^(1/2)]
#define VCB_AVG (double)(0.92 * VCB_RMS)
//JBM:the avg value of vcb is slightly smaller than vrms in a MB distribution
#define DO_ANISOTROPY (0)
//JBM:whether to calculate P(k,\mu) or not. Slower if yes.




/*******  END USER CHANGABLE DEFINITIONS ****************/
/******* THINGS BELOW SHOULD NOT REQUIRE CHANGING *******/
/********************************************************/


#define CLASS_FILENAME (const char *) "../External_tables/Transfers_z0.dat"
#define CLASS_LENGTH (150) //length of the CLASS transfer function


#define N_nu (1.0) /* # of heavy neutrinos (for EH trasfer function) */
#define BODE_e (0.361) /* Epsilon parameter in Bode et al. 2000 trans. funct.*/
#define BODE_n (5.0) /* Eda parameter in Bode et al. 2000 trans. funct.*/
#define BODE_v (1.2) /* Nu parameter in Bode et al. 2000 trans. funct.*/
/* SHETH params gotten from barkana et al. 2001, pg. 487 */
#define SHETH_a (0.73) /* Sheth and Tormen a parameter (from Jenkins et al. 2001) */
#define SHETH_p (0.175) /* Sheth and Tormen p parameter (from Jenkins et al. 2001) */
#define SHETH_A (0.353) /* Sheth and Tormen A parameter (from Jenkins et al. 2001) */
/********************************************************/


/* CONSTANTS */
#define LN10 (double) (2.30258509299)
#define SIGMAT (double) (6.6524e-25)  /* Thomson scattering cross section in cm^-2 */
#define SIGMA_HI (double) (6.3e-18)  /* HI ionization  cross section at 13.6 eV in cm^-2 */
#define E (double) (2.71828182846)
#define PI (double) (3.14159265358979323846264338327)
#define TWOPI (double) (2.0*PI)
#define FOURPI (double) (4.0*PI)
#define G (double) 6.67259e-8 /* cm^3 g^-1 s^-2*/
#define hplank (double) 6.62606896e-27 /* erg s */
#define TINY (double) 1e-30
#define FRACT_FLOAT_ERR (double) 1e-7 /* fractional floating point error */
#define f_alpha (float) 0.4162 /* oscillator strength of Lya */
#define Ly_alpha_HZ  (double ) 2.46606727e15  /* frequency of Lyalpha */
#define C  (double) 29979245800.0  /*  speed of light  (cm/s)  */
#define alphaA_10k (double) 4.18e-13 /* taken from osterbrock for T=10000 */
#define alphaB_10k (double) 2.59e-13 /* taken from osterbrock for T=10000 */
#define alphaB_20k (double) 2.52e-13 /* taken from osterbrock for T=20000 */
#define Ly_alpha_ANG (double) 1215.67
#define Ly_beta_ANG (double) 1025.18
#define Ly_gamma_ANG (double) 972.02
#define NV_ANG (double) 1240.81 /* NV line center */
#define CMperMPC (double) 3.086e24 /* cm/Mpc */
#define SperYR (double) 31556925.9747 /* s/yr */
#define Msun (double) 1.989e33 /* g */
#define Rsun (double) 6.9598e10 /* cm */
#define Lsun (double) 3.90e33 /* erg/s */
#define T_cmb (double) 2.728 /* K */
#define k_B (double) 1.380658e-16 /* erg / K */
#define m_p (double) 1.6726231e-24 /* proton mass (g) */
#define m_e (double) 9.10938188e-28 /* electron mass (g) */
#define e_charge (double) 4.8033e-10 /* elemetary charge (esu=g^1/2 cm^3/2 s^-1*/
#define SQDEG_ALLSKY (double) ((360.0*360.0)/PI) /* Square degrees in all sky */
#define G_AB_Jy (double) 3631.0 /* AB mag constant in Jy */
#define NU_over_EV (double) (1.60217646e-12 / hplank)
#define NUIONIZATION (double) (13.60*NU_over_EV)  /* ionization frequency of H */
#define HeII_NUIONIZATION (double) (NUIONIZATION*4) /* ionization frequency of HeII */
#define HeI_NUIONIZATION (double) (24.59*NU_over_EV) /* ionization frequency of HeI */
#define Ho  (double) (hlittle*3.2407e-18) /* s^-1 at z=0 */
#define RHOcrit (double) ( (3.0*Ho*Ho / (8.0*PI*G)) * (CMperMPC*CMperMPC*CMperMPC)/Msun) /* Msun Mpc^-3 */ /* at z=0 */
#define RHOcrit_cgs (double) (3.0*Ho*Ho / (8.0*PI*G)) /* g pcm^-3 */ /* at z=0 */
#define No  (double) (RHOcrit_cgs*OMb*(1-Y_He)/m_p)  /*  current hydrogen number density estimate  (#/cm^3)  ~1.92e-7*/
#define He_No (double) (RHOcrit_cgs*OMb*Y_He/(4.0*m_p)) /*  current helium number density estimate */
#define N_b0 (double) (No+He_No) /* present-day baryon num density, H + He */
#define f_H (double) (No/(No+He_No))  /* hydrogen number fraction */
#define f_He (double) (He_No/(No+He_No))  /* helium number fraction */
#define T21 (double) 0.0628 /* temperature corresponding to the 21cm photon */
#define A10_HYPERFINE (double) 2.85e-15 /* spontaneous emission coefficient in s^-1 */

#define C_KMS (double) C/1e5  /*JBM:  speed of light in km/s  */


#endif
