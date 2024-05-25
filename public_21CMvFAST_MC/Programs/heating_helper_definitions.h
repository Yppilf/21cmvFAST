#ifndef _HEATING_HELPER_H_
#define _HEATING_HELPER_H_

#define NSPEC_MAX (int) 23
#define RECFAST_NPTS (int) 501
#define KAPPA_10_NPTS (int) 27
#define KAPPA_10_elec_NPTS (int) 20
#define KAPPA_10_pH_NPTS (int) 17

#define KAPPA_10_NPTS_Spline (int) 30
#define KAPPA_10_elec_NPTS_Spline (int) 30
#define KAPPA_10_pH_NPTS_Spline (int) 30

#define X_RAY_Tvir_POINTS (int) 100

#define zpp_interp_points (int) (400)                  /* Number of interpolation points for the interpolation table for z'' */
#define dens_Ninterp (int) (400)                       /* Number of interpolation points for the interpolation table for the value of the density field */

#define ZINT_MIN (double) 5.0
#define ZINT_MAX (double) 50.0

#define NZINT (int) 500
#define ZINT_STEP (double) (ZINT_MAX-ZINT_MIN)/(NZINT-1.0)

#define VINT_MIN (double) 0.0
#define VINT_MAX (double) 90.0
#define NVINT (int) 100
#define VINT_STEP (double) (VINT_MAX-VINT_MIN)/(NVINT-1.0)

#define FCOLLMIN (double) 1e-13 //if Fcoll is smaller than this we set it to it to avoid numerical issues

#define FILENAME_FCOLL_VCB  "../External_tables/Fcollapse_table"
#define FILENAME_SIGMACOOL_VCB "../External_tables/sigmacool_table"

#endif