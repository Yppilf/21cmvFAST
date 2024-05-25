// Header file for inclusions in drive_21cmMC_streamlined.c to prevent double inclusion errors
#ifndef _DRIVE_H_
#define _DRIVE_H_

#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/Variables.h"
#include "init_21cmMC_Ts_arrays.c"
#include "bubble_helper_progs.c"
#include "heating_helper_progs.c"
#include "gsl/gsl_sf_erf.h"
#include "filter.c"
#include <sys/stat.h>

// TS: Additional functions for dynamic freeing of memory during allocation
void* allocate_memory_fftw(size_t size, void** ptr);
void* allocate_memory_calloc(size_t num, size_t size, void** ptr);
void free_pointer_fftw(fftwf_complex** ptr);
void free_pointer_float(float** ptr);
void free_2D_float_array(float*** ptr, size_t dim1, size_t dim2);
void free_pointer_double(double** ptr);
void free_2D_double_array(double*** ptr, size_t dim1, size_t dim2);
void free_3D_double_array(double ****ptr, size_t dim1, size_t dim2, size_t dim3);
void free_pointer_short(short** ptr);
void free_2D_short_array(short*** ptr, size_t dim1, size_t dim2);
void free_memory();

#endif