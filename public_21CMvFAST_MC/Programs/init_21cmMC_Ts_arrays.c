#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/Variables.h"
#include "bubble_helper_progs.c"
#include "heating_helper_progs.c"
#include "gsl/gsl_sf_erf.h"
#include "filter.c"
#include <sys/stat.h>

// TS: Adding checks to make sure memory is properly allocated
void* allocate_memory_fftw(size_t size, void** ptr) {
    *ptr = malloc(size);
    if (*ptr == NULL) {
        perror("Failed to allocate memory");
    }
    return *ptr;
}

void* allocate_memory_calloc(size_t num, size_t size, void** ptr) {
    *ptr = calloc(num, size);
    if (*ptr == NULL) {
        perror("Failed to allocate memory");
    }
    return *ptr;
}

// TS: Make sure memory is properly freed if issues arise
void free_pointer_fftw(fftwf_complex** ptr) {
    if (*ptr != NULL) {
        fftwf_free(*ptr);
        *ptr = NULL;
    }
}

void free_pointer_float(float** ptr) {
    if (*ptr != NULL) {
        free(*ptr);
        *ptr = NULL;
    }
}

void free_2D_float_array(float*** ptr, size_t dim1, size_t dim2) {
    if (ptr == NULL || *ptr == NULL) {
        return; // Pointer is already NULL, nothing to free
    }

    // Free memory for each row
    for (size_t i = 0; i < dim1; i++) {
        free_pointer_float(&((*ptr)[i]));
    }

    // Free memory for the array of pointers itself
    free(*ptr);
    *ptr = NULL;
}

void free_pointer_double(double** ptr) {
    if (*ptr != NULL) {
        free(*ptr);
        *ptr = NULL;
    }
}

void free_2D_double_array(double*** ptr, size_t dim1, size_t dim2) {
    if (ptr == NULL || *ptr == NULL) {
        return; // Pointer is already NULL, nothing to free
    }

    // Free memory for each row
    for (size_t i = 0; i < dim1; i++) {
        free_pointer_double(&((*ptr)[i]));
    }

    // Free memory for the array of pointers itself
    free(*ptr);
    *ptr = NULL;
}

void free_3D_double_array(double ****ptr, size_t dim1, size_t dim2, size_t dim3) {
    if (ptr == NULL || *ptr == NULL) {
        return;
    }

    // Free memory for each row and each 2D array
    for (size_t i = 0; i < dim1; i++) {
        // Free memory for the 2D array of pointers
        free_2D_double_array(&((*ptr)[i]), dim2, dim3);
    }

    // Free memory for the array of pointers itself
    free(*ptr);
    *ptr = NULL;
}

void free_pointer_short(short** ptr) {
    if (*ptr != NULL) {
        free(*ptr);
        *ptr = NULL;
    }
}

void free_2D_short_array(short*** ptr, size_t dim1, size_t dim2) {
    if (ptr == NULL || *ptr == NULL) {
        return; // Pointer is already NULL, nothing to free
    }

    // Free memory for each row
    for (size_t i = 0; i < dim1; i++) {
        free_pointer_short(&((*ptr)[i]));
    }

    // Free memory for the array of pointers itself
    free(*ptr);
    *ptr = NULL;
}

void free_memory() {
    free_pointer_fftw(&box);
    free_pointer_fftw(&unfiltered_box);
    free_pointer_fftw(&box_vcb);
    free_pointer_fftw(&unfiltered_vcb_box);
    free_pointer_float(&Tk_box);
    free_pointer_float(&x_e_box);
    free_pointer_float(&Ts);
    free(inverse_diff);
    free_pointer_float(&zpp_growth);
    free_3D_double_array(&fcoll_R_grid, NUM_FILTER_STEPS_FOR_Ts, zpp_interp_points, dens_Ninterp);
    free_3D_double_array(&dfcoll_dz_grid, NUM_FILTER_STEPS_FOR_Ts, zpp_interp_points, dens_Ninterp);
    free(fcoll_R_array);
    free_pointer_double(&Sigma_Tmin_grid);
    free_2D_double_array(&grid_dens, NUM_FILTER_STEPS_FOR_Ts, dens_Ninterp);
    free_2D_double_array(&density_gridpoints, dens_Ninterp, NUM_FILTER_STEPS_FOR_Ts);
    free_pointer_double(&ST_over_PS_arg_grid);
    free_2D_double_array(&logFcoll_vcb, NZINT, NVINT);
    free_2D_double_array(&sigmacool_vcb, NZINT, NVINT);
    free_2D_short_array(&dens_grid_int_vals, HII_TOT_NUM_PIXELS, NUM_FILTER_STEPS_FOR_Ts);
    free_2D_float_array(&delNL0_rev, HII_TOT_NUM_PIXELS, NUM_FILTER_STEPS_FOR_Ts);
    free_2D_float_array(&vcb_rev, HII_TOT_NUM_PIXELS, NUM_FILTER_STEPS_FOR_Ts);
    free_2D_double_array(&fcoll_interp1, dens_Ninterp, NUM_FILTER_STEPS_FOR_Ts);
    free_2D_double_array(&fcoll_interp2, dens_Ninterp, NUM_FILTER_STEPS_FOR_Ts);
    free_2D_double_array(&dfcoll_interp1, dens_Ninterp, NUM_FILTER_STEPS_FOR_Ts);
    free_2D_double_array(&dfcoll_interp2, dens_Ninterp, NUM_FILTER_STEPS_FOR_Ts);
    free_pointer_double(&zpp_edge);
    free_pointer_double(&sigma_atR);
    free_pointer_double(&sigma_Tmin);
    free_pointer_double(&ST_over_PS);
    free_pointer_double(&sum_lyn);
    free_pointer_float(&zpp_for_evolve_list);
    free_pointer_float(&R_values);
    free_pointer_float(&SingleVal_float);
    free_pointer_float(&delNL0_bw);
    free_pointer_float(&delNL0_Offset);
    free_pointer_float(&delNL0_LL);
    free_pointer_float(&delNL0_UL);
    free_pointer_float(&delNL0_ibw);
    free_pointer_float(&log10delNL0_diff);
    free_pointer_float(&log10delNL0_diff_UL);
    free_2D_double_array(&freq_int_heat_tbl,x_int_NXHII,NUM_FILTER_STEPS_FOR_Ts);
    free_2D_double_array(&freq_int_ion_tbl,x_int_NXHII,NUM_FILTER_STEPS_FOR_Ts);
    free_2D_double_array(&freq_int_lya_tbl,x_int_NXHII,NUM_FILTER_STEPS_FOR_Ts);
    free_2D_double_array(&freq_int_heat_tbl_diff,x_int_NXHII,NUM_FILTER_STEPS_FOR_Ts);
    free_2D_double_array(&freq_int_ion_tbl_diff,x_int_NXHII,NUM_FILTER_STEPS_FOR_Ts);
    free_2D_double_array(&freq_int_lya_tbl_diff,x_int_NXHII,NUM_FILTER_STEPS_FOR_Ts);
    free_pointer_double(&dstarlya_dt_prefactor);
    free_pointer_short(&SingleVal_int);
    // Repeat for other global allocations...
}

// TS: New version of the function, originally from drive_21cmMC_streamlined.c
// This checks if memory is properly allocated
void init_21cmMC_Ts_arrays() {
    int i,j;

    // Initialize the pointers to NULL to ensure a known state before allocation attempt
    box = NULL;
    unfiltered_box = NULL;
    box_vcb = NULL;
    unfiltered_vcb_box = NULL;
    Tk_box = NULL;
    x_e_box = NULL;
    Ts = NULL;
    inverse_diff = NULL;
    zpp_growth = NULL;
    fcoll_R_grid = NULL;
    dfcoll_dz_grid = NULL;
    fcoll_R_array = NULL;
    Sigma_Tmin_grid = NULL;
    grid_dens = NULL;
    density_gridpoints = NULL;
    ST_over_PS_arg_grid = NULL;
    logFcoll_vcb = NULL;
    sigmacool_vcb = NULL;
    dens_grid_int_vals = NULL;
    delNL0_rev = NULL;
    vcb_rev = NULL;
    fcoll_interp1 = NULL;
    fcoll_interp2 = NULL;
    dfcoll_interp1 = NULL;
    dfcoll_interp2 = NULL;
    zpp_edge = NULL;
    sigma_atR = NULL;
    sigma_Tmin = NULL;
    ST_over_PS = NULL;
    sum_lyn = NULL;
    zpp_for_evolve_list = NULL;
    R_values = NULL;
    SingleVal_float = NULL;
    delNL0_bw = NULL;
    delNL0_Offset = NULL;
    delNL0_LL = NULL;
    delNL0_UL = NULL;
    delNL0_ibw = NULL;
    log10delNL0_diff = NULL;
    log10delNL0_diff_UL = NULL;
    freq_int_heat_tbl = NULL;
    freq_int_ion_tbl = NULL;
    freq_int_lya_tbl = NULL;
    freq_int_heat_tbl_diff = NULL;
    freq_int_ion_tbl_diff = NULL;
    freq_int_lya_tbl_diff = NULL;
    dstarlya_dt_prefactor = NULL;
    SingleVal_int = NULL;

    // box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    if (allocate_memory_fftw(sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS, (void**)&box) == NULL) {
        free_memory();
        return;
    }

    if (allocate_memory_fftw(sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS, (void**)&unfiltered_box) == NULL) {
        free_memory();
        return;
    }
    //JBM:
    if (allocate_memory_fftw(sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS, (void**)&box_vcb) == NULL) {
        free_memory();
        return;
    }

    if (allocate_memory_fftw(sizeof(fftwf_complex) * HII_KSPACE_NUM_PIXELS, (void**)&unfiltered_vcb_box) == NULL) {
        free_memory();
        return;
    }

    if (allocate_memory_calloc(HII_TOT_NUM_PIXELS, sizeof(float), (void**)&Tk_box) == NULL){
        free_memory();
        return;
    }

    if (allocate_memory_calloc(HII_TOT_NUM_PIXELS, sizeof(float), (void**)&x_e_box) == NULL){
        free_memory();
        return;
    }

    if (allocate_memory_calloc(HII_TOT_NUM_PIXELS, sizeof(float), (void**)&Ts) == NULL){
        free_memory();
        return;
    }

    // TS: This one originally was allocated in a different way, so I'll handle it separately
    inverse_diff = calloc(x_int_NXHII,sizeof(float));
    if (inverse_diff == NULL) {
        free_memory();
        return;
    }

    if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(float), (void**)&zpp_growth) == NULL){
        free_memory();
        return;
    }

    // Dynamically allocate the 3D arrays and check at each step
    if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(double**), (void**)&fcoll_R_grid) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<NUM_FILTER_STEPS_FOR_Ts;i++){
        if (allocate_memory_calloc(zpp_interp_points, sizeof(double*), (void**)&fcoll_R_grid[i]) == NULL){
            free_memory();
            return;
        }
        for(j=0;j<zpp_interp_points;j++) {
            if (allocate_memory_calloc(dens_Ninterp, sizeof(double), (void**)&fcoll_R_grid[i][j]) == NULL){
                free_memory();
                return;
            }
        }
    }

    if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(double**), (void**)&dfcoll_dz_grid) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<NUM_FILTER_STEPS_FOR_Ts;i++){
        if (allocate_memory_calloc(zpp_interp_points, sizeof(double*), (void**)&dfcoll_dz_grid[i]) == NULL){
            free_memory();
            return;
        }
        for(j=0;j<zpp_interp_points;j++) {
            if (allocate_memory_calloc(dens_Ninterp, sizeof(double), (void**)&dfcoll_dz_grid[i][j]) == NULL){
                free_memory();
                return;
            }
        }
    }

    // TS: This one originally was allocated in a different way, so I'll handle it separately
    fcoll_R_array = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
    if (fcoll_R_array == NULL) {
        free_memory();
        return;
    }

    if (allocate_memory_calloc(zpp_interp_points, sizeof(double), (void**)&Sigma_Tmin_grid) == NULL){
        free_memory();
        return;
    }

    if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(double*), (void**)&grid_dens) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<NUM_FILTER_STEPS_FOR_Ts;i++) {
        if (allocate_memory_calloc(dens_Ninterp, sizeof(double), (void**)&grid_dens[i]) == NULL){
            free_memory();
            return;
        }
    }

    if (allocate_memory_calloc(dens_Ninterp, sizeof(double*), (void**)&density_gridpoints) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<dens_Ninterp;i++) {
        if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(double), (void**)&density_gridpoints[i]) == NULL){
            free_memory();
            return;
        }
    }

    if (allocate_memory_calloc(zpp_interp_points, sizeof(double), (void**)&ST_over_PS_arg_grid) == NULL){
        free_memory();
        return;
    }

    //JBM: we define the F_ST(v) 2D array.
    if (allocate_memory_calloc(NZINT, sizeof(double*), (void**)&logFcoll_vcb) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<NZINT;i++) {
        if (allocate_memory_calloc(NVINT, sizeof(double), (void**)&logFcoll_vcb[i]) == NULL){
            free_memory();
            return;
        }
    }

    //JBM: and sigma(z,v,M_cool(z,v)).
    if (allocate_memory_calloc(NZINT, sizeof(double*), (void**)&sigmacool_vcb) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<NZINT;i++) {
        if (allocate_memory_calloc(NVINT, sizeof(double), (void**)&sigmacool_vcb[i]) == NULL){
            free_memory();
            return;
        }
    }

    if (allocate_memory_calloc(HII_TOT_NUM_PIXELS, sizeof(short*), (void**)&dens_grid_int_vals) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<HII_TOT_NUM_PIXELS;i++) {
        if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(short), (void**)&dens_grid_int_vals[i]) == NULL){
            free_memory();
            return;
        }
    }

    if (allocate_memory_calloc(HII_TOT_NUM_PIXELS, sizeof(float*), (void**)&delNL0_rev) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<HII_TOT_NUM_PIXELS;i++) {
        if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(float), (void**)&delNL0_rev[i]) == NULL){
            free_memory();
            return;
        }
    }

    //JBM:values of velocity smoothed over some scale R
    if (allocate_memory_calloc(HII_TOT_NUM_PIXELS, sizeof(float*), (void**)&vcb_rev) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<HII_TOT_NUM_PIXELS;i++) {
        if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(float), (void**)&vcb_rev[i]) == NULL){
            free_memory();
            return;
        }
    }
    
    if (allocate_memory_calloc(dens_Ninterp, sizeof(double*), (void**)&fcoll_interp1) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<HII_TOT_NUM_PIXELS;i++) {
        if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(double), (void**)&fcoll_interp1[i]) == NULL){
            free_memory();
            return;
        }
    }

    if (allocate_memory_calloc(dens_Ninterp, sizeof(double*), (void**)&fcoll_interp2) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<HII_TOT_NUM_PIXELS;i++) {
        if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(double), (void**)&fcoll_interp2[i]) == NULL){
            free_memory();
            return;
        }
    }

    if (allocate_memory_calloc(dens_Ninterp, sizeof(double*), (void**)&dfcoll_interp1) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<HII_TOT_NUM_PIXELS;i++) {
        if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(double), (void**)&dfcoll_interp1[i]) == NULL){
            free_memory();
            return;
        }
    }

    if (allocate_memory_calloc(dens_Ninterp, sizeof(double*), (void**)&dfcoll_interp2) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<HII_TOT_NUM_PIXELS;i++) {
        if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(double), (void**)&dfcoll_interp2[i]) == NULL){
            free_memory();
            return;
        }
    }

    zpp_edge = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
    if (zpp_edge == NULL) {
        free_memory();
        return;
    }

    sigma_atR = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
    if (sigma_atR == NULL) {
        free_memory();
        return;
    }

    sigma_Tmin = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
    if (sigma_Tmin == NULL) {
        free_memory();
        return;
    }

    ST_over_PS = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
    if (ST_over_PS == NULL) {
        free_memory();
        return;
    }

    sum_lyn = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
    if (sum_lyn == NULL) {
        free_memory();
        return;
    }

    zpp_for_evolve_list = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    if (zpp_for_evolve_list == NULL) {
        free_memory();
        return;
    }

    R_values = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    if (R_values == NULL) {
        free_memory();
        return;
    }

    SingleVal_float = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    if (SingleVal_float == NULL) {
        free_memory();
        return;
    }

    delNL0_bw = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    if (delNL0_bw == NULL) {
        free_memory();
        return;
    }

    delNL0_Offset = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    if (delNL0_Offset == NULL) {
        free_memory();
        return;
    }

    delNL0_LL = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    if (delNL0_LL == NULL) {
        free_memory();
        return;
    }

    delNL0_UL = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    if (delNL0_UL == NULL) {
        free_memory();
        return;
    }

    delNL0_ibw = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    if (delNL0_ibw == NULL) {
        free_memory();
        return;
    }

    log10delNL0_diff = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    if (log10delNL0_diff == NULL) {
        free_memory();
        return;
    }

    log10delNL0_diff_UL = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
    if (log10delNL0_diff_UL == NULL) {
        free_memory();
        return;
    }

    if (allocate_memory_calloc(x_int_NXHII, sizeof(double*), (void**)&freq_int_heat_tbl) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<x_int_NXHII;i++) {
        if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(double), (void**)&freq_int_heat_tbl[i]) == NULL){
            free_memory();
            return;
        }
    }

    if (allocate_memory_calloc(x_int_NXHII, sizeof(double*), (void**)&freq_int_ion_tbl) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<x_int_NXHII;i++) {
        if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(double), (void**)&freq_int_ion_tbl[i]) == NULL){
            free_memory();
            return;
        }
    }

    if (allocate_memory_calloc(x_int_NXHII, sizeof(double*), (void**)&freq_int_lya_tbl) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<x_int_NXHII;i++) {
        if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(double), (void**)&freq_int_lya_tbl[i]) == NULL){
            free_memory();
            return;
        }
    }

    if (allocate_memory_calloc(x_int_NXHII, sizeof(double*), (void**)&freq_int_heat_tbl_diff) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<x_int_NXHII;i++) {
        if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(double), (void**)&freq_int_heat_tbl_diff[i]) == NULL){
            free_memory();
            return;
        }
    }

    if (allocate_memory_calloc(x_int_NXHII, sizeof(double*), (void**)&freq_int_ion_tbl_diff) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<x_int_NXHII;i++) {
        if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(double), (void**)&freq_int_ion_tbl_diff[i]) == NULL){
            free_memory();
            return;
        }
    }

    if (allocate_memory_calloc(x_int_NXHII, sizeof(double*), (void**)&freq_int_lya_tbl_diff) == NULL){
        free_memory();
        return;
    }
    for(i=0;i<x_int_NXHII;i++) {
        if (allocate_memory_calloc(NUM_FILTER_STEPS_FOR_Ts, sizeof(double), (void**)&freq_int_lya_tbl_diff[i]) == NULL){
            free_memory();
            return;
        }
    }

    dstarlya_dt_prefactor = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
    if (dstarlya_dt_prefactor == NULL) {
        free_memory();
        return;
    }

    SingleVal_int = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(short));
    if (SingleVal_int == NULL) {
        free_memory();
        return;
    }
}