// #include "../Parameter_files/INIT_PARAMS.H"
// #include "../Parameter_files/ANAL_PARAMS.H"
// #include "../Parameter_files/Variables.h"
// #include "init_21cmMC_Ts_arrays.c"
// #include "bubble_helper_progs.c"
// #include "heating_helper_progs.c"
// #include "gsl/gsl_sf_erf.h"
// #include "filter.c"
// #include <sys/stat.h>
#include "drive_21cmMC_streamlined.h"


/*

 This is the main file for 21CMMC. This combines Ts.c, find_HII_bubbles.c, delta_T.c and redshift_interpolate_boxes.c (additionally, it includes init.c and perturb_field.c for varying the cosmology)
 from 21cmFAST

 Author: Brad Greig (July 7th 2017). This is effectively an entire re-write of 21cmFAST.

 It is called from command line, with a fixed number of arguments (order is important). There is basically no error checking, as it would be too complicated to use that within the MCMC.

 An example command line call: ./drive_21cmMC_streamlined 1.000000 1.000000 0 1 0 6.0

 First two indices are required for opening the Walker_ID1_ID2.txt and WalkerCosmology_ID1_ID2.txt which contain all the cosmology and astrophysical parameters (example included)

 Third argument: (0 or N), contains the number of co-eval redshifts to be sampled (will be ignored if the light-cone option is set)

 Fourth argument: (0 or 1), calculate the light-cone (1) or co-eval (0) redshifts

 Fifth argument: (0 or 1), whether or not to include a power-law index for the ionising efficiency (this option is only valid for the co-eval boxes. Does not work for calculating the IGM spin temperature
                    or the light-cone box). Have to update the Ts.c part of the code to allow this for the light-cone

 Sixth argument: Redshift to which Ts.c is evolved down to

*/

/* Throughout this and other 21cmMC drivers, the use of Deltac is not for checking against
 the z=0 collapsed fraction, but rather, it is used as a threshold for the validity of the
 collapse fraction expression. This expression is only valid up to Deltac
 */

// For storing the 21cm PS light cone filenames to be able to write them to file to be read by the MCMC sampler
char lightcone_box_names[1000][500];

float REDSHIFT;

void init_21cmMC_Ts_arrays();
void init_21cmMC_HII_arrays();
void init_21cmMC_TsSaveBoxes_arrays();

void ComputeBoxesForFile();
void ComputeTsBoxes();
void free_TS_memory(); // Added this to avoid double free errors until after the data is obtained
void ComputeIonisationBoxes(int sample_index, float REDSHIFT_SAMPLE, float PREV_REDSHIFT);

void adj_complex_conj();
void ComputeInitialConditions();
void ComputePerturbField(float REDSHIFT_SAMPLE);
void GeneratePS(int CO_EVAL, double AverageTb);
void GeneratePS_aniso(int CO_EVAL, double AverageTb, char* filename);



void GeneratePS_vel();//JBM:test function to output the power spectrum of the v_cb box.
float *smoothed_box_vcb;


void ReadFcollTable();

void destroy_21cmMC_Ts_arrays();
void destroy_21cmMC_HII_arrays();
void destroy_21cmMC_TsSaveBoxes_arrays();

// This, and the one below are functions for determining the correct cell positions for direction of the light-cone. Tested this for the z-direction, but should be valid for all.
// Note that there is no option for FLIP_BOXES as we want to mimic the observed light-cone
unsigned long long coeval_box_pos(int LOS_dir,int xi,int yi,int zi){
    unsigned long long position;

    switch(LOS_dir) {
        case 0:
            position = HII_R_INDEX(zi, xi, yi);
            break;
        case 1:
            position = HII_R_INDEX(xi, zi, yi);
            break;
        case 2:
            position = HII_R_INDEX(xi, yi, zi);
            break;
    }
    return position;
}

unsigned long long coeval_box_pos_FFT(int LOS_dir,int xi,int yi,int zi){
    unsigned long long position;

    switch(LOS_dir) {
        case 0:
            position = HII_R_FFT_INDEX(zi, xi, yi);
            break;
        case 1:
            position = HII_R_FFT_INDEX(xi, zi, yi);
            break;
        case 2:
            position = HII_R_FFT_INDEX(xi, yi, zi);
            break;
    }
    return position;
}

int file_exists(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file) {
        fclose(file);
        return 1;
    }
    return 0;
}

// TS: This file had issues with double freeing of memory
// Helper macro to safely free memory and set pointer to NULL
#define SAFE_FREE(ptr) if(ptr) { free(ptr); ptr = NULL; }
#define SAFE_FFTWF_FREE(ptr) if(ptr) { fftwf_free(ptr); ptr = NULL; }

// TS: Note memory is not properly freed and an error can be expected at freeing, after data is obtained
int main(int argc, char ** argv){
     fprintf(stderr, "\n----------------------------\nStarting drive_21cmMC_streamlined.c\n----------------------------\n");


    // The standard build of 21cmFAST requires openmp for the FFTs. 21CMMC does not, however, for some computing architectures, I found it important to include this
    omp_set_num_threads(1);

    char filename[500];
    char dummy_string[500];
    FILE *F;

    sprintf(dummy_string,"mkdir %s",OUTPUT_FOLDER);
    struct stat st = {0};
    if (stat(OUTPUT_FOLDER, &st) == -1) {
        system(dummy_string);
    }

    LC_BOX_PADDING = (int)ceil(LC_BOX_PADDING_IN_MPC/((float)BOX_LEN*(float)HII_DIM));

    int i,j,k,temp_int,temp_int2, counter;
    float z_prime,prev_z_prime;


    unsigned long long ct;

    // Setting it to its maximum value
    INHOMO_RECO_R_BUBBLE_MAX = 50.0;

    // All parameters, redshifts and flag options are written to file by the python MCMC driver

    // Take as input:
    // 1) Random ID for MCMC walker
    // 2) Number of parameters set by the user
    // 3) Number of redshifts required (number of co-eval cubes to be used. Will be overwritten if the lightcone option is set)

    // Assign random walker ID (takes two values as a fail-safe. First is a random number generated by the MCMC algorithm, the second is the first parameter of the specific walker)
    INDIVIDUAL_ID = atof(argv[1]);
    INDIVIDUAL_ID_2 = atof(argv[2]);

    // Number of user defined redshifts for which find_HII_bubbles will be called
    N_USER_REDSHIFT = atof(argv[3]);

    // Flag set to 1 if light cone boxes are to be used (feature has yet to be added)
    USE_LIGHTCONE = atof(argv[4]);

    /****** NOTE: Need to add in a flag here, which toggles including alpha as a useable parameter ******/
    /****** In doing it here, it enables the majority of the remaining code to be relatively straight-forward (i.e. doesn't change existing text-file structure etc.) ******/
    /****** This hasn't been rigorously checked yet. Need to look into this at some point... *******/
    INCLUDE_ZETA_PL = atof(argv[5]);

    // Redshift for which Ts.c is evolved down to, i.e. z'
    REDSHIFT = atof(argv[6]);

    // Determines the lenght of the walker file, given the values set by TOTAL_AVAILABLE_PARAMS in Variables.h and the number of redshifts
    if(USE_LIGHTCONE) {
        WALKER_FILE_LENGTH = TOTAL_AVAILABLE_PARAMS + 1;
    }
    else {
        WALKER_FILE_LENGTH = N_USER_REDSHIFT + TOTAL_AVAILABLE_PARAMS + 1;
    }

    // Create arrays to read in all the parameter data from the two separate walker files
    double *PARAM_COSMOLOGY_VALS = calloc(TOTAL_COSMOLOGY_FILEPARAMS,sizeof(double));
    double *PARAM_VALS = calloc(TOTAL_AVAILABLE_PARAMS,sizeof(double));

    /////////////////   Read in the cosmological parameter data     /////////////////

    sprintf(filename,"%sWalkerCosmology_%1.6lf_%1.6lf.txt",WALKER_FOLDER,INDIVIDUAL_ID,INDIVIDUAL_ID_2);
    fprintf(stderr, "Printing to %s\n", filename);
    F = fopen(filename,"rt");

    for(i=0;i<TOTAL_COSMOLOGY_FILEPARAMS;i++) {
        fscanf(F,"%s\t%lf\n",&dummy_string,&PARAM_COSMOLOGY_VALS[i]);
    }
    fclose(F);

    // Assign these values. Hard-coded, so order is important
    RANDOM_SEED = (unsigned long long)PARAM_COSMOLOGY_VALS[0];
    SIGMA8 = (float)PARAM_COSMOLOGY_VALS[1];
    hlittle = (float)PARAM_COSMOLOGY_VALS[2];
    OMm = (float)PARAM_COSMOLOGY_VALS[3];
    OMl = (float)PARAM_COSMOLOGY_VALS[4];
    OMb = (float)PARAM_COSMOLOGY_VALS[5];
    POWER_INDEX = (float)PARAM_COSMOLOGY_VALS[6]; //power law on the spectral index, ns




    /////////////////   Read in the astrophysical parameter data     /////////////////

    // Determine length of parameter file to read in
    // All available parameters to be varied in the MCMC are always listed, but are toggled on/off using 1/0
    // The MCMC sets the toggle, this C file reads the toggle and uses/sets the parameter values appropriately

    sprintf(filename,"%sWalker_%1.6lf_%1.6lf.txt",WALKER_FOLDER, INDIVIDUAL_ID,INDIVIDUAL_ID_2);
    fprintf(stderr, "\n----------------------------\nReading parameter data from %s\n----------------------------\n", filename);
    F = fopen(filename,"rt");

    if(!INHOMO_RECO||!USE_LIGHTCONE) {
        redshifts = calloc(N_USER_REDSHIFT,sizeof(double));
    }

    temp_int = 0;
    temp_int2 = 0;
    for(i=0;i<WALKER_FILE_LENGTH;i++) {
        if(i==0) {
            fscanf(F,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",&dummy_string,&GenerateNewICs,&SUBCELL_RSD,&USE_FCOLL_IONISATION_TABLE,&SHORTEN_FCOLL,&USE_TS_FLUCT,&INHOMO_RECO,&STORE_DATA);
        }
        else if (i > 0 && i <= TOTAL_AVAILABLE_PARAMS) {
            fscanf(F,"%s\t%lf\n",&dummy_string,&PARAM_VALS[temp_int]);
            temp_int += 1;
        }
        else {
            if(!INHOMO_RECO) {
                fscanf(F,"%s\t%lf\n",&dummy_string,&redshifts[temp_int2]);
                temp_int2 += 1;
            }
            else {
                // Since the INHOMO_RECO flag has been set, we need to take all the redshifts that are used for the Ts.c part of the calculation (in order to have the right
                // time-steps for tracking the recombinations.
                continue;
            }
        }
    }
    fclose(F);


    // GenerateNewICs: Whether to create a new density field at each sampling (i.e. new initial conditions). Must use if the cosmology is being varied
    // SUBCELL_RSD: Whether to include redshift space distortions along the line-of-sight (z-direction only).
    // USE_FCOLL_IONISATION_TABLE: Whether to use an interpolation for the collapsed fraction for the find_HII_bubbles part of the computation
    // SHORTEN_FCOLL: Whether to use an interpolation for the collapsed fraction for the Ts.c computation
    // USE_TS_FLUCT: Whether to perform the full evolution of the IGM spin temperature, or just assume the saturated spin temperature limit
    // INHOMO_RECO: Whether to include inhomogeneous recombinations into the calculation of the ionisation fraction
    // STORE_DATA: Whether to output the global data for the IGM neutral fraction and average temperature brightness (used for the global signal)


    // Initialise the power spectrum data, and relevant functions etc., for the entire file here (i.e. it is only done once here)
    init_ps();

    // If the USE_LIGHTCONE option is set, need to determing the size of the entire line-of-sight dimension for storing the slice indexes and corresponding reshifts per slice
    dR = (BOX_LEN / (double) HII_DIM) * CMperMPC; // size of cell (in comoving cm)

    if(USE_LIGHTCONE||INHOMO_RECO) {
        // Determine the number of redshifts within the Ts.c calculation to set N_USER_REDSHIFT for the light-cone version of the computation.

        counter = 0;

        z_prime = REDSHIFT*1.0001; //higher for rounding
        while (z_prime < Z_HEAT_MAX)
            z_prime = ((1.+z_prime)*ZPRIME_STEP_FACTOR - 1.);
        prev_z_prime = Z_HEAT_MAX;
        z_prime = ((1.+z_prime)/ ZPRIME_STEP_FACTOR - 1.);

        while (z_prime > REDSHIFT){
            counter += 1;

            prev_z_prime = z_prime;
            z_prime = ((1.+prev_z_prime) / ZPRIME_STEP_FACTOR - 1.);

        }

        N_USER_REDSHIFT = counter;
        // Number of redshifts for boxes used to construct the light-cone. Light-cone ends at final redshift, final box of light-cone is the penultimate box linear interpolated to the final redshift
        N_USER_REDSHIFT_LC = counter - 1;

        if(USE_LIGHTCONE) {
            redshifts_LC = calloc(N_USER_REDSHIFT_LC,sizeof(double));
            start_index_LC = calloc(N_USER_REDSHIFT_LC,sizeof(int));
            end_index_LC = calloc(N_USER_REDSHIFT_LC,sizeof(int));
        }
    }

    // Hard coded to 100,000. Should never excede this, unless very high resolution boxes are being used! (200^3, from z_min = 6 to z_max (z = 35) corresponds to 2232 indices).
    full_index_LC = calloc(100000,sizeof(int));
    slice_redshifts = calloc(100000,sizeof(double));

    if(INHOMO_RECO||USE_LIGHTCONE) {
        redshifts = calloc(N_USER_REDSHIFT,sizeof(double));
    }


    // Some very rudimentary conditionals (the python script should catch them. But, add them here in the C code for testing purposes
    if(GenerateNewICs==1 && (USE_FCOLL_IONISATION_TABLE==1 || SHORTEN_FCOLL==1)) {
        printf("\n");
        printf("Cannot use interpolation tables when generating new initial conditions on the fly\n");
        printf("(Interpolation tables are only valid for a single cosmology/initial condition)\n");
        printf("\n");
        printf("Exiting...");
        printf("\n");
        // Probably should do free the memory properly here...
        return 0;
    }

    if(USE_TS_FLUCT==1 && INCLUDE_ZETA_PL==1) {
        printf("\n");
        printf("Cannot use a non-constant ionising efficiency (zeta) in conjuction with the IGM spin temperature part of the code.\n");
        printf("(This will be changed in future)\n");
        printf("\n");
        printf("Exiting...");
        printf("\n");
        // Probably should do free the memory properly here...
        return 0;
    }

    if(USE_FCOLL_IONISATION_TABLE==1 && INHOMO_RECO==1) {
        printf("\n");
        printf("Cannot use the f_coll interpolation table for find_hii_bubbles with inhomogeneous recombinations\n");
        printf("\n");
        printf("Exiting...");
        printf("\n");
        // Probably should do free the memory properly here...
        return 0;
    }

    if(INHOMO_RECO==1 && USE_TS_FLUCT==0) {
        printf("\n");
        printf("Inhomogeneous recombinations have been set, but the spin temperature is not being computed.\n");
        printf("Inhomogeneous recombinations can only be used in combination with the spin temperature calculation (different from 21cmFAST).\n");
        printf("\n");
        printf("Exiting...");
        printf("\n");
        // Probably should do free the memory properly here...
        return 0;
    }

    ///////////////// Hard coded assignment of parameters, but can't do much about it (problem of merging C and Python code) //////////////////////////////////
    EFF_FACTOR_PL_INDEX = PARAM_VALS[0];
    HII_EFF_FACTOR = PARAM_VALS[1];

    // If inhomogeneous recombinations are set, need to switch to an upper limit on the maximum bubble horizon (this is set above).
    // The default choice is chosen to be 50 Mpc, as is default in 21cmFAST.
    if(INHOMO_RECO) {
        R_BUBBLE_MAX = INHOMO_RECO_R_BUBBLE_MAX;
    }
    else {
        R_BUBBLE_MAX = PARAM_VALS[2];
    }

    ION_Tvir_MIN = pow(10.,PARAM_VALS[3]);
    L_X = pow(10.,PARAM_VALS[4]);
    NU_X_THRESH = PARAM_VALS[5];
    NU_X_BAND_MAX = PARAM_VALS[6];
    NU_X_MAX = PARAM_VALS[7];
    X_RAY_SPEC_INDEX = PARAM_VALS[8];
    X_RAY_Tvir_MIN = pow(10.,PARAM_VALS[9]);
    X_RAY_Tvir_LOWERBOUND = PARAM_VALS[10];
    X_RAY_Tvir_UPPERBOUND = PARAM_VALS[11];
    F_STAR = PARAM_VALS[12];
    t_STAR = PARAM_VALS[13];
    N_RSD_STEPS = (int)PARAM_VALS[14];
    LOS_direction = (int)PARAM_VALS[15];

    // Converts the variables from eV into the quantities that Ts.c is familiar with
    NU_X_THRESH *= NU_over_EV;
    NU_X_BAND_MAX *= NU_over_EV;
    NU_X_MAX *= NU_over_EV;

    /////////////////   Populating requisite arrays for the construction of the light-cone box (including the indexing and individual slice redshifts etc.     /////////////////

    // Direction for LOS for light-cone (0 = x, 1 = y, 2 = z). Store the original starting LOS_direction
    Original_LOS_direction = LOS_direction;
    // This is set to 2 (z-direction) for the case when USE_FCOLL_IONISATION_TABLE but the full box is required anyway (i.e. cannot reduce the computation for the light-cone code).
    Default_LOS_direction = 2;

    Stored_LOS_direction_state_1 = LOS_direction;

    if(USE_LIGHTCONE||INHOMO_RECO) {

        counter = 0;

        z_prime = REDSHIFT*1.0001; //higher for rounding
        while (z_prime < Z_HEAT_MAX)
            z_prime = ((1.+z_prime)*ZPRIME_STEP_FACTOR - 1.);
        prev_z_prime = Z_HEAT_MAX;
        z_prime = ((1.+z_prime)/ ZPRIME_STEP_FACTOR - 1.);

        while (z_prime > REDSHIFT){
            redshifts[counter] = z_prime;

            counter += 1;

            prev_z_prime = z_prime;
            z_prime = ((1.+prev_z_prime) / ZPRIME_STEP_FACTOR - 1.);
        }

        if(USE_LIGHTCONE) {
            // For determining parameters for the interpolation of the boxes for the construction of the light-cone
            // redshift_interpolate_boxes.c required redshifts in increasing order, the array redshifts is in decreasing order. Therefore just invert to leave as much code the same as possible
            z1_LC = redshifts[N_USER_REDSHIFT-1];

            z_LC = start_z = z1_LC;
            slice_ct = 0;
            total_slice_ct = 0;
            num_boxes_interp = 0;
            i = 0;

            while(z1_LC < redshifts[0]) {
                z2_LC = redshifts[N_USER_REDSHIFT-2-i];
                // now do the interpolation
                while (z_LC < z2_LC){ // until we move to the next set of boxes
                    slice_redshifts[total_slice_ct] = z_LC;
                    full_index_LC[total_slice_ct] = total_slice_ct;
                    // check if we filled-up our array and write-out
                    if (slice_ct == HII_DIM){
                        end_z = z_LC;
                        num_boxes_interp += 1;

                        // update quantities
                        slice_ct=0;
                        start_z = end_z;

                    } // we are now continuing with a new interpolation box

                    slice_ct++;
                    total_slice_ct++;
                    z_LC -= dR / drdz(z_LC);
                } // done with this pair of boxes, moving on to the next redshift pair
                redshifts_LC[i] = z1_LC;
                if(i==0) {
                    start_index_LC[i] = 0;
                    end_index_LC[i] = slice_ct; //This value not inclusive
                }
                else {
                    start_index_LC[i] = end_index_LC[i-1]; // Inclusive of this value
                    end_index_LC[i] = slice_ct;
                }
                z1_LC = z2_LC;
                i += 1;
            }

            total_num_boxes = num_boxes_interp;
            remainder_LC = total_slice_ct - num_boxes_interp*HII_DIM;

            final_z = z_LC;

            box_z1 = (float *)calloc(HII_TOT_NUM_PIXELS,sizeof(float));
            box_z2 = (float *)calloc(HII_TOT_NUM_PIXELS,sizeof(float));
            box_interpolate = (float *)calloc(HII_TOT_NUM_PIXELS,sizeof(float));
            box_interpolate_remainder = (float *)calloc((float)HII_DIM*(float)HII_DIM*(float)remainder_LC,sizeof(float));

        }
    }
    Stored_LOS_direction_state_2 = LOS_direction;




    /////////////////   Read the requisite data if the USE_FCOLL_IONISATION_TABLE option is set for the USE_LIGHTCONE option    /////////////////

    /* The python script "Create_ionisaton_fcoll_tables.py" which calls the subsequent C file "Createfcoll_ionisation_LC" are used to create this interpolation table.
     Effectively, this creates an interpolation table for the collapsed fraction for find_HII_bubbles, to avoid having to calculate the full cubic box for each co-eval redshift.
     This can yield > 20 per cent in computational efficiency at the cost of a small decrease in the accuracy of the neutral fraction, and ionisation fraction per voxel. Check this
     accuracy, and whether it is sufficent for your purposes before using this interpolation table. */

    if(USE_LIGHTCONE) {

        if(USE_FCOLL_IONISATION_TABLE) {
            // Note the below code does not yet work for ALPHA != 0
            ReadFcollTable();

        }
    }

    if((R_BUBBLE_MAX > R_MFP_UB)&&(USE_FCOLL_IONISATION_TABLE==1)) {
        printf("\n");
        printf("The interpolation table for the ionisation box collapse fraction does not have the requisite sampling of the R_MFP\n");
        printf("(The selected R_MFP exceeds the upper bound on R_MFP set in the interpolation table)\n");
        printf("(Either reduce R_MFP or re-run the creation of the interpolation table with a sufficiently large upper bound)\n");
        printf("\n");
        printf("Exiting...");
        printf("\n");
        // Probably should do free the memory properly here...
        return 0;
    }

    fprintf(stderr, "\nAllocating memory...\n");
    // Allocate memory for the IGM spin temperature and electron fraction which is stored globally to be taken from Ts.c and used with find_HII_bubbles.c
    Ts_z = (float *)calloc(HII_TOT_NUM_PIXELS,sizeof(float));
    x_e_z = (float *)calloc(HII_TOT_NUM_PIXELS,sizeof(float));

    // Setup an interpolation table for the error function, helpful for calcluating the collapsed fraction (only for the default model, i.e. mass-independent ionising efficiency)
    erfc_arg_min = -15.0;
    erfc_arg_max = 15.0;

    ERFC_NUM_POINTS = 10000;

    ERFC_VALS = calloc(ERFC_NUM_POINTS,sizeof(double));
    ERFC_VALS_DIFF = calloc(ERFC_NUM_POINTS,sizeof(double));

    ArgBinWidth = (erfc_arg_max - erfc_arg_min)/((double)ERFC_NUM_POINTS - 1.);
    InvArgBinWidth = 1./ArgBinWidth;

    for(i=0;i<ERFC_NUM_POINTS;i++) {

        erfc_arg_val = erfc_arg_min + ArgBinWidth*(double)i;

        ERFC_VALS[i] = splined_erfc(erfc_arg_val);
    }

    for(i=0;i<(ERFC_NUM_POINTS-1);i++) {
        ERFC_VALS_DIFF[i] = ERFC_VALS[i+1] - ERFC_VALS[i];
    }



    // Allocate memory for storing the global history of the IGM neutral fraction and brightness temperature contrast
    if(STORE_DATA) {
        aveNF = calloc(N_USER_REDSHIFT,sizeof(double));
        aveTb = calloc(N_USER_REDSHIFT,sizeof(double));
    //JBM:we also save Tkin and Tspin
        aveTkin = calloc(N_USER_REDSHIFT,sizeof(double));
        aveTspin = calloc(N_USER_REDSHIFT,sizeof(double));
        aveXalpha = calloc(N_USER_REDSHIFT,sizeof(double));

        aveTspin_inv = calloc(N_USER_REDSHIFT,sizeof(double));
        aveTkin_inv_sq = calloc(N_USER_REDSHIFT,sizeof(double));
    }   

    // if GenerateNewICs == 1, generate the new initial conditions. This calculates the initial conditions in fourier space, and stores the relevant boxes in memory only (nothing is written to file)
    // At the same time, calculate the density field for calculating the IGM spin temperature.
    // This option must be set if the cosmology is to be varied.
    if(GenerateNewICs) {
        HIRES_density = (float *) fftwf_malloc(sizeof(float)*TOT_FFT_NUM_PIXELS);
        ComputeInitialConditions();
        LOWRES_density_REDSHIFT = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
        LOWRES_velocity_REDSHIFT = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    }

    fprintf(stderr, "\nInitializing MHR: %d\n", INHOMO_RECO);
    if (INHOMO_RECO) {
        
        init_MHR();
    }

    fprintf(stderr, "\nAllocating additional boxes... %d\n", INHOMO_RECO);

    // ALLOCATE AND INITIALIZE ADDITIONAL BOXES NEEDED TO KEEP TRACK OF RECOMBINATIONS (Sobacchi & Mesinger 2014; NEW IN v1.3)
    if (INHOMO_RECO){ //  flag in ANAL_PARAMS.H to determine whether to compute recombinations or not
        z_re = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS); // the redshift at which the cell is ionized
        Gamma12 = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS);  // stores the ionizing backgroud
        N_rec_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS); // cumulative number of recombinations
        N_rec_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

        for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++) {
            z_re[ct] = -1.0;
        }

        // initialize N_rec
        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    *((float *)N_rec_unfiltered + HII_R_FFT_INDEX(i,j,k)) = 0.0;
                }
            }
        }
    } //  end if INHOMO_RECO


    /////////////////   Calculate the filtering scales for all the relevant smoothing scales for the HII_BUBBLES excursion set formalism    /////////////////

    fprintf(stderr, "\nComputing spin temperature fluctuation boxes, %d\n", USE_TS_FLUCT);
    ///////////////////////////////// Decide whether or not the spin temperature fluctuations are to be computed (Ts.c) /////////////////////////////////
    if(USE_TS_FLUCT) {
        ///////////////////////////////// Perform 'Ts.c' /////////////////////////////////

        // This will perform the Ts.c part of the computation, and at the end of each redshift sampling, it will calculate the ionisation and brightness temperature boxes.
        ComputeTsBoxes();
    }
    else {

        // If here, the spin temperature is assumed to be saturated, and thus only the HII_BUBBLES part of the code is performed.
        // User must provide the co-eval redshifts in this case.

        // Reversed the order for iterating to facilitate the code for generating the LC. Here, doesn't matter as assuming Ts >> Tcmb
        for(i=N_USER_REDSHIFT;i--;) {
            // Note, that INHOMO_RECO cannot be set when entering here.
            // INHOMO_RECO must be set with USE_TS_FLUCT
            // This is because no data is stored
            ComputeIonisationBoxes(i,redshifts[i],redshifts[i]+0.2);
        }
    }

    fprintf(stderr, "Storing data: %d\n", STORE_DATA);
    // Storing the global history of the IGM neutral fraction and brightness temperature contrast into a text-file
    if(STORE_DATA) {
        sprintf(filename, "%s/AveData_%f_%f.txt",OUTPUT_FOLDER,INDIVIDUAL_ID,INDIVIDUAL_ID_2);
        if (OVERWRITE_OUTPUT || !file_exists(filename)){
            F=fopen(filename, "wt");
            for(i=0;i<N_USER_REDSHIFT;i++) {
                fprintf(F,"%f\t%e\t%e\n",redshifts[i],aveNF[i],aveTb[i]);
            }
            fclose(F);
        }
    }

    //JBM:also save the spin and gas temperatures
    if(STORE_DATA) {
        sprintf(filename, "%s/AveTemps_%f_%f.txt",OUTPUT_FOLDER,INDIVIDUAL_ID,INDIVIDUAL_ID_2);
        if (OVERWRITE_OUTPUT || !file_exists(filename) ){
            F=fopen(filename, "wt");
            for(i=0;i<N_USER_REDSHIFT;i++) {
                fprintf(F,"%f\t%e\t%e\t%e\t%e\n",redshifts[i],aveTkin[i],aveTspin[i],aveTspin_inv[i], aveXalpha[i]);
            }
            fclose(F);
        }
    }

    // Output the text-file containing the file names of all the 21cm PS calculated from the light-cone boxes
    fprintf(stderr, "Storing 21cm PS filenames: %d\n", USE_LIGHTCONE);
    if(USE_LIGHTCONE) {
        sprintf(filename, "%s/delTps_lightcone_filenames_%f_%f.txt",OUTPUT_FOLDER,INDIVIDUAL_ID,INDIVIDUAL_ID_2);
        if (OVERWRITE_OUTPUT || !file_exists(filename) ){
            fprintf(stderr, "Outputting spin temperature for lightcone to %s", filename);
            F=fopen(filename, "w");
            for(i=0;i<total_num_boxes;i++) {
                fprintf(F,"%s\n",lightcone_box_names[i]);
            }
            fclose(F);
        }
    }


    free_TS_memory();

    // De-allocate all arrays etc. that have been allocated and used
    SAFE_FREE(Ts_z);
    SAFE_FREE(x_e_z);

    SAFE_FREE(PARAM_VALS);

    SAFE_FREE(ERFC_VALS);
    SAFE_FREE(ERFC_VALS_DIFF);

    if(STORE_DATA) {
        SAFE_FREE(aveNF);
        SAFE_FREE(aveTb);

        SAFE_FREE(aveTkin);
        SAFE_FREE(aveTspin_inv);
        SAFE_FREE(aveTkin_inv_sq);
        SAFE_FREE(aveTspin);
        SAFE_FREE(aveXalpha);
    }
    SAFE_FREE(redshifts);

    SAFE_FREE(full_index_LC);
    SAFE_FREE(slice_redshifts);

    if(USE_LIGHTCONE) {
        SAFE_FREE(box_z1);
        SAFE_FREE(box_z2);
        SAFE_FREE(box_interpolate);
        SAFE_FREE(box_interpolate_remainder);
        SAFE_FREE(redshifts_LC);
        SAFE_FREE(start_index_LC);
        SAFE_FREE(end_index_LC);
    }

    free_ps();

    if (INHOMO_RECO) {
        free_MHR();
    }
    SAFE_FFTWF_FREE(N_rec_unfiltered);
    SAFE_FFTWF_FREE(N_rec_filtered);
    SAFE_FFTWF_FREE(z_re);
    SAFE_FFTWF_FREE(Gamma12);

    return 0;
}

void ComputeTsBoxes() {
  printf("COMPUTING TS!!\n");

    /* This is an entire re-write of Ts.c from 21cmFAST. You can refer back to Ts.c in 21cmFAST if this become a little obtuse. The computation has remained the same */

    /////////////////// Defining variables for the computation of Ts.c //////////////
    char filename[500];
    FILE *F, *OUT;
    fftwf_plan plan;

    unsigned long long ct, FCOLL_SHORT_FACTOR;

    int R_ct,i,ii,j,k,i_z,COMPUTE_Ts,x_e_ct,m_xHII_low,m_xHII_high,n_ct, zpp_gridpoint1_int, zpp_gridpoint2_int,zpp_evolve_gridpoint1_int, zpp_evolve_gridpoint2_int;

    short dens_grid_int;

    double Tk_inv_ave, Tk_inv_sq_ave;
    double Tk_ave, J_alpha_ave, xalpha_ave, J_alpha_tot, Xheat_ave, Xion_ave, nuprime, Ts_ave, T21_ave, lower_int_limit,Luminosity_converstion_factor,T_inv_TS_fast_inv;
    double dadia_dzp, dcomp_dzp, dxheat_dt, dxion_source_dt, dxion_sink_dt, T, x_e, dxe_dzp, n_b, dspec_dzp, dxheat_dzp, dxlya_dt, dstarlya_dt, fcoll_R;
    double Trad_fast,xc_fast,xc_inverse,TS_fast,TSold_fast,xa_tilde_fast,TS_prefactor,xa_tilde_prefactor,T_inv,T_inv_sq,xi_power,xa_tilde_fast_arg,Trad_fast_inv,TS_fast_inv,dcomp_dzp_prefactor;

    float growth_factor_z, inverse_growth_factor_z, R, R_factor, zp, mu_for_Ts, filling_factor_of_HI_zp, dzp, prev_zp, zpp, prev_zpp, prev_R, Tk_BC, xe_BC;
    float xHII_call, curr_xalpha, TK, TS, xe, deltax_highz;
    float zpp_for_evolve,dzpp_for_evolve;

    float determine_zpp_max, determine_zpp_min, zpp_grid, zpp_gridpoint1, zpp_gridpoint2,zpp_evolve_gridpoint1, zpp_evolve_gridpoint2, grad1, grad2, grad3, grad4, zpp_bin_width, delNL0_bw_val;
    float OffsetValue, DensityValueLow, min_density, max_density;

    double curr_delNL0, inverse_val,prefactor_1,prefactor_2,dfcoll_dz_val, density_eval1, density_eval2, grid_sigmaTmin, grid_dens_val, dens_grad, dens_width;
    double curr_vbc; //jbm


    float M_MIN_WDM =  M_J_WDM();

    double total_time, total_time2, total_time3, total_time4;

    int Tvir_min_int,Numzp_for_table,counter;
    double X_RAY_Tvir_BinWidth;

    X_RAY_Tvir_BinWidth = (X_RAY_Tvir_UPPERBOUND - X_RAY_Tvir_LOWERBOUND)/( (double)X_RAY_Tvir_POINTS - 1. );

    // Can speed up computation (~20%) by pre-sampling the fcoll field as a function of X_RAY_TVIR_MIN (performed by calling CreateFcollTable.
    // Can be helpful when HII_DIM > ~128, otherwise its easier to just do the full box
    // This table can be created using "CreateFcollTable.c". See this file for further details.
    if(SHORTEN_FCOLL) {

        Tvir_min_int = (int)floor( (log10(X_RAY_Tvir_MIN) - X_RAY_Tvir_LOWERBOUND)/X_RAY_Tvir_BinWidth );

        sprintf(filename, "FcollTvirTable_Numzp_ZPRIME_FACTOR%0.2f_logTvirmin%0.6f_logTvirmax%0.6f_XRAY_POINTS%d_z_end%06.6f_%0.2fMpc_%d.dat",ZPRIME_STEP_FACTOR,X_RAY_Tvir_LOWERBOUND,X_RAY_Tvir_UPPERBOUND,X_RAY_Tvir_POINTS,REDSHIFT,BOX_LEN,HII_DIM);
        F = fopen(filename, "rb");
        fread(&Numzp_for_table, sizeof(int),1,F);
        fclose(F);
    }
    else {
        // Need to take on some number for the memory allocation
        Numzp_for_table = 1;
    }

    // Allocate the memory for this interpolation table
    double ***Fcoll_R_Table = (double ***)calloc(Numzp_for_table,sizeof(double **));
    for(i=0;i<Numzp_for_table;i++) {
        Fcoll_R_Table[i] = (double **)calloc(X_RAY_Tvir_POINTS,sizeof(double *));
        for(j=0;j<X_RAY_Tvir_POINTS;j++) {
            Fcoll_R_Table[i][j] = (double *)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
        }
    }

    if(SHORTEN_FCOLL) {

        sprintf(filename, "FcollTvirTable_ZPRIME_FACTOR%0.2f_logTvirmin%0.6f_logTvirmax%0.6f_XRAY_POINTS%d_z_end%06.6f_%0.2fMpc_%d.dat",ZPRIME_STEP_FACTOR,X_RAY_Tvir_LOWERBOUND,X_RAY_Tvir_UPPERBOUND,X_RAY_Tvir_POINTS,REDSHIFT,BOX_LEN,HII_DIM);
        F = fopen(filename, "rb");
        for(i=0;i<Numzp_for_table;i++) {
            for(j=0;j<X_RAY_Tvir_POINTS;j++) {
                fread(Fcoll_R_Table[i][j], sizeof(double),NUM_FILTER_STEPS_FOR_Ts,F);
            }
        }
    }


    double factor_relative_velocities; //JBM. 1 if not using rel vel, and F_coll (v) if you are.

    double fcoll_cell;
    double dfcoll_cell;
    double dlogPS_v;
    double dlog_ST_PS;


    // Initialise arrays to be used for the Ts.c computation //
    init_21cmMC_Ts_arrays();
    fprintf(stderr, "Ts arrays initialized\n");



    ///////////////////////////////  BEGIN INITIALIZATION   //////////////////////////////
    growth_factor_z = dicke(REDSHIFT);
    inverse_growth_factor_z = 1./growth_factor_z;

    if (X_RAY_Tvir_MIN < 9.99999e3) // neutral IGM
        mu_for_Ts = 1.22;
    else // ionized IGM
        mu_for_Ts = 0.6;

    //set the minimum ionizing source mass
    M_MIN_at_z = get_M_min_ion(REDSHIFT);

    // Initialize some interpolation tables
    init_heat();
    fprintf(stderr, "Finished init_heat()\n");

    //JBM:initialize the vcb effect: read the interpolation tables
    char ending_filename [10];
    int len = sprintf(ending_filename,"_%d.dat",FEEDBACK_LEVEL);


    char filename_fcollapse [200];
    strcpy (filename_fcollapse,FILENAME_FCOLL_VCB);
    strcat (filename_fcollapse, ending_filename);


    F = fopen(filename_fcollapse,"r");
    double trash;

    for(i=0;i<NZINT;i++) {
        for(j=0;j<NVINT;j++) {
            fscanf(F, "%le %le %le" , &trash, &trash, &logFcoll_vcb[i][j]); //it has z, v, Fcoll, but we don't save z or v.
            // printf("%d %d %le \n",i,j,exp(logFcoll_vcb[i][j]));
        }
    }

    //JBM:and the sigma (M_cool(v,z),z,v) too: read the interpolation table
    char filename_sigmacool [200];
    strcpy (filename_sigmacool,FILENAME_SIGMACOOL_VCB);
    strcat (filename_sigmacool, ending_filename);
    F = fopen(filename_sigmacool,"r");

    for(i=0;i<NZINT;i++) {
        for(j=0;j<NVINT;j++) {
            fscanf(F, "%le %le %le" , &trash, &trash, &sigmacool_vcb[i][j]); //it has z, v, Fcoll, but we don't save z or v.
            // printf("%d %d %le \n",i,j,exp(logFcoll_vcb[i][j]));
        }
    }

    // check if we are in the really high z regime before the first stars; if so, simple
    if (REDSHIFT > Z_HEAT_MAX){

        /**** NOTE: THIS NEEDS TO CHANGE. Though it'll only cause problems if this condition is met (which should never happen) ****/

        //(FgtrM(REDSHIFT, FMAX(TtoM(REDSHIFT, X_RAY_Tvir_MIN, mu_for_Ts),  M_MIN_WDM)) < 1e-15 ){
        xe = xion_RECFAST(REDSHIFT,0);
        TK = T_RECFAST(REDSHIFT,0);

        // open input
        sprintf(filename, "%s/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc",BOXES_INPUT_FOLDER, REDSHIFT, HII_DIM, BOX_LEN);
        F = fopen(filename, "rb");

        // open output
        sprintf(filename, "%s/Ts_z%06.2f_zetaX%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc",BOXES_INPUT_FOLDER, REDSHIFT, HII_EFF_FACTOR, X_RAY_SPEC_INDEX, X_RAY_Tvir_MIN, R_BUBBLE_MAX, Pop, HII_DIM, BOX_LEN);

        // read file
        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    fread(&deltax_highz, sizeof(float), 1, F);

                    // compute the spin temperature
                    TS = get_Ts(REDSHIFT, deltax_highz, TK, xe, 0, &curr_xalpha);

                    // and print it out
                    fwrite(&TS, sizeof(float), 1, OUT);
                }
            }
        }

        destruct_heat(); fclose(F); fclose(OUT);
    }
    else {

        // set boundary conditions for the evolution equations->  values of Tk and x_e at Z_HEAT_MAX
        if (XION_at_Z_HEAT_MAX > 0) // user has opted to use his/her own value
            xe_BC = XION_at_Z_HEAT_MAX;
        else// will use the results obtained from recfast
            xe_BC = xion_RECFAST(Z_HEAT_MAX,0);
        if (TK_at_Z_HEAT_MAX > 0)
            Tk_BC = TK_at_Z_HEAT_MAX;
        else
            Tk_BC = T_RECFAST(Z_HEAT_MAX,0);

        /////////////// Create the z=0 non-linear density fields smoothed on scale R to be used in computing fcoll //////////////
        R = L_FACTOR*BOX_LEN/(float)HII_DIM;
        R_factor = pow(R_XLy_MAX/R, 1/(float)NUM_FILTER_STEPS_FOR_Ts);
        // R_factor = pow(E, log(HII_DIM)/(float)NUM_FILTER_STEPS_FOR_Ts);


        ///////////////////  Read in density box at z-prime  ///////////////
        if(GenerateNewICs) {

            // If GenerateNewICs == 1, we are generating a new set of initial conditions and density field. Hence, calculate the density field to be used for Ts.c

            ComputePerturbField(REDSHIFT);

            for (i=0; i<HII_DIM; i++){
                for (j=0; j<HII_DIM; j++){
                    for (k=0; k<HII_DIM; k++){
                        *((float *)unfiltered_box + HII_R_FFT_INDEX(i,j,k)) = LOWRES_density_REDSHIFT[HII_R_INDEX(i,j,k)];
                    }
                }
            }
        }
        else {

            // Read in a pre-computed density field which is stored in the "Boxes" folder

            // allocate memory for the nonlinear density field and open file
            sprintf(filename, "%s/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc",BOXES_INPUT_FOLDER, REDSHIFT, HII_DIM, BOX_LEN);
            F = fopen(filename, "rb");
            for (i=0; i<HII_DIM; i++){
                for (j=0; j<HII_DIM; j++){
                    for (k=0; k<HII_DIM; k++){
                        fread((float *)unfiltered_box + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F);
                    }
                }
            }
            fclose(F);
        }


        //JBM:we also read the velocity field at z_dec, also in the /Boxes folder.
        fprintf(stderr, "Use relative velocities: %d\n", USE_RELATIVE_VELOCITIES);
        if(USE_RELATIVE_VELOCITIES){
        sprintf(filename, "%s/smoothed_vcb_x_z0.00_%i_%.0fMpc",BOXES_INPUT_FOLDER, HII_DIM, BOX_LEN);
            F = fopen(filename, "rb");
            for (i=0; i<HII_DIM; i++){
                for (j=0; j<HII_DIM; j++){
                    for (k=0; k<HII_DIM; k++){
                        fread((float *)unfiltered_vcb_box + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F);
                    }
                }
            }
            fclose(F);
        }


        fprintf(stderr, "Transform unfiltered box to k-space to prepare for filtering\n");
        ////////////////// Transform unfiltered box to k-space to prepare for filtering /////////////////
        plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)unfiltered_box, (fftwf_complex *)unfiltered_box, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
        fftwf_cleanup();

	    //JBM:transform the vcb box as well
        plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)unfiltered_vcb_box, (fftwf_complex *)unfiltered_vcb_box, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
        fftwf_cleanup();


        // remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
        // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
        for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
          unfiltered_box[ct] /= (float)HII_TOT_NUM_PIXELS;
          unfiltered_vcb_box[ct] /= (float)HII_TOT_NUM_PIXELS;
        }


        fprintf(stderr, "Smooth the density field\n");
        // Smooth the density field, at the same time store the minimum and maximum densities for their usage in the interpolation tables
        for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){

            R_values[R_ct] = R;
            sigma_atR[R_ct] = sigma_z0(RtoM(R));

            // copy over unfiltered box
            memcpy(box, unfiltered_box, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
            memcpy(box_vcb, unfiltered_vcb_box, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

            if (R_ct > 0){ // don't filter on cell size
                HII_filter(box, HEAT_FILTER, R);
                HII_filter(box_vcb, HEAT_FILTER, R);
            }
            // now fft back to real space
            plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)box, (float *)box, FFTW_ESTIMATE);
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);
            fftwf_cleanup();

            plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)box_vcb, (float *)box_vcb, FFTW_ESTIMATE);
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);
            fftwf_cleanup();

            min_density = 0.0;
            max_density = 0.0;

            // copy over the values
            for (i=HII_DIM; i--;){
                for (j=HII_DIM; j--;){
                    for (k=HII_DIM; k--;){
                        curr_delNL0 = *((float *) box + HII_R_FFT_INDEX(i,j,k));

                        if (curr_delNL0 < -1){ // correct for alliasing in the filtering step
                            curr_delNL0 = -1+FRACT_FLOAT_ERR;
                        }

                        // and linearly extrapolate to z=0
                        curr_delNL0 *= inverse_growth_factor_z;

                        delNL0_rev[HII_R_INDEX(i,j,k)][R_ct] = curr_delNL0;

                        if(curr_delNL0 < min_density) {
                            min_density = curr_delNL0;
                        }
                        if(curr_delNL0 > max_density) {
                            max_density = curr_delNL0;
                        }

                        //JBM:we also save the velocity:
                        curr_vbc = *((float *) box_vcb + HII_R_FFT_INDEX(i,j,k));
                        if(curr_vbc < VBCMIN) {
                            curr_vbc = VBCMIN;
                        }
                        if(curr_vbc > VBCMAX) {
                            curr_vbc = VBCMAX;
                        }
                        //JBM:we cap the velocity to avoid unusually large values for simplicity.

                        vcb_rev[HII_R_INDEX(i,j,k)][R_ct] = curr_vbc;

                        if(USE_RELATIVE_VELOCITIES && ZERO_REL_VELOCITY){
                          vcb_rev[HII_R_INDEX(i,j,k)][R_ct] = 0.0;
                        }
                        if(USE_RELATIVE_VELOCITIES && AVG_REL_VELOCITY){
                          vcb_rev[HII_R_INDEX(i,j,k)][R_ct] = VCB_AVG;
                        }
                    }
                }
            }
            if(min_density < 0.0) {
                delNL0_LL[R_ct] = min_density*1.001;
                delNL0_Offset[R_ct] = 1.e-6 - (delNL0_LL[R_ct]);
            }
            else {
                delNL0_LL[R_ct] = min_density*0.999;
                delNL0_Offset[R_ct] = 1.e-6 + (delNL0_LL[R_ct]);
            }
            if(max_density < 0.0) {
                delNL0_UL[R_ct] = max_density*0.999;
            }
            else {
                delNL0_UL[R_ct] = max_density*1.001;
            }

            R *= R_factor;

        } //end for loop through the filter scales R

        //moved above since we have 2 boxes now.
        // fftwf_destroy_plan(plan);
        // fftwf_cleanup();

        // and initialize to the boundary values at Z_HEAT_END
        for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
            Tk_box[ct] = Tk_BC;
            x_e_box[ct] = xe_BC;
        }
        x_e_ave = xe_BC;
        Tk_ave = Tk_BC;


        ////////////////////////////    END INITIALIZATION   /////////////////////////////

        fprintf(stderr, "Main integral\n");
        // main trapezoidal integral over z' (see eq. ? in Mesinger et al. 2009)
        zp = REDSHIFT*1.0001; //higher for rounding
        while (zp < Z_HEAT_MAX)
            zp = ((1+zp)*ZPRIME_STEP_FACTOR - 1);
        prev_zp = Z_HEAT_MAX;
        zp = ((1+zp)/ ZPRIME_STEP_FACTOR - 1);
        dzp = zp - prev_zp;
        COMPUTE_Ts = 0;

        determine_zpp_min = REDSHIFT*0.999;

        for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
            if (R_ct==0){
                prev_zpp = zp;
                prev_R = 0;
            }
            else{
                prev_zpp = zpp_edge[R_ct-1];
                prev_R = R_values[R_ct-1];
            }
            zpp_edge[R_ct] = prev_zpp - (R_values[R_ct] - prev_R)*CMperMPC / drdz(prev_zpp); // cell size
            zpp = (zpp_edge[R_ct]+prev_zpp)*0.5; // average redshift value of shell: z'' + 0.5 * dz''
        }

        determine_zpp_max = zpp*1.001;


        ////////////////////////////    Create and fill interpolation tables to be used by Ts.c   /////////////////////////////


        fprintf(stderr, "Initiating interpolation tables\n");
        // An interpolation table for f_coll (delta vs redshift)
        init_FcollTable(determine_zpp_min,determine_zpp_max);

        zpp_bin_width = (determine_zpp_max - determine_zpp_min)/((float)zpp_interp_points-1.0);

        dens_width = 1./((double)dens_Ninterp - 1.);

        // Determine the sampling of the density values, for the various interpolation tables
        for(ii=0;ii<NUM_FILTER_STEPS_FOR_Ts;ii++) {
            log10delNL0_diff_UL[ii] = log10( delNL0_UL[ii] + delNL0_Offset[ii] );
            log10delNL0_diff[ii] = log10( delNL0_LL[ii] + delNL0_Offset[ii] );
            delNL0_bw[ii] = ( log10delNL0_diff_UL[ii] - log10delNL0_diff[ii] )*dens_width;
            delNL0_ibw[ii] = 1./delNL0_bw[ii];
        }



        // Gridding the density values for the interpolation tables
        for(ii=0;ii<NUM_FILTER_STEPS_FOR_Ts;ii++) {
            for(j=0;j<dens_Ninterp;j++) {
                grid_dens[ii][j] = log10delNL0_diff[ii] + ( log10delNL0_diff_UL[ii] - log10delNL0_diff[ii] )*dens_width*(double)j;
                grid_dens[ii][j] = pow(10,grid_dens[ii][j]) - delNL0_Offset[ii];
            }
        }


        // Calculate the sigma_z and Fgtr_M values for each point in the interpolation table
        for(i=0;i<zpp_interp_points;i++) {
            zpp_grid = determine_zpp_min + (determine_zpp_max - determine_zpp_min)*(float)i/((float)zpp_interp_points-1.0);

            Sigma_Tmin_grid[i] = sigma_z0(FMAX(TtoM(zpp_grid, X_RAY_Tvir_MIN, mu_for_Ts),  M_MIN_WDM));
            ST_over_PS_arg_grid[i] = FgtrM_st(zpp_grid, FMAX(TtoM(zpp_grid, X_RAY_Tvir_MIN, mu_for_Ts),  M_MIN_WDM));

        }

        // Create the interpolation tables for the derivative of the collapsed fraction and the collapse fraction itself
          for(ii=0;ii<NUM_FILTER_STEPS_FOR_Ts;ii++) {
              for(i=0;i<zpp_interp_points;i++) {

                  zpp_grid = determine_zpp_min + (determine_zpp_max - determine_zpp_min)*(float)i/((float)zpp_interp_points-1.0);
                  grid_sigmaTmin = Sigma_Tmin_grid[i];

                  if(USE_RELATIVE_VELOCITIES){
                    grid_sigmaTmin = interpol_linear_2D(ZINT_MIN, ZINT_STEP, NZINT, VINT_MIN, VINT_STEP, NVINT,
                          	   sigmacool_vcb, zpp_grid, VCB_AVG*(1.0-ZERO_REL_VELOCITY))/dicke(zpp_grid);
                  }
                  //JBM: we take the \sigma(M) of molecular cooling at v=vrms for Fcoll(avg), since <F(v)>~F(<v>) to good approximation.


                  for(j=0;j<dens_Ninterp;j++) {

                      grid_dens_val = grid_dens[ii][j];
                      if(!SHORTEN_FCOLL) {
                          fcoll_R_grid[ii][i][j] = sigmaparam_FgtrM_bias(zpp_grid, grid_sigmaTmin, grid_dens_val, sigma_atR[ii]);
                      }
                      dfcoll_dz_grid[ii][i][j] = dfcoll_dz(zpp_grid, grid_sigmaTmin, grid_dens_val, sigma_atR[ii]);
                  }
              }
          }



        fprintf(stderr, "Determining and evaluating grid point locations for solving interpolation tables\n");
        // Determine the grid point locations for solving the interpolation tables
        for (box_ct=HII_TOT_NUM_PIXELS; box_ct--;){
            for (R_ct=NUM_FILTER_STEPS_FOR_Ts; R_ct--;){
                SingleVal_int[R_ct] = (short)floor( ( log10(delNL0_rev[box_ct][R_ct] + delNL0_Offset[R_ct]) - log10delNL0_diff[R_ct] )*delNL0_ibw[R_ct]);
            }
            memcpy(dens_grid_int_vals[box_ct],SingleVal_int,sizeof(short)*NUM_FILTER_STEPS_FOR_Ts);
        }

        // Evaluating the interpolated density field points (for using the interpolation tables for fcoll and dfcoll_dz)
        for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
            OffsetValue = delNL0_Offset[R_ct];
            DensityValueLow = delNL0_LL[R_ct];
            delNL0_bw_val = delNL0_bw[R_ct];

            for(i=0;i<dens_Ninterp;i++) {
                density_gridpoints[i][R_ct] = pow(10.,( log10( DensityValueLow + OffsetValue) + delNL0_bw_val*((float)i) )) - OffsetValue;
            }
        }

        counter = 0;

        fprintf(stderr, "Starting calculation for IGM spin temperature\n");
        // This is the main loop for calculating the IGM spin temperature. Structure drastically different from Ts.c in 21cmFAST, however algorithm and computation remain the same.
        while (zp > REDSHIFT){
            // check if we will next compute the spin temperature (i.e. if this is the final zp step)
            if (Ts_verbose || (((1+zp) / ZPRIME_STEP_FACTOR) < (REDSHIFT+1)) )
                COMPUTE_Ts = 1;

            // check if we are in the really high z regime before the first stars..
            if (FgtrM(zp, FMAX(TtoM(zp, X_RAY_Tvir_MIN, mu_for_Ts),  M_MIN_WDM)) < 1e-15 )
                NO_LIGHT = 1;
            else
                NO_LIGHT = 0;

            M_MIN_at_zp = get_M_min_ion(zp);
            filling_factor_of_HI_zp = 1 - HII_EFF_FACTOR * FgtrM_st(zp, M_MIN_at_zp) / (1.0 - x_e_ave);
            if (filling_factor_of_HI_zp > 1) filling_factor_of_HI_zp=1;

            // let's initialize an array of redshifts (z'') corresponding to the
            // far edge of the dz'' filtering shells
            // and the corresponding minimum halo scale, sigma_Tmin,
            // as well as an array of the frequency integrals
            for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
                if (R_ct==0){
                    prev_zpp = zp;
                    prev_R = 0;
                }
                else{
                    prev_zpp = zpp_edge[R_ct-1];
                    prev_R = R_values[R_ct-1];
                }

                zpp_edge[R_ct] = prev_zpp - (R_values[R_ct] - prev_R)*CMperMPC / drdz(prev_zpp); // cell size
                zpp = (zpp_edge[R_ct]+prev_zpp)*0.5; // average redshift value of shell: z'' + 0.5 * dz''

                // Determining values for the evaluating the interpolation table
                zpp_gridpoint1_int = (int)floor((zpp - determine_zpp_min)/zpp_bin_width);
                zpp_gridpoint2_int = zpp_gridpoint1_int + 1;

                zpp_gridpoint1 = determine_zpp_min + zpp_bin_width*(float)zpp_gridpoint1_int;
                zpp_gridpoint2 = determine_zpp_min + zpp_bin_width*(float)zpp_gridpoint2_int;

                grad1 = ( zpp_gridpoint2 - zpp )/( zpp_gridpoint2 - zpp_gridpoint1 );
                grad2 = ( zpp - zpp_gridpoint1 )/( zpp_gridpoint2 - zpp_gridpoint1 );

                sigma_Tmin[R_ct] = Sigma_Tmin_grid[zpp_gridpoint1_int] + grad2*( Sigma_Tmin_grid[zpp_gridpoint2_int] - Sigma_Tmin_grid[zpp_gridpoint1_int] );

                // let's now normalize the total collapse fraction so that the mean is the
                // Sheth-Torman collapse fraction

                zpp_for_evolve_list[R_ct] = zpp;
                if (R_ct==0){
                    dzpp_for_evolve = zp - zpp_edge[0];
                }
                else{
                    dzpp_for_evolve = zpp_edge[R_ct-1] - zpp_edge[R_ct];
                }
                zpp_growth[R_ct] = dicke(zpp);

                  // Evaluating the interpolation table for the collapse fraction and its derivative
                  for(i=0;i<(dens_Ninterp-1);i++) {
                      dens_grad = 1./( density_gridpoints[i+1][R_ct] - density_gridpoints[i][R_ct] );

                      if(!SHORTEN_FCOLL) {
                          fcoll_interp1[i][R_ct] = ( ( fcoll_R_grid[R_ct][zpp_gridpoint1_int][i] )*grad1 + ( fcoll_R_grid[R_ct][zpp_gridpoint2_int][i] )*grad2 )*dens_grad;
                          fcoll_interp2[i][R_ct] = ( ( fcoll_R_grid[R_ct][zpp_gridpoint1_int][i+1] )*grad1 + ( fcoll_R_grid[R_ct][zpp_gridpoint2_int][i+1] )*grad2 )*dens_grad;
                      }

                      dfcoll_interp1[i][R_ct] = ( ( dfcoll_dz_grid[R_ct][zpp_gridpoint1_int][i] )*grad1 + ( dfcoll_dz_grid[R_ct][zpp_gridpoint2_int][i] )*grad2 )*dens_grad;
                      dfcoll_interp2[i][R_ct] = ( ( dfcoll_dz_grid[R_ct][zpp_gridpoint1_int][i+1] )*grad1 + ( dfcoll_dz_grid[R_ct][zpp_gridpoint2_int][i+1] )*grad2 )*dens_grad;
                  }


                fcoll_R_array[R_ct] = 0.0;

                // Using the interpolated values to update arrays of relevant quanties for the IGM spin temperature calculation
                ST_over_PS[R_ct] = dzpp_for_evolve * pow(1+zpp, -X_RAY_SPEC_INDEX);

                if(!(USE_RELATIVE_VELOCITIES) || zpp<ZLYMANWERNER){ //JBM: if we use vcb do not add ST like this.
                  ST_over_PS[R_ct] *= ( ST_over_PS_arg_grid[zpp_gridpoint1_int] + grad2*( ST_over_PS_arg_grid[zpp_gridpoint2_int] - ST_over_PS_arg_grid[zpp_gridpoint1_int] ) );
                }



                lower_int_limit = FMAX(nu_tau_one_approx(zp, zpp, x_e_ave, filling_factor_of_HI_zp), NU_X_THRESH);

                if (filling_factor_of_HI_zp < 0) filling_factor_of_HI_zp = 0; // for global evol; nu_tau_one above treats negative (post_reionization) inferred filling factors properly

                // set up frequency integral table for later interpolation for the cell's x_e value
                for (x_e_ct = 0; x_e_ct < x_int_NXHII; x_e_ct++){
                    freq_int_heat_tbl[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit, 0);
                    freq_int_ion_tbl[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit, 1);
                    if (COMPUTE_Ts)
                        freq_int_lya_tbl[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit, 2);
                }



                // and create the sum over Lya transitions from direct Lyn flux
                sum_lyn[R_ct] = 0;
                for (n_ct=NSPEC_MAX; n_ct>=2; n_ct--){
                    if (zpp > zmax(zp, n_ct))
                        continue;

                    nuprime = nu_n(n_ct)*(1+zpp)/(1.0+zp);
                    sum_lyn[R_ct] += frecycle(n_ct) * spectral_emissivity(nuprime, 0);
                }
            } // end loop over R_ct filter steps

            // Calculate fcoll for each smoothing radius

            // Can speed up computation (~20%) by pre-sampling the fcoll field as a function of X_RAY_TVIR_MIN (performed by calling CreateFcollTable.
            // Can be helpful when HII_DIM > ~128, otherwise its easier to just do the full box
            if(SHORTEN_FCOLL) {
                for(R_ct=0;R_ct<NUM_FILTER_STEPS_FOR_Ts;R_ct++) {
                    fcoll_R = Fcoll_R_Table[counter][Tvir_min_int][R_ct] + ( log10(X_RAY_Tvir_MIN) - ( X_RAY_Tvir_LOWERBOUND + (double)Tvir_min_int*X_RAY_Tvir_BinWidth ) )*( Fcoll_R_Table[counter][Tvir_min_int+1][R_ct] - Fcoll_R_Table[counter][Tvir_min_int][R_ct] )/X_RAY_Tvir_BinWidth;

                    ST_over_PS[R_ct] = ST_over_PS[R_ct]/fcoll_R;
                }
            }
            else {
                for (box_ct=HII_TOT_NUM_PIXELS; box_ct--;){
                    for (R_ct=NUM_FILTER_STEPS_FOR_Ts; R_ct--;){
                        fcoll_R_array[R_ct] += ( fcoll_interp1[dens_grid_int_vals[box_ct][R_ct]][R_ct]*( density_gridpoints[dens_grid_int_vals[box_ct][R_ct] + 1][R_ct] - delNL0_rev[box_ct][R_ct] ) + fcoll_interp2[dens_grid_int_vals[box_ct][R_ct]][R_ct]*( delNL0_rev[box_ct][R_ct] - density_gridpoints[dens_grid_int_vals[box_ct][R_ct]][R_ct] ) );
                    }
                }
                for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
                    ST_over_PS[R_ct] = ST_over_PS[R_ct]/(fcoll_R_array[R_ct]/(double)HII_TOT_NUM_PIXELS);
                }
            }
            // scroll through each cell and update the temperature and residual ionization fraction
            growth_factor_zp = dicke(zp);
            dgrowth_factor_dzp = ddicke_dz(zp);
            dt_dzp = dtdz(zp);



            // Conversion of the input bolometric luminosity to a ZETA_X, as used to be used in Ts.c
            // Conversion here means the code otherwise remains the same as the original Ts.c
            if(fabs(X_RAY_SPEC_INDEX - 1.0) < 0.000001) {
                Luminosity_converstion_factor = NU_X_THRESH * log( NU_X_BAND_MAX/NU_X_THRESH );
                Luminosity_converstion_factor = 1./Luminosity_converstion_factor;
            }
            else {
                Luminosity_converstion_factor = pow( NU_X_BAND_MAX , 1. - X_RAY_SPEC_INDEX ) - pow( NU_X_THRESH , 1. - X_RAY_SPEC_INDEX ) ;
                Luminosity_converstion_factor = 1./Luminosity_converstion_factor;
                Luminosity_converstion_factor *= pow( NU_X_THRESH, - X_RAY_SPEC_INDEX )*(1 - X_RAY_SPEC_INDEX);
            }
            // Finally, convert to the correct units. NU_over_EV*hplank as only want to divide by eV -> erg (owing to the definition of Luminosity)
            Luminosity_converstion_factor *= (3.1556226e7)/(hplank);

            // Leave the original 21cmFAST code for reference. Refer to Greig & Mesinger (2017) for the new parameterisation.
            const_zp_prefactor = ( L_X * Luminosity_converstion_factor ) / NU_X_THRESH * C * F_STAR * OMb * RHOcrit * pow(CMperMPC, -3) * pow(1+zp, X_RAY_SPEC_INDEX+3);
            // This line below is kept purely for reference w.r.t to the original 21cmFAST
            // const_zp_prefactor = ZETA_X * X_RAY_SPEC_INDEX / NU_X_THRESH * C * F_STAR * OMb * RHOcrit * pow(CMperMPC, -3) * pow(1+zp, X_RAY_SPEC_INDEX+3);

            //////////////////////////////  LOOP THROUGH BOX //////////////////////////////

            J_alpha_ave = xalpha_ave = Xheat_ave = Xion_ave = 0.;

            // Extra pre-factors etc. are defined here, as they are independent of the density field, and only have to be computed once per z' or R_ct, rather than each box_ct
            for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
                dstarlya_dt_prefactor[R_ct]  = ( pow(1+zp,2)*(1+zpp_for_evolve_list[R_ct]) * sum_lyn[R_ct] )/( pow(1+zpp_for_evolve_list[R_ct], -X_RAY_SPEC_INDEX) );
            }

            // Required quantities for calculating the IGM spin temperature
            // Note: These used to be determined in evolveInt (and other functions). But I moved them all here, into a single location.
            Trad_fast = T_cmb*(1.0+zp);
            Trad_fast_inv = 1.0/Trad_fast;
            TS_prefactor = pow(1.0e-7*(1.342881e-7 / hubble(zp))*No*pow(1+zp,3),1./3.);
            xa_tilde_prefactor = 1.66e11/(1.0+zp);

            xc_inverse =  pow(1.0+zp,3.0)*T21/( Trad_fast*A10_HYPERFINE );

            dcomp_dzp_prefactor = (-1.51e-4)/(hubble(zp)/Ho)/hlittle*pow(Trad_fast,4.0)/(1.0+zp);

            prefactor_1 = N_b0 * pow(1+zp, 3);
            prefactor_2 = F_STAR * C * N_b0 / FOURPI;

            x_e_ave = 0; Tk_ave = 0; Ts_ave = 0;
            T21_ave=0; Tk_inv_ave = 0; Tk_inv_sq_ave = 0;//JBM also created that just to print it.

            // Note: I have removed the call to evolveInt, as is default in the original Ts.c. Removal of evolveInt and moving that computation below, removes unneccesary repeated computations
            // and allows for the interpolation tables that are now used to be more easily computed

            // Can precompute these quantities, independent of the density field (i.e. box_ct)
            for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
                for (i=0; i<(x_int_NXHII-1); i++) {
                    m_xHII_low = i;
                    m_xHII_high = m_xHII_low + 1;

                    inverse_diff[i] = 1./(x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
                    freq_int_heat_tbl_diff[i][R_ct] = freq_int_heat_tbl[m_xHII_high][R_ct] - freq_int_heat_tbl[m_xHII_low][R_ct];
                    freq_int_ion_tbl_diff[i][R_ct] = freq_int_ion_tbl[m_xHII_high][R_ct] - freq_int_ion_tbl[m_xHII_low][R_ct];
                    freq_int_lya_tbl_diff[i][R_ct] = freq_int_lya_tbl[m_xHII_high][R_ct] - freq_int_lya_tbl[m_xHII_low][R_ct];
                }
            }
            // Main loop over the entire box for the IGM spin temperature and relevant quantities.


            for (box_ct=HII_TOT_NUM_PIXELS; box_ct--;){
                if (!COMPUTE_Ts && (Tk_box[box_ct] > MAX_TK)) //just leave it alone and go to next value
                    continue;

                x_e = x_e_box[box_ct];
                T = Tk_box[box_ct];

                xHII_call = x_e;

                // Check if ionized fraction is within boundaries; if not, adjust to be within
                if (xHII_call > x_int_XHII[x_int_NXHII-1]*0.999) {
                    xHII_call = x_int_XHII[x_int_NXHII-1]*0.999;
                } else if (xHII_call < x_int_XHII[0]) {
                    xHII_call = 1.001*x_int_XHII[0];
                }
                //interpolate to correct nu integral value based on the cell's ionization state

                m_xHII_low = locate_xHII_index(xHII_call);

                inverse_val = (xHII_call - x_int_XHII[m_xHII_low])*inverse_diff[m_xHII_low];

                // First, let's do the trapEzoidal integration over zpp
                dxheat_dt = 0;
                dxion_source_dt = 0;
                dxlya_dt = 0;
                dstarlya_dt = 0;

                curr_delNL0 = delNL0_rev[box_ct][0];



                if (!NO_LIGHT){
                    // Now determine all the differentials for the heating/ionisation rate equations
                    for (R_ct=NUM_FILTER_STEPS_FOR_Ts; R_ct--;){

                      double zpp_curr=zpp_for_evolve_list[R_ct];


                      //JBM:the velocity part, we add the ST as a function of velocity here.
                      factor_relative_velocities = 1.0;
                      if(zpp_curr>=ZLYMANWERNER && USE_RELATIVE_VELOCITIES){
                         factor_relative_velocities = fmax(FCOLLMIN, exp( interpol_linear_2D(ZINT_MIN, ZINT_STEP, NZINT, VINT_MIN, VINT_STEP, NVINT,
                          	   logFcoll_vcb, zpp_curr, vcb_rev[box_ct][R_ct]) ));

                         double currentR = R_values[R_ct];
                         double currentgrowth = zpp_growth[R_ct];
                         double sigma_large = sigma_atR[R_ct]; //at z=0.
                         double sigmacool_curr = interpol_linear_2D(ZINT_MIN, ZINT_STEP, NZINT, VINT_MIN, VINT_STEP, NVINT,
                                    sigmacool_vcb, zpp_curr, vcb_rev[box_ct][R_ct])/currentgrowth;

                         double argument = (Deltac/currentgrowth - delNL0_rev[box_ct][R_ct])/
                                   sqrt(2.0 * (sigmacool_curr*sigmacool_curr - sigma_large*sigma_large));

                        double epsilon_derivative = 0.01; //relative change in z to find derivative.

                         double dcurrentgrowth_dz = (dicke(zpp_curr*(1.0+epsilon_derivative)) - dicke(zpp_curr*(1.0-epsilon_derivative)))/
                                                    (2.0*zpp_curr*epsilon_derivative);

                         double fcolldot_EPS = -2.0/sqrt(PI) * exp(-argument*argument) * (-Deltac*dcurrentgrowth_dz)/(currentgrowth*currentgrowth)
                                          /sqrt(2.0 * (sigmacool_curr*sigmacool_curr - sigma_large*sigma_large));


                          dfcoll_dz_val = factor_relative_velocities*ST_over_PS[R_ct]*(1.+delNL0_rev[box_ct][R_ct]*zpp_growth[R_ct])*
                                            fcolldot_EPS;


                      }
                      else{//standard 21cmmc
                          dfcoll_dz_val = ST_over_PS[R_ct]*(1.+delNL0_rev[box_ct][R_ct]*zpp_growth[R_ct])*
                          	( dfcoll_interp1[dens_grid_int_vals[box_ct][R_ct]][R_ct]
                            *(density_gridpoints[dens_grid_int_vals[box_ct][R_ct] + 1][R_ct] - delNL0_rev[box_ct][R_ct]) +
                            dfcoll_interp2[dens_grid_int_vals[box_ct][R_ct]][R_ct]
                            *(delNL0_rev[box_ct][R_ct] - density_gridpoints[dens_grid_int_vals[box_ct][R_ct]][R_ct]) );
                      }

                        dxheat_dt += dfcoll_dz_val * ( (freq_int_heat_tbl_diff[m_xHII_low][R_ct])*inverse_val + freq_int_heat_tbl[m_xHII_low][R_ct] );
                        dxion_source_dt += dfcoll_dz_val * ( (freq_int_ion_tbl_diff[m_xHII_low][R_ct])*inverse_val + freq_int_ion_tbl[m_xHII_low][R_ct] );

                        if (COMPUTE_Ts){
                            dxlya_dt += dfcoll_dz_val * ( (freq_int_lya_tbl_diff[m_xHII_low][R_ct])*inverse_val + freq_int_lya_tbl[m_xHII_low][R_ct] );
                            dstarlya_dt += dfcoll_dz_val*dstarlya_dt_prefactor[R_ct];
                        }
                    }
                }

                // add prefactors
                dxheat_dt *= const_zp_prefactor;
                dxion_source_dt *= const_zp_prefactor;
                if (COMPUTE_Ts){
                    dxlya_dt *= const_zp_prefactor*prefactor_1 * (1.+curr_delNL0*growth_factor_zp);
                    dstarlya_dt *= prefactor_2;
                }

                // Now we can solve the evolution equations  //

                // First let's do dxe_dzp //
                dxion_sink_dt = alpha_A(T) * CLUMPING_FACTOR * x_e*x_e * f_H * prefactor_1 * (1.+curr_delNL0*growth_factor_zp);
                dxe_dzp = dt_dzp*(dxion_source_dt - dxion_sink_dt );

                // Next, let's get the temperature components //
                // first, adiabatic term
                dadia_dzp = 3/(1.0+zp);
                if (fabs(curr_delNL0) > FRACT_FLOAT_ERR) // add adiabatic heating/cooling from structure formation
                    dadia_dzp += dgrowth_factor_dzp/(1.0/curr_delNL0+growth_factor_zp);

                dadia_dzp *= (2.0/3.0)*T;

                // next heating due to the changing species
                dspec_dzp = - dxe_dzp * T / (1+x_e);

                // next, Compton heating
                //                dcomp_dzp = dT_comp(zp, T, x_e);
                dcomp_dzp = dcomp_dzp_prefactor*(x_e/(1.0+x_e+f_He))*( Trad_fast - T );

                // lastly, X-ray heating
                dxheat_dzp = dxheat_dt * dt_dzp * 2.0 / 3.0 / k_B / (1.0+x_e);

                //update quantities
                x_e += ( dxe_dzp ) * dzp; // remember dzp is negative
                if (x_e > 1) // can do this late in evolution if dzp is too large
                    x_e = 1 - FRACT_FLOAT_ERR;
                else if (x_e < 0)
                    x_e = 0;
                if (T < MAX_TK) {
                    T += ( dxheat_dzp + dcomp_dzp + dspec_dzp + dadia_dzp ) * dzp;
                }

                if (T<=0){ // spurious bahaviour of the trapazoidalintegrator. generally overcooling in underdensities
                    T = T_cmb*(1+zp);
                }

                x_e_box[box_ct] = x_e;
                Tk_box[box_ct] = T;

                if (COMPUTE_Ts){
                    J_alpha_tot = ( dxlya_dt + dstarlya_dt ); //not really d/dz, but the lya flux

                    // Note: to make the code run faster, the get_Ts function call to evaluate the spin temperature was replaced with the code below.
                    // Algorithm is the same, but written to be more computationally efficient
                    T_inv = pow(T,-1.);
                    T_inv_sq = pow(T,-2.);

                    xc_fast = (1.0+curr_delNL0*growth_factor_zp)*xc_inverse*( (1.0-x_e)*No*kappa_10_float(T,0) + x_e*N_b0*kappa_10_elec_float(T,0) + x_e*No*kappa_10_pH_float(T,0) );
                    xi_power = TS_prefactor * pow((1.0+curr_delNL0*growth_factor_zp)*(1.0-x_e)*T_inv_sq, 1.0/3.0);
                    xa_tilde_fast_arg = xa_tilde_prefactor*J_alpha_tot*pow( 1.0 + 2.98394*xi_power + 1.53583*pow(xi_power,2.) + 3.85289*pow(xi_power,3.), -1. );

                    if (J_alpha_tot > 1.0e-20) { // Must use WF effect
                        TS_fast = Trad_fast;
                        TSold_fast = 0.0;
                        while (fabs(TS_fast-TSold_fast)/TS_fast > 1.0e-3) {

                            TSold_fast = TS_fast;

                            xa_tilde_fast = ( 1.0 - 0.0631789*T_inv + 0.115995*T_inv_sq - 0.401403*T_inv*pow(TS_fast,-1.) + 0.336463*T_inv_sq*pow(TS_fast,-1.) )*xa_tilde_fast_arg;

                            TS_fast = (1.0+xa_tilde_fast+xc_fast)*pow(Trad_fast_inv+xa_tilde_fast*( T_inv + 0.405535*T_inv*pow(TS_fast,-1.) - 0.405535*T_inv_sq ) + xc_fast*T_inv,-1.);
                        }
                    } else { // Collisions only
                        TS_fast = (1.0 + xc_fast)/(Trad_fast_inv + xc_fast*T_inv);
                        xa_tilde_fast = 0.0;
                    }
                    if(TS_fast < 0.) {
                        // It can very rarely result in a negative spin temperature. If negative, it is a very small number. Take the absolute value, the optical depth can deal with very large numbers, so ok to be small
                        //printf("Tiny Ts=%.1e \n",TS_fast);
                        TS_fast = fabs(TS_fast);
                    }

                    Ts[box_ct] = TS_fast;

                    if(OUTPUT_AVE || STORE_DATA) {
                        J_alpha_ave += J_alpha_tot;
                        xalpha_ave += xa_tilde_fast;
                        Xheat_ave += ( dxheat_dzp );
                        Xion_ave += ( dt_dzp*dxion_source_dt );

                        Ts_ave += TS_fast;
                        Tk_ave += T;
                        Tk_inv_ave += 1.0/T;
                        Tk_inv_sq_ave += 1.0/TS_fast;

                        
                    }
                }
                x_e_ave += x_e;
            }
            // For this redshift snapshot, we now determine the ionisation field and subsequently the 21cm brightness temperature map (also the 21cm PS)
            // Note the relatively small tolerance for zp and the input redshift. The user needs to be careful to provide the correct redshifts for evaluating this to high precision.
            // If the light-cone option is set, this criterion should automatically be met
            for(i_z=0;i_z<N_USER_REDSHIFT;i_z++) {
                if(fabs(redshifts[i_z] - zp)<0.001) {

                    memcpy(Ts_z,Ts,sizeof(float)*HII_TOT_NUM_PIXELS);
                    memcpy(x_e_z,x_e_box,sizeof(float)*HII_TOT_NUM_PIXELS);

                    if(i_z==0) {
                        // If in here, it doesn't matter what PREV_REDSHIFT is set to
                        // as the recombinations will not be calculated
                        ComputeIonisationBoxes(i_z,redshifts[i_z],redshifts[i_z]+0.2);
                    }
                    else {
                        ComputeIonisationBoxes(i_z,redshifts[i_z],redshifts[i_z-1]);
                    }
                }
            }

            /////////////////////////////  END LOOP ////////////////////////////////////////////

            // compute new average values
            x_e_ave /= (double)HII_TOT_NUM_PIXELS;


            Ts_ave /= (double)HII_TOT_NUM_PIXELS;
            Tk_ave /= (double)HII_TOT_NUM_PIXELS;
            Tk_inv_ave /= (double)HII_TOT_NUM_PIXELS;
            Tk_inv_sq_ave /= (double)HII_TOT_NUM_PIXELS;
            J_alpha_ave /= (double)HII_TOT_NUM_PIXELS;
            xalpha_ave /= (double)HII_TOT_NUM_PIXELS;
            Xheat_ave /= (double)HII_TOT_NUM_PIXELS;
            Xion_ave /= (double)HII_TOT_NUM_PIXELS;
            T21_ave /= (double)HII_TOT_NUM_PIXELS;
              if(OUTPUT_AVE) {
                printf("zp = %e T21_avg= %e [mK] Ts_ave = %e x_e_ave = %e Tk_ave = %e J_alpha_ave = %e xalpha_ave = %e Xheat_ave = %e Xion_ave = %e\n",
                        zp, T21_ave ,Ts_ave,x_e_ave,Tk_ave,J_alpha_ave,xalpha_ave,Xheat_ave,Xion_ave);
              }
            //JBM:also save Tspin and Tkin averages
            if(STORE_DATA) {
              aveTkin[counter] = Tk_ave;
              aveTspin_inv[counter] = Tk_inv_ave;
              aveTkin_inv_sq[counter] = Tk_inv_sq_ave;
              aveTspin[counter] = Ts_ave;
              aveXalpha[counter] = xalpha_ave;
            }

            prev_zp = zp;
            zp = ((1+prev_zp) / ZPRIME_STEP_FACTOR - 1);
            dzp = zp - prev_zp;

            counter += 1;
        } // end main integral loop over z'

        // This memory freeing is problematic so is moved until after the data is saved
        // fprintf(stderr, "Start freeing Ts memory\n");
        // destroy_21cmMC_Ts_arrays();
        // fprintf(stderr, "Freed Ts arrays, Starting destruct_heat()\n");
        // destruct_heat();
    }

    fprintf(stderr, "Start freeing interpolation table memory\n");
    for(i=0;i<Numzp_for_table;i++) {
        for(j=0;j<X_RAY_Tvir_POINTS;j++) {
            SAFE_FREE(Fcoll_R_Table[i][j]);
        }
        SAFE_FREE(Fcoll_R_Table[i]);
    }
    SAFE_FREE(Fcoll_R_Table);
}

void free_TS_memory() {
    fprintf(stderr, "Start freeing Ts memory\n");
    destroy_21cmMC_Ts_arrays();

    fprintf(stderr, "Freed Ts arrays, Starting destruct_heat()\n");
    destruct_heat();
}

void ComputeIonisationBoxes(int sample_index, float REDSHIFT_SAMPLE, float PREV_REDSHIFT) {

    /* This is an entire re-write of find_HII_bubbles.c from 21cmFAST. Refer back to that code if this become a little confusing, the computation and algorithm are the same.
     Additionally, the code here includes delta_T.c for calculating the 21cm PS, and also redshift_interpolate_boxes.c for calculating the lightcones. */

    char filename[500];
    FILE *F;
    fftwf_plan plan;

    // Other parameters used in the code
    int i,j,k,ii, x,y,z, N_min_cell, LAST_FILTER_STEP, short_completely_ionised,skip_deallocate,first_step_R;
    int n_x, n_y, n_z,counter,LOOP_INDEX;
    unsigned long long ct;

    float growth_factor,MFEEDBACK, R, pixel_mass, cell_length_factor, ave_N_min_cell, M_MIN, nf;
    float f_coll_crit, erfc_denom, erfc_denom_cell, res_xH, Splined_Fcoll, sqrtarg, xHI_from_xrays, curr_dens, stored_R, massofscaleR;

    double global_xH, global_step_xH, ST_over_PS, mean_f_coll_st, f_coll, f_coll_temp, f_coll_from_table, f_coll_from_table_1, f_coll_from_table_2;

    double t_ast, dfcolldt, Gamma_R_prefactor, rec, dNrec;
    float growth_factor_dz, fabs_dtdz, ZSTEP, Gamma_R, z_eff;
    const float dz = 0.01;

    const gsl_rng_type * T;
    gsl_rng * r;

    skip_deallocate = 0;

    // Choice of DIM is arbitrary, just needs to be a value larger than HII_DIM. DIM should be sufficient as it shouldn't exceeded DIM (provided DIM > HII_DIM by a factor of at least ~3)
    int *LOS_index = calloc(DIM,sizeof(int));
    int *slice_index = calloc(DIM,sizeof(int));

    int total_in_z = 0;

    float d1_low, d1_high, d2_low, d2_high, gradient_component, min_gradient_component, subcell_width, x_val1, x_val2, subcell_displacement;
    float RSD_pos_new, RSD_pos_new_boundary_low,RSD_pos_new_boundary_high, fraction_within, fraction_outside, cell_distance;

    int min_slice_index,slice_index_reducedLC;

    min_slice_index = HII_DIM + 1;
    slice_index_reducedLC = 0;

    // For recombinations
    if(INHOMO_RECO) {
        ZSTEP = PREV_REDSHIFT - REDSHIFT_SAMPLE;

        for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++) {
            Gamma12[ct] = 0.0;
        }
    }
    else {
        ZSTEP = 0.2;
    }
    fabs_dtdz = fabs(dtdz(REDSHIFT_SAMPLE));
    t_ast = t_STAR * t_hubble(REDSHIFT_SAMPLE);
    growth_factor_dz = dicke(REDSHIFT_SAMPLE-dz);

    // if USE_FCOLL_IONISATION_TABLE == 1, we are only calculating the ionisation fraction within a smaller volume (however much is required to linearly interpolate the fields for the light-cone).
    // To know which slices to keep and discard, need to store the relevant indices etc. for the given redshift. The below code does this.
    if(USE_FCOLL_IONISATION_TABLE) {

        total_in_z = 0;

        // LOS_direction is used in two separate locations (iterated), so we store the state for each individucal usage
        LOS_direction = Stored_LOS_direction_state_2;

        // First store the requisite slices for this redshift snapshot.
        if(start_index_LC[N_USER_REDSHIFT-sample_index-1] > end_index_LC[N_USER_REDSHIFT-sample_index-1]) {
            for(ii=0;ii<end_index_LC[N_USER_REDSHIFT-sample_index-1];ii++) {
                LOS_index[total_in_z] = LOS_direction;
                slice_index[total_in_z] = ii;
                total_in_z += 1;
            }
            for(ii=start_index_LC[N_USER_REDSHIFT-sample_index-1];ii<HII_DIM;ii++) {
                LOS_index[total_in_z] = LOS_direction;
                slice_index[total_in_z] = ii;
                total_in_z += 1;
            }
        }
        else {
            for(ii=start_index_LC[N_USER_REDSHIFT-sample_index-1];ii<end_index_LC[N_USER_REDSHIFT-sample_index-1];ii++) {
                LOS_index[total_in_z] = LOS_direction;
                slice_index[total_in_z] = ii;
                total_in_z += 1;
            }
        }

        if(sample_index<(N_USER_REDSHIFT-1)) {

            // Now store the requisite slices for the next redshift (lower redshift) snapshot
            if(start_index_LC[N_USER_REDSHIFT-(sample_index+1)-1] > end_index_LC[N_USER_REDSHIFT-(sample_index+1)-1]) {
                for(ii=0;ii<end_index_LC[N_USER_REDSHIFT-(sample_index+1)-1];ii++) {
                    LOS_index[total_in_z] = LOS_direction;
                    slice_index[total_in_z] = ii;
                    total_in_z += 1;
                }

                for(ii=start_index_LC[N_USER_REDSHIFT-(sample_index+1)-1];ii<HII_DIM;ii++) {
                    LOS_index[total_in_z] = LOS_direction;
                    slice_index[total_in_z] = ii;
                    total_in_z += 1;
                }
            }
            else {
                for(ii=start_index_LC[N_USER_REDSHIFT-(sample_index+1)-1];ii<end_index_LC[N_USER_REDSHIFT-(sample_index+1)-1];ii++) {
                    LOS_index[total_in_z] = LOS_direction;
                    slice_index[total_in_z] = ii;
                    total_in_z += 1;
                }
            }
        }
        if(total_in_z > HII_DIM) {
            // The entire box is used (i.e. we haven't reduced the computation), since the number of required slices between two subsequent redshifts exceeds HII_DIM (i.e. boxes too small).
            total_in_z = HII_DIM;
            for(ii=0;ii<HII_DIM;ii++) {
                slice_index[ii] = ii;
                LOS_index[ii] = Default_LOS_direction;
            }
        }

        if((N_USER_REDSHIFT-(sample_index)-1)>0) {
            min_slice_index = start_index_LC[N_USER_REDSHIFT-(sample_index+1)-1];
        }
        else {
            min_slice_index = start_index_LC[N_USER_REDSHIFT-sample_index-1];
        }

        if(SUBCELL_RSD) {
            // Add in the padding to the calculation of the ionisation field to account for the fact that cells can enter/exit the reduced LC box.
            min_slice_index = min_slice_index - LC_BOX_PADDING;
            if(min_slice_index < 0) {
                min_slice_index = min_slice_index + HII_DIM;
            }
            total_in_z += 2*LC_BOX_PADDING;

            if(total_in_z > HII_DIM) {
                // if after adding in the padding, the total is larger than the box length, then only need to sample at the box length (HII_DIM) and not greater
                total_in_z = HII_DIM;
            }
        }
    }

    /////////////////////////////////   BEGIN INITIALIZATION   //////////////////////////////////

    // perform a very rudimentary check to see if we are underresolved and not using the linear approx
    if ((BOX_LEN > DIM) && !EVOLVE_DENSITY_LINEARLY){
    printf("perturb_field.c: WARNING: Resolution is likely too low for accurate evolved density fields\n It Is recommended that you either increase the resolution (DIM/Box_LEN) or set the EVOLVE_DENSITY_LINEARLY flag to 1\n");
    }

    // initialize power spectrum
    growth_factor = dicke(REDSHIFT_SAMPLE);

    init_21cmMC_HII_arrays();
    if(GenerateNewICs) {

        // Calculate the density field for this redshift if the initial conditions/cosmology are changing
        ComputePerturbField(REDSHIFT_SAMPLE);

        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    *((float *)deltax_unfiltered + HII_R_FFT_INDEX(i,j,k)) = LOWRES_density_REDSHIFT[HII_R_INDEX(i,j,k)];
                }
            }
        }

    }
    else {
        // Read the desnity field of this redshift from file
        sprintf(filename, "%s/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc",BOXES_INPUT_FOLDER, REDSHIFT_SAMPLE, HII_DIM, BOX_LEN);
        F = fopen(filename, "rb");
        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    if (fread((float *)deltax_unfiltered + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
                        printf("Read error occured while reading deltax box.\n");
                    }
                }
            }
        }
        fclose(F);
    }

    // keep the unfiltered density field in an array, to save it for later
    memcpy(deltax_unfiltered_original, deltax_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

    i=0;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    pixel_mass = RtoM(L_FACTOR*BOX_LEN/(float)HII_DIM);
//    f_coll_crit = 1/HII_EFF_FACTOR;
    cell_length_factor = L_FACTOR;

    //set the minimum source mass
    if (ION_Tvir_MIN > 0){ // use the virial temperature for Mmin
        if (ION_Tvir_MIN < 9.99999e3) // neutral IGM
            M_MIN = TtoM(REDSHIFT_SAMPLE, ION_Tvir_MIN, 1.22);
        else // ionized IGM
            M_MIN = TtoM(REDSHIFT_SAMPLE, ION_Tvir_MIN, 0.6);
    }
    else if (ION_Tvir_MIN < 0){ // use the mass
        M_MIN = ION_M_MIN;
    }
    // check for WDM

    if (P_CUTOFF && ( M_MIN < M_J_WDM())){
        printf( "The default Jeans mass of %e Msun is smaller than the scale supressed by the effective pressure of WDM.\n", M_MIN);
        M_MIN = M_J_WDM();
        printf( "Setting a new effective Jeans mass from WDM pressure supression of %e Msun\n", M_MIN);
    }

    for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
        xH[ct] = 1.;
    }

    // lets check if we are going to bother with computing the inhmogeneous field at all...

    global_xH = 0.0;

    MFEEDBACK = M_MIN;

    if(EFF_FACTOR_PL_INDEX != 0.) {
      mean_f_coll_st = FgtrM_st_PL(REDSHIFT_SAMPLE,M_MIN,MFEEDBACK,EFF_FACTOR_PL_INDEX);

    }
    else {

        mean_f_coll_st = FgtrM_st(REDSHIFT_SAMPLE, M_MIN);
    }
    if (mean_f_coll_st/(1./HII_EFF_FACTOR) < HII_ROUND_ERR){ // way too small to ionize anything...
//        printf( "The ST mean collapse fraction is %e, which is much smaller than the effective critical collapse fraction of %e\n I will just declare everything to be neutral\n", mean_f_coll_st, f_coll_crit);

        // find the neutral fraction
        if(USE_TS_FLUCT) {
            for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
                xH[ct] = 1-x_e_z[ct]; // convert from x_e to xH
                global_xH += xH[ct];
            }
            global_xH /= (double)HII_TOT_NUM_PIXELS;
        }
        else {
            init_heat();
            global_xH = 1. - xion_RECFAST(REDSHIFT_SAMPLE, 0);
            destruct_heat();
            for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
                xH[ct] = global_xH;
            }
        }
        skip_deallocate = 1;

    }
    else {


        // Take the ionisation fraction from the X-ray ionisations from Ts.c (only if the calculate spin temperature flag is set)
        if(USE_TS_FLUCT) {
            for (i=0; i<HII_DIM; i++){
                for (j=0; j<HII_DIM; j++){
                    for (k=0; k<HII_DIM; k++){
                        *((float *)xe_unfiltered + HII_R_FFT_INDEX(i,j,k)) = x_e_z[HII_R_INDEX(i,j,k)];
                    }
                }
            }
        }

        if(USE_TS_FLUCT) {
            plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)xe_unfiltered, (fftwf_complex *)xe_unfiltered, FFTW_ESTIMATE);
            fftwf_execute(plan);
        }

        if (INHOMO_RECO){
            plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)N_rec_unfiltered, (fftwf_complex *)N_rec_unfiltered, FFTW_ESTIMATE);
            fftwf_execute(plan);
        }

        plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deltax_unfiltered, (fftwf_complex *)deltax_unfiltered, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
        fftwf_cleanup();
        // remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from
        //  real space to k-space
        // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below

        for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
            deltax_unfiltered[ct] /= (HII_TOT_NUM_PIXELS+0.0);
        }

        if(USE_TS_FLUCT) {
            for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
                xe_unfiltered[ct] /= (double)HII_TOT_NUM_PIXELS;
            }
        }

        if (INHOMO_RECO){
            for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
                N_rec_unfiltered[ct] /= (double)HII_TOT_NUM_PIXELS;
            }
        }

        printf("boxes FFT'd!\n ");

        /*************************************************************************************/
        /***************** LOOP THROUGH THE FILTER RADII (in Mpc)  ***************************/
        /*************************************************************************************/
        // set the max radius we will use, making sure we are always sampling the same values of radius
        // (this avoids aliasing differences w redshift)

        int determine_R_intermediate;

        determine_R_intermediate = 0;

        short_completely_ionised = 0;
        // loop through the filter radii (in Mpc)
        erfc_denom_cell=1; //dummy value

        R=fmax(R_BUBBLE_MIN, (cell_length_factor*BOX_LEN/(float)HII_DIM));

        int N_RSTEPS, counter_R;

        counter = 0;
        while ((R - fmin(R_BUBBLE_MAX, L_FACTOR*BOX_LEN)) <= FRACT_FLOAT_ERR ) {
            R*= DELTA_R_HII_FACTOR;
            determine_R_intermediate += 1;
            if(R >= fmin(R_BUBBLE_MAX, L_FACTOR*BOX_LEN)) {
                stored_R = R/DELTA_R_HII_FACTOR;
            }
            counter += 1;
        }

        N_RSTEPS = counter;

        counter_R = N_RSTEPS;

        determine_R_intermediate = determine_R_intermediate - 2;

        R=fmin(R_BUBBLE_MAX, L_FACTOR*BOX_LEN);

        LAST_FILTER_STEP = 0;

        initialiseSplinedSigmaM(M_MIN,1e16);

        first_step_R = 1;

        counter = 0;

        while (!LAST_FILTER_STEP && (M_MIN < RtoM(R)) ){

            // Check if we are the last filter step
            if ( ((R/DELTA_R_HII_FACTOR - cell_length_factor*BOX_LEN/(float)HII_DIM) <= FRACT_FLOAT_ERR) || ((R/DELTA_R_HII_FACTOR - R_BUBBLE_MIN) <= FRACT_FLOAT_ERR) ) {
                LAST_FILTER_STEP = 1;
                R = fmax(cell_length_factor*BOX_LEN/(double)HII_DIM, R_BUBBLE_MIN);
            }

            // Copy all relevant quantities from memory into new arrays to be smoothed and FFT'd.
            if(USE_TS_FLUCT) {
                memcpy(xe_filtered, xe_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
            }
            if (INHOMO_RECO){
                memcpy(N_rec_filtered, N_rec_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
            }
            memcpy(deltax_filtered, deltax_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

            if (!LAST_FILTER_STEP || (R > cell_length_factor*BOX_LEN/(double)HII_DIM) ){
                if(USE_TS_FLUCT) {
                    HII_filter(xe_filtered, HII_FILTER, R);
                }
                if (INHOMO_RECO){
                    HII_filter(N_rec_filtered, HII_FILTER, R);
                }
                HII_filter(deltax_filtered, HII_FILTER, R);
            }

            // Perform FFTs
            plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)deltax_filtered, (float *)deltax_filtered, FFTW_ESTIMATE);
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);
            fftwf_cleanup();

            if (USE_TS_FLUCT) {
                plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)xe_filtered, (float *)xe_filtered, FFTW_ESTIMATE);
                fftwf_execute(plan);
            }

            if (INHOMO_RECO){
                plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)N_rec_filtered, (float *)N_rec_filtered, FFTW_ESTIMATE);
                fftwf_execute(plan);
            }

            // printf("More FFTs! \n\n");


            // Check if this is the last filtering scale.  If so, we don't need deltax_unfiltered anymore.
            // We will re-read it to get the real-space field, which we will use to set the residual neutral fraction
            ST_over_PS = 0;
            f_coll = 0;
            massofscaleR = RtoM(R);

            erfc_denom = 2.*(pow(sigma_z0(M_MIN), 2) - pow(sigma_z0(massofscaleR), 2) );
            if (erfc_denom < 0) { // our filtering scale has become too small
                break;
            }
            erfc_denom = sqrt(erfc_denom);
            erfc_denom = 1./( growth_factor * erfc_denom );


            if(EFF_FACTOR_PL_INDEX!=0.) {
                initialiseGL_Fcoll(NGLlow,NGLhigh,M_MIN,massofscaleR);
                initialiseFcoll_spline(REDSHIFT_SAMPLE,M_MIN,massofscaleR,massofscaleR,MFEEDBACK,EFF_FACTOR_PL_INDEX);
            }

            if(!USE_FCOLL_IONISATION_TABLE) {

                // Determine the global averaged f_coll for the overall normalisation

                // renormalize the collapse fraction so that the mean matches ST,
                // since we are using the evolved (non-linear) density field
                for (x=0; x<HII_DIM; x++){
                    for (y=0; y<HII_DIM; y++){
                        for (z=0; z<HII_DIM; z++){

                            // delta cannot be less than -1
                            *((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z)) = FMAX(*((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z)) , -1.+FRACT_FLOAT_ERR);

                            // <N_rec> cannot be less than zero
                            if (INHOMO_RECO){
                                *((float *)N_rec_filtered + HII_R_FFT_INDEX(x,y,z)) = FMAX(*((float *)N_rec_filtered + HII_R_FFT_INDEX(x,y,z)) , 0.0);
                            }

                            // x_e has to be between zero and unity
                            if (USE_TS_IN_21CM){
                                *((float *)xe_filtered + HII_R_FFT_INDEX(x,y,z)) = FMAX(*((float *)xe_filtered + HII_R_FFT_INDEX(x,y,z)) , 0.);
                                *((float *)xe_filtered + HII_R_FFT_INDEX(x,y,z)) = FMIN(*((float *)xe_filtered + HII_R_FFT_INDEX(x,y,z)) , 0.999);
                            }

                            curr_dens = *((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z));

                            if(EFF_FACTOR_PL_INDEX!=0.) {
                                // Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
                                // However, such densities should always be collapsed, so just set f_coll to unity. Additionally, the fraction of points in this regime relative
                                // to the entire simulation volume is extremely small.
                                if(curr_dens <= 0.99*Deltac) {
                                    FcollSpline(curr_dens,&(Splined_Fcoll));
                                }
                                else { // the entrire cell belongs to a collpased halo...  this is rare...
                                    Splined_Fcoll =  1.0;
                                }
                            }
                            else {

                                erfc_arg_val = (Deltac - curr_dens)*erfc_denom;
                                if( erfc_arg_val < erfc_arg_min || erfc_arg_val > erfc_arg_max ) {
                                    Splined_Fcoll = splined_erfc(erfc_arg_val);
                                }
                                else {
                                    erfc_arg_val_index = (int)floor(( erfc_arg_val - erfc_arg_min )*InvArgBinWidth);
                                    Splined_Fcoll = ERFC_VALS[erfc_arg_val_index] + (erfc_arg_val - (erfc_arg_min + ArgBinWidth*(double)erfc_arg_val_index))*ERFC_VALS_DIFF[erfc_arg_val_index]*InvArgBinWidth;
                                }
                            }

                            // save the value of the collasped fraction into the Fcoll array
                            Fcoll[HII_R_INDEX(x,y,z)] = Splined_Fcoll;
                            f_coll += Splined_Fcoll;
                        }
                    }
                } //  end loop through Fcoll box

                f_coll /= (double) HII_TOT_NUM_PIXELS;
            }
            else {

                // Evaluate the interpolation table of the global average f_coll

                R_MFP_INT_1 = (int)floor((R - R_MFP_MIN)/R_MFP_BINWIDTH);
                R_MFP_INT_2 = R_MFP_INT_1 + 1;

                R_MFP_VAL_1 = R_MFP_MIN + (R_MFP_UB - R_MFP_MIN)*(float)R_MFP_INT_1/((float)R_MFP_STEPS - 1.);
                R_MFP_VAL_2 = R_MFP_MIN + (R_MFP_UB - R_MFP_MIN)*(float)R_MFP_INT_2/((float)R_MFP_STEPS - 1.);

                if (LAST_FILTER_STEP){
                    f_coll = Ionisation_fcoll_table_final[TVIR_INT_1 + TVIR_STEPS*sample_index] + ( log10(ION_Tvir_MIN) - TVIR_VAL_1 )*( Ionisation_fcoll_table_final[TVIR_INT_2 + TVIR_STEPS*sample_index] - Ionisation_fcoll_table_final[TVIR_INT_1 + TVIR_STEPS*sample_index] )/( TVIR_VAL_2 - TVIR_VAL_1 );
                }
                else {

                    f_coll_from_table_1 = ( R_MFP_VAL_2 - R )*Ionisation_fcoll_table[TVIR_INT_1 + TVIR_STEPS*( R_MFP_INT_1 + R_MFP_STEPS*sample_index )];
                    f_coll_from_table_1 += ( R - R_MFP_VAL_1 )*Ionisation_fcoll_table[TVIR_INT_1 + TVIR_STEPS*( R_MFP_INT_2 + R_MFP_STEPS*sample_index )];
                    f_coll_from_table_1 /= ( R_MFP_VAL_2 - R_MFP_VAL_1 );

                    f_coll_from_table_2 = ( R_MFP_VAL_2 - R )*Ionisation_fcoll_table[TVIR_INT_2 + TVIR_STEPS*( R_MFP_INT_1 + R_MFP_STEPS*sample_index )];
                    f_coll_from_table_2 += ( R - R_MFP_VAL_1 )*Ionisation_fcoll_table[TVIR_INT_2 + TVIR_STEPS*( R_MFP_INT_2 + R_MFP_STEPS*sample_index )];
                    f_coll_from_table_2 /= ( R_MFP_VAL_2 - R_MFP_VAL_1 );

                    f_coll = ( ( TVIR_VAL_2 - log10(ION_Tvir_MIN) )*f_coll_from_table_1 + ( log10(ION_Tvir_MIN) - TVIR_VAL_1 )*f_coll_from_table_2 )/( TVIR_VAL_2 - TVIR_VAL_1 );
                }
            }


            ST_over_PS = mean_f_coll_st/f_coll;

            //////////////////////////////  MAIN LOOP THROUGH THE BOX ///////////////////////////////////
            // now lets scroll through the filtered box

            rec = 0.;

            xHI_from_xrays = 1;
            Gamma_R_prefactor = pow(1+REDSHIFT_SAMPLE, 2) * (R*CMperMPC) * SIGMA_HI * ALPHA_UVB / (ALPHA_UVB+2.75) * N_b0 * HII_EFF_FACTOR / 1.0e-12;

            if(!USE_FCOLL_IONISATION_TABLE) {
                LOOP_INDEX = HII_DIM;
            }
            else {
                LOOP_INDEX = total_in_z;
            }

            for (x=0; x<HII_DIM; x++){
                for (y=0; y<HII_DIM; y++){
                    for (z=0; z<LOOP_INDEX; z++){

                        if(USE_FCOLL_IONISATION_TABLE) {
                            if((min_slice_index + z) >= HII_DIM) {
                                slice_index_reducedLC = (min_slice_index + z) - HII_DIM;
                            }
                            else {
                                slice_index_reducedLC = (min_slice_index + z);
                            }

                            // delta cannot be less than -1
                            *((float *)deltax_filtered + coeval_box_pos_FFT(Default_LOS_direction,x,y,slice_index_reducedLC)) = FMAX(*((float *)deltax_filtered + coeval_box_pos_FFT(Default_LOS_direction,x,y,slice_index_reducedLC)) , -1.+FRACT_FLOAT_ERR);

                            // <N_rec> cannot be less than zero
                            if (INHOMO_RECO){
                                *((float *)N_rec_filtered + coeval_box_pos_FFT(Default_LOS_direction,x,y,slice_index_reducedLC)) = FMAX(*((float *)N_rec_filtered + coeval_box_pos_FFT(Default_LOS_direction,x,y,slice_index_reducedLC)) , 0.0);
                            }

                            // x_e has to be between zero and unity
                            if (USE_TS_IN_21CM){
                                *((float *)xe_filtered + coeval_box_pos_FFT(Default_LOS_direction,x,y,slice_index_reducedLC)) = FMAX(*((float *)xe_filtered + coeval_box_pos_FFT(Default_LOS_direction,x,y,slice_index_reducedLC)) , 0.);
                                *((float *)xe_filtered + coeval_box_pos_FFT(Default_LOS_direction,x,y,slice_index_reducedLC)) = FMIN(*((float *)xe_filtered + coeval_box_pos_FFT(Default_LOS_direction,x,y,slice_index_reducedLC)) , 0.999);
                            }
                        }
                        else {
                            slice_index_reducedLC = z;
                        }

                        curr_dens = *((float *)deltax_filtered + coeval_box_pos_FFT(Default_LOS_direction,x,y,slice_index_reducedLC));

                        if(USE_FCOLL_IONISATION_TABLE) {

                            if(EFF_FACTOR_PL_INDEX!=0.) {
                                if(curr_dens <= 0.99*Deltac) {
                                    // This is here as the interpolation tables have some issues very close
                                    // to Deltac. So lets just assume these voxels collapse anyway.
                                    FcollSpline(curr_dens,&(Splined_Fcoll));
                                    // Using the fcoll ionisation table, fcoll is not computed in each cell.
                                    // For this option, fcoll must be calculated.
                                }
                                else {
                                    Splined_Fcoll = 1.;
                                }
                            }
                            // check for aliasing which can occur for small R and small cell sizes,
                            // since we are using the analytic form of the window function for speed and simplicity
                            else {

                                erfc_arg_val = (Deltac - curr_dens)*erfc_denom;
                                if( erfc_arg_val < erfc_arg_min || erfc_arg_val > erfc_arg_max ) {
                                    Splined_Fcoll = splined_erfc(erfc_arg_val);
                                }
                                else {
                                    erfc_arg_val_index = (int)floor(( erfc_arg_val - erfc_arg_min )*InvArgBinWidth);
                                    Splined_Fcoll = ERFC_VALS[erfc_arg_val_index] + (erfc_arg_val - (erfc_arg_min + ArgBinWidth*(double)erfc_arg_val_index))*ERFC_VALS_DIFF[erfc_arg_val_index]*InvArgBinWidth;
                                }
                            }
                        }
                        else {

                            Splined_Fcoll = Fcoll[coeval_box_pos(Default_LOS_direction,x,y,slice_index_reducedLC)];
                        }

                        f_coll = ST_over_PS * Splined_Fcoll;

                        if (INHOMO_RECO){
                            dfcolldt = f_coll / t_ast;
                            Gamma_R = Gamma_R_prefactor * dfcolldt;
                            rec = (*((float *)N_rec_filtered + coeval_box_pos_FFT(Default_LOS_direction,x,y,slice_index_reducedLC))); // number of recombinations per mean baryon
                            rec /= (1. + curr_dens); // number of recombinations per baryon inside <R>
                        }

                        // adjust the denominator of the collapse fraction for the residual electron fraction in the neutral medium
                        if (USE_TS_IN_21CM){
                            xHI_from_xrays = (1. - *((float *)xe_filtered + coeval_box_pos_FFT(Default_LOS_direction,x,y,slice_index_reducedLC)));
                        }

                        // check if fully ionized!
                        if ( (f_coll > (xHI_from_xrays/HII_EFF_FACTOR)*(1.0+rec)) ){ //IONIZED!!

                            // if this is the first crossing of the ionization barrier for this cell (largest R), record the gamma
                            // this assumes photon-starved growth of HII regions...  breaks down post EoR
                            if (INHOMO_RECO && (xH[coeval_box_pos(Default_LOS_direction,x,y,slice_index_reducedLC)] > FRACT_FLOAT_ERR) ){
                                Gamma12[coeval_box_pos(Default_LOS_direction,x,y,slice_index_reducedLC)] = Gamma_R;
                            }

                            // keep track of the first time this cell is ionized (earliest time)
                            if (INHOMO_RECO && (z_re[coeval_box_pos(Default_LOS_direction,x,y,slice_index_reducedLC)] < 0)){
                                z_re[coeval_box_pos(Default_LOS_direction,x,y,slice_index_reducedLC)] = REDSHIFT_SAMPLE;
                            }

                            // FLAG CELL(S) AS IONIZED
                            if (FIND_BUBBLE_ALGORITHM == 2) // center method
                                xH[coeval_box_pos(Default_LOS_direction,x,y,slice_index_reducedLC)] = 0;
                            else if (FIND_BUBBLE_ALGORITHM == 1) // sphere method
                                update_in_sphere(xH, HII_DIM, R/BOX_LEN, x/(HII_DIM+0.0), y/(HII_DIM+0.0), slice_index_reducedLC/(HII_DIM+0.0));
                            else{
                                printf( "Incorrect choice of find bubble algorithm: %i\nAborting...", FIND_BUBBLE_ALGORITHM);
                                xH[coeval_box_pos(Default_LOS_direction,x,y,slice_index_reducedLC)] = 0;
                            }
                        } // end ionized
                        // If not fully ionized, then assign partial ionizations
                        else if (LAST_FILTER_STEP && (xH[coeval_box_pos(Default_LOS_direction,x,y,slice_index_reducedLC)] > TINY)){

                            if (f_coll>1) f_coll=1;

                            ave_N_min_cell = f_coll * pixel_mass*(1. + curr_dens) / M_MIN; // ave # of M_MIN halos in cell

                            if (ave_N_min_cell < N_POISSON){
                                N_min_cell = (int) gsl_ran_poisson(r, ave_N_min_cell);
                                f_coll = N_min_cell * M_MIN / (pixel_mass*(1. + curr_dens));
                            }

                            if (f_coll>1) f_coll=1;
                            res_xH = xHI_from_xrays - f_coll * HII_EFF_FACTOR;
                            // and make sure fraction doesn't blow up for underdense pixels
                            if (res_xH < 0)
                                res_xH = 0;
                            else if (res_xH > 1)
                                res_xH = 1;

                            xH[coeval_box_pos(Default_LOS_direction,x,y,slice_index_reducedLC)] = res_xH;
                        } // end partial ionizations at last filtering step
                    } // k
                } // j
            } // i

            if(!USE_FCOLL_IONISATION_TABLE) {

                global_step_xH = 0;
                for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
                    global_step_xH += xH[ct];
                }
                global_step_xH /= (float)HII_TOT_NUM_PIXELS;

                if(global_step_xH==0.0) {
                    short_completely_ionised = 1;
                    break;
                }
            }

            if(first_step_R) {
                R = stored_R;
                first_step_R = 0;
            }
            else {
                R /= DELTA_R_HII_FACTOR;
            }
            counter_R -= 1;

        }
        if(!USE_FCOLL_IONISATION_TABLE) {
            // find the neutral fraction
            global_xH = 0;

            for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
                global_xH += xH[ct];
            }
            global_xH /= (float)HII_TOT_NUM_PIXELS;
        }
        else {
            // Estimate the neutral fraction from the reduced box. Can be handy to have, but shouldn't be trusted for anything more as only a fraction of the co-eval box is being used
            global_xH = 0;
            for (x=0; x<HII_DIM; x++){
                for (y=0; y<HII_DIM; y++){
                    for (z=0; z<total_in_z; z++){

                        if((min_slice_index + z) >= HII_DIM) {
                            slice_index_reducedLC = (min_slice_index + z) - HII_DIM;
                        }
                        else {
                            slice_index_reducedLC = (min_slice_index + z);
                        }
                        global_xH += xH[coeval_box_pos(Default_LOS_direction,x,y,slice_index_reducedLC)];
                    }
                }
            }
            global_xH /= ((float)HII_DIM*(float)HII_DIM*(float)total_in_z);
        }

        // update the N_rec field
        if (INHOMO_RECO){

            //fft to get the real N_rec  and delta fields
            plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)N_rec_unfiltered, (float *)N_rec_unfiltered, FFTW_ESTIMATE);
            fftwf_execute(plan);
            plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)deltax_unfiltered, (float *)deltax_unfiltered, FFTW_ESTIMATE);

            fftwf_execute(plan);
            fftwf_destroy_plan(plan);

            fftwf_cleanup();

            for (x=0; x<HII_DIM; x++){
                for (y=0; y<HII_DIM; y++){
                    for (z=0; z<HII_DIM; z++){

                        curr_dens = 1.0 + (*((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z)));
                        z_eff = (1+REDSHIFT_SAMPLE) * pow(curr_dens, 1.0/3.0) - 1;
                        dNrec = splined_recombination_rate(z_eff, Gamma12[HII_R_INDEX(x,y,z)]) * fabs_dtdz * ZSTEP * (1 - xH[HII_R_INDEX(x,y,z)]);
                        *((float *)N_rec_unfiltered + HII_R_FFT_INDEX(x,y,z)) += dNrec;
                    }
                }
            }
        }
    }


    // deallocate
    gsl_rng_free (r);

    nf = global_xH;
    if(STORE_DATA) {
        aveNF[sample_index] = nf;
    }

    if(!USE_LIGHTCONE) {
        if(PRINT_FILES) {
            sprintf(filename, "%s/neutralFraction_%f_%f_%f.txt",OUTPUT_FOLDER,INDIVIDUAL_ID,INDIVIDUAL_ID_2,REDSHIFT_SAMPLE);
            F=fopen(filename, "wt");
            fprintf(F, "%lf\n",nf);
            fclose(F);
        }
    }
    ///////////////////////////////// End of perform 'find_HII_bubbles.c' ///////////////////////////////

    ///////////////////////////////////// Perform 'delta_T.c' ///////////////////////////////////////////

    float fz1, fz2, fz, z_slice;
    double dvdx, ave, max_v_deriv;
    float k_x, k_y, k_z, k_mag, const_factor, T_rad, pixel_Ts_factor, pixel_x_HI, pixel_deltax, H;

    ////////////////////////////////////  BEGIN INITIALIZATION //////////////////////////////////////////
    ave = 0;

    if(GenerateNewICs) {

        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    *((float *)v + HII_R_FFT_INDEX(i,j,k)) = LOWRES_velocity_REDSHIFT[HII_R_INDEX(i,j,k)];
                }
            }
        }
    }
    else {
        switch(VELOCITY_COMPONENT){
            case 1:  sprintf(filename, "%s/updated_vx_z%06.2f_%i_%.0fMpc", BOXES_INPUT_FOLDER, REDSHIFT_SAMPLE, HII_DIM, BOX_LEN);
                break;
            case 3:  sprintf(filename, "%s/updated_vz_z%06.2f_%i_%.0fMpc", BOXES_INPUT_FOLDER, REDSHIFT_SAMPLE, HII_DIM, BOX_LEN);
                break;
            default: sprintf(filename, "%s/updated_vy_z%06.2f_%i_%.0fMpc", BOXES_INPUT_FOLDER, REDSHIFT_SAMPLE, HII_DIM, BOX_LEN);
        }
        F=fopen(filename, "rb");
        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    if (fread((float *)v + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
                        printf("Read error occured while reading velocity box.\n");
                        fclose(F);
                    }
                }
            }
        }
        fclose(F);
    }

    T_rad = T_cmb*(1+REDSHIFT_SAMPLE);
    H = hubble(REDSHIFT_SAMPLE);
    const_factor = 27 * (OMb*hlittle*hlittle/0.023) *
    sqrt( (0.15/OMm/hlittle/hlittle) * (1+REDSHIFT_SAMPLE)/10.0 );

    memcpy(deltax, deltax_unfiltered_original, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

    ///////////////////////////////  END INITIALIZATION /////////////////////////////////////////////

    // ok, lets fill the delta_T box; which will be the same size as the bubble box

    for (i=0; i<HII_DIM; i++){
        for (j=0; j<HII_DIM; j++){
            for (k=0; k<HII_DIM; k++){

                pixel_deltax = deltax[HII_R_FFT_INDEX(i,j,k)];
                pixel_x_HI = xH[HII_R_INDEX(i,j,k)];

                delta_T[HII_R_INDEX(i,j,k)] = const_factor*pixel_x_HI*(1+pixel_deltax);

                if (USE_TS_FLUCT) {

                    if(SUBCELL_RSD) {
                        // Converting the prefactors into the optical depth, tau. Factor of 1000 is the conversion of spin temperature from K to mK
                        delta_T[HII_R_INDEX(i,j,k)] *= (1. + REDSHIFT_SAMPLE)/(1000.*Ts_z[HII_R_INDEX(i,j,k)]);
                    }
                    else {
                        pixel_Ts_factor = (1 - T_rad / Ts_z[HII_R_INDEX(i,j,k)]);
                        delta_T[HII_R_INDEX(i,j,k)] *= pixel_Ts_factor;
                    }
                }
                ave += delta_T[HII_R_INDEX(i,j,k)];
            }
        }
    }
    ave /= (float)HII_TOT_NUM_PIXELS;

    x_val1 = 0.;
    x_val2 = 1.;

    subcell_width = (BOX_LEN/(float)HII_DIM)/(float)N_RSD_STEPS;

    float max_cell_distance;

    max_cell_distance = 0.;

    // now write out the delta_T box
    if (T_USE_VELOCITIES){
        ave = 0.;

        memcpy(vel_gradient, v, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

        plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)vel_gradient, (fftwf_complex *)vel_gradient, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        for (n_x=0; n_x<HII_DIM; n_x++){
            if (n_x>HII_MIDDLE)
                k_x =(n_x-HII_DIM) * DELTA_K;  // wrap around for FFT convention
            else
                k_x = n_x * DELTA_K;

            for (n_y=0; n_y<HII_DIM; n_y++){
                if (n_y>HII_MIDDLE)
                    k_y =(n_y-HII_DIM) * DELTA_K;
                else
                    k_y = n_y * DELTA_K;

                for (n_z=0; n_z<=HII_MIDDLE; n_z++){
                    k_z = n_z * DELTA_K;

                    // take partial deriavative along the line of sight
                    switch(VELOCITY_COMPONENT){
                        case 1:
                            *((fftwf_complex *) vel_gradient + HII_C_INDEX(n_x,n_y,n_z)) *= k_x*I/(float)HII_TOT_NUM_PIXELS;
                            break;
                        case 3:
                            *((fftwf_complex *) vel_gradient + HII_C_INDEX(n_x,n_y,n_z)) *= k_z*I/(float)HII_TOT_NUM_PIXELS;
                            break;
                        default:
                            *((fftwf_complex *) vel_gradient + HII_C_INDEX(n_x,n_y,n_z)) *= k_y*I/(float)HII_TOT_NUM_PIXELS;
                    }
                }
            }
        }

        plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)vel_gradient, (float *)vel_gradient, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        if(SUBCELL_RSD) {

            // now add the velocity correction to the delta_T maps
            min_gradient_component = 1.0;

            for (i=0; i<HII_DIM; i++){
                for (j=0; j<HII_DIM; j++){

                    if(!USE_FCOLL_IONISATION_TABLE) {
                        LOOP_INDEX = HII_DIM;
                    }
                    else {
                        LOOP_INDEX = total_in_z;
                    }

                    for (k=0; k<LOOP_INDEX; k++){

                        if(!USE_FCOLL_IONISATION_TABLE) {
                            slice_index_reducedLC = k;
                        }
                        else {
                            // We are using the reduced box, with padding, so need to only fill the arrays in the necessary cells
                            if((min_slice_index + k) >= HII_DIM) {
                                slice_index_reducedLC = (min_slice_index + k) - HII_DIM;
                            }
                            else {
                                slice_index_reducedLC = (min_slice_index + k);
                            }
                        }

                        gradient_component = fabs(vel_gradient[HII_R_FFT_INDEX(i,j,slice_index_reducedLC)]/H + 1.0);

                        // Calculate the brightness temperature, using the optical depth
                        if(gradient_component < FRACT_FLOAT_ERR) {
                            // Gradient component goes to zero, optical depth diverges. But, since we take exp(-tau), this goes to zero and (1 - exp(-tau)) goes to unity.
                            // Again, factors of 1000. are conversions from K to mK
                            delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)] = 1000.*(Ts_z[HII_R_INDEX(i,j,slice_index_reducedLC)] - T_rad)/(1. + REDSHIFT_SAMPLE);
                        }
                        else {
                            delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)] = (1. - exp(- delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/gradient_component ))*1000.*(Ts_z[HII_R_INDEX(i,j,slice_index_reducedLC)] - T_rad)/(1. + REDSHIFT_SAMPLE);
                        }
                    }
                }
            }

            // normalised units of cell length. 0 equals beginning of cell, 1 equals end of cell
            // These are the sub-cell central positions (x_pos_offset), and the corresponding normalised value (x_pos) between 0 and 1
            for(ii=0;ii<N_RSD_STEPS;ii++) {
                x_pos_offset[ii] = subcell_width*(float)ii + subcell_width/2.;
                x_pos[ii] = x_pos_offset[ii]/( BOX_LEN/(float)HII_DIM );
            }

            // Note to convert the velocity v, to a displacement in redshift space, convert from s -> r + (1+z)*v/H(z)
            // To convert the velocity within the array v to km/s, it is a*dD/dt*delta. Where the scale factor a comes from the continuity equation
            // The array v as defined in 21cmFAST is (ik/k^2)*dD/dt*delta, as it is defined as a comoving quantity (scale factor is implicit).
            // However, the conversion between real and redshift space also picks up a scale factor, therefore the scale factors drop out and therefore
            // the displacement of the sub-cells is purely determined from the array, v and the Hubble factor: v/H.

            for (i=0; i<HII_DIM; i++){
                for (j=0; j<HII_DIM; j++){

                    // Generate the optical-depth for the specific line-of-sight with R.S.D
                    for(k=0;k<HII_DIM;k++) {
                        delta_T_RSD_LOS[k] = 0.0;
                    }

                    for (k=0; k<LOOP_INDEX; k++){

                        if(!USE_FCOLL_IONISATION_TABLE) {
                            slice_index_reducedLC = k;
                        }
                        else {
                            if((min_slice_index + k) >= HII_DIM) {
                                slice_index_reducedLC = (min_slice_index + k) - HII_DIM;
                            }
                            else {
                                slice_index_reducedLC = (min_slice_index + k);
                            }
                        }

                        if((fabs(delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]) >= FRACT_FLOAT_ERR) && (xH[HII_R_INDEX(i,j,slice_index_reducedLC)] >= FRACT_FLOAT_ERR)) {

                            if(slice_index_reducedLC==0) {
                                d1_low = v[HII_R_FFT_INDEX(i,j,HII_DIM-1)]/H;
                                d2_low = v[HII_R_FFT_INDEX(i,j,slice_index_reducedLC)]/H;
                            }
                            else {
                                d1_low = v[HII_R_FFT_INDEX(i,j,slice_index_reducedLC-1)]/H;
                                d2_low = v[HII_R_FFT_INDEX(i,j,slice_index_reducedLC)]/H;
                            }
                            // Displacements (converted from velocity) for the original cell centres straddling half of the sub-cells (cell after)
                            if(slice_index_reducedLC==(HII_DIM-1)) {
                                d1_high = v[HII_R_FFT_INDEX(i,j,slice_index_reducedLC)]/H;
                                d2_high = v[HII_R_FFT_INDEX(i,j,0)]/H;
                            }
                            else {
                                d1_high = v[HII_R_FFT_INDEX(i,j,slice_index_reducedLC)]/H;
                                d2_high = v[HII_R_FFT_INDEX(i,j,slice_index_reducedLC+1)]/H;
                            }

                            for(ii=0;ii<N_RSD_STEPS;ii++) {

                                // linearly interpolate the displacements to determine the corresponding displacements of the sub-cells
                                // Checking of 0.5 is for determining if we are left or right of the mid-point of the original cell (for the linear interpolation of the displacement)
                                // to use the appropriate cell

                                if(x_pos[ii] <= 0.5) {
                                    subcell_displacement = d1_low + ( (x_pos[ii] + 0.5 ) - x_val1)*( d2_low - d1_low )/( x_val2 - x_val1 );
                                }
                                else {
                                    subcell_displacement = d1_high + ( (x_pos[ii] - 0.5 ) - x_val1)*( d2_high - d1_high )/( x_val2 - x_val1 );
                                }

                                // The new centre of the sub-cell post R.S.D displacement. Normalised to units of cell width for determining it's displacement
                                RSD_pos_new = (x_pos_offset[ii] + subcell_displacement)/( BOX_LEN/(float)HII_DIM );
                                // The sub-cell boundaries of the sub-cell, for determining the fractional contribution of the sub-cell to neighbouring cells when
                                // the sub-cell straddles two cell positions
                                RSD_pos_new_boundary_low = RSD_pos_new - (subcell_width/2.)/( BOX_LEN/(float)HII_DIM );
                                RSD_pos_new_boundary_high = RSD_pos_new + (subcell_width/2.)/( BOX_LEN/(float)HII_DIM );

                                if(RSD_pos_new_boundary_low >= 0.0 && RSD_pos_new_boundary_high < 1.0) {
                                    // sub-cell has remained in the original cell (just add it back to the original cell)

                                    delta_T_RSD_LOS[slice_index_reducedLC] += delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                }
                                else if(RSD_pos_new_boundary_low < 0.0 && RSD_pos_new_boundary_high < 0.0) {
                                    // sub-cell has moved completely into a new cell (toward the observer)

                                    // determine how far the sub-cell has moved in units of original cell boundary
                                    cell_distance = ceil(fabs(RSD_pos_new_boundary_low))-1.;

                                    // Determine the location of the sub-cell relative to the original cell binning
                                    if(fabs(RSD_pos_new_boundary_high) > cell_distance) {
                                        // sub-cell is entirely contained within the new cell (just add it to the new cell)

                                        // check if the new cell position is at the edge of the box. If so, periodic boundary conditions
                                        if(slice_index_reducedLC<((int)cell_distance+1)) {
                                            delta_T_RSD_LOS[slice_index_reducedLC-((int)cell_distance+1) + HII_DIM] += delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        }
                                        else {
                                            delta_T_RSD_LOS[slice_index_reducedLC-((int)cell_distance+1)] += delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        }
                                    }
                                    else {
                                        // sub-cell is partially contained within the cell

                                        // Determine the fraction of the sub-cell which is in either of the two original cells
                                        fraction_outside = (fabs(RSD_pos_new_boundary_low) - cell_distance)/(subcell_width/( BOX_LEN/(float)HII_DIM ));
                                        fraction_within = 1. - fraction_outside;

                                        // Check if the first part of the sub-cell is at the box edge
                                        if(slice_index_reducedLC<(((int)cell_distance))) {
                                            delta_T_RSD_LOS[slice_index_reducedLC-((int)cell_distance) + HII_DIM] += fraction_within*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        }
                                        else {
                                            delta_T_RSD_LOS[slice_index_reducedLC-((int)cell_distance)] += fraction_within*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        }
                                        // Check if the second part of the sub-cell is at the box edge
                                        if(slice_index_reducedLC<(((int)cell_distance + 1))) {
                                            delta_T_RSD_LOS[slice_index_reducedLC-((int)cell_distance+1) + HII_DIM] += fraction_outside*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        }
                                        else {
                                            delta_T_RSD_LOS[slice_index_reducedLC-((int)cell_distance+1)] += fraction_outside*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        }
                                    }
                                }
                                else if(RSD_pos_new_boundary_low < 0.0 && (RSD_pos_new_boundary_high > 0.0 && RSD_pos_new_boundary_high < 1.0)) {
                                    // sub-cell has moved partially into a new cell (toward the observer)

                                    // Determine the fraction of the sub-cell which is in either of the two original cells
                                    fraction_within = RSD_pos_new_boundary_high/(subcell_width/( BOX_LEN/(float)HII_DIM ));
                                    fraction_outside = 1. - fraction_within;

                                    // Check the periodic boundaries conditions and move the fraction of each sub-cell to the appropriate new cell
                                    if(slice_index_reducedLC==0) {
                                        delta_T_RSD_LOS[HII_DIM-1] += fraction_outside*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        delta_T_RSD_LOS[slice_index_reducedLC] += fraction_within*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                    }
                                    else {
                                        delta_T_RSD_LOS[slice_index_reducedLC-1] += fraction_outside*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        delta_T_RSD_LOS[slice_index_reducedLC] += fraction_within*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                    }
                                }
                                else if((RSD_pos_new_boundary_low >= 0.0 && RSD_pos_new_boundary_low < 1.0) && (RSD_pos_new_boundary_high >= 1.0)) {
                                    // sub-cell has moved partially into a new cell (away from the observer)

                                    // Determine the fraction of the sub-cell which is in either of the two original cells
                                    fraction_outside = (RSD_pos_new_boundary_high - 1.)/(subcell_width/( BOX_LEN/(float)HII_DIM ));
                                    fraction_within = 1. - fraction_outside;

                                    // Check the periodic boundaries conditions and move the fraction of each sub-cell to the appropriate new cell
                                    if(slice_index_reducedLC==(HII_DIM-1)) {
                                        delta_T_RSD_LOS[slice_index_reducedLC] += fraction_within*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        delta_T_RSD_LOS[0] += fraction_outside*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                    }
                                    else {
                                        delta_T_RSD_LOS[slice_index_reducedLC] += fraction_within*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        delta_T_RSD_LOS[slice_index_reducedLC+1] += fraction_outside*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                    }
                                }
                                else {
                                    // sub-cell has moved completely into a new cell (away from the observer)

                                    // determine how far the sub-cell has moved in units of original cell boundary
                                    cell_distance = floor(fabs(RSD_pos_new_boundary_high));

                                    if(RSD_pos_new_boundary_low >= cell_distance) {
                                        // sub-cell is entirely contained within the new cell (just add it to the new cell)

                                        // check if the new cell position is at the edge of the box. If so, periodic boundary conditions
                                        if(slice_index_reducedLC>(HII_DIM - 1 - (int)cell_distance)) {
                                            delta_T_RSD_LOS[slice_index_reducedLC+(int)cell_distance - HII_DIM] += delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        }
                                        else {
                                            delta_T_RSD_LOS[slice_index_reducedLC+(int)cell_distance] += delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        }
                                    }
                                    else {
                                        // sub-cell is partially contained within the cell

                                        // Determine the fraction of the sub-cell which is in either of the two original cells
                                        fraction_outside = (RSD_pos_new_boundary_high - cell_distance)/(subcell_width/( BOX_LEN/(float)HII_DIM ));
                                        fraction_within = 1. - fraction_outside;

                                        // Check if the first part of the sub-cell is at the box edge
                                        if(slice_index_reducedLC>(HII_DIM - 1 - ((int)cell_distance-1))) {
                                            delta_T_RSD_LOS[slice_index_reducedLC+(int)cell_distance-1 - HII_DIM] += fraction_within*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        }
                                        else {
                                            delta_T_RSD_LOS[slice_index_reducedLC+(int)cell_distance-1] += fraction_within*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        }
                                        // Check if the second part of the sub-cell is at the box edge
                                        if(slice_index_reducedLC>(HII_DIM - 1 - ((int)cell_distance))) {
                                            delta_T_RSD_LOS[slice_index_reducedLC+(int)cell_distance - HII_DIM] += fraction_outside*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        }
                                        else {
                                            delta_T_RSD_LOS[slice_index_reducedLC+(int)cell_distance] += fraction_outside*delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)]/(float)N_RSD_STEPS;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if(!USE_FCOLL_IONISATION_TABLE) {
                        for(k=0;k<HII_DIM;k++) {
                            delta_T[HII_R_INDEX(i,j,k)] = delta_T_RSD_LOS[k];

                            ave += delta_T_RSD_LOS[k];
                        }
                    }
                    else {

                        for(k=0;k<total_in_z;k++) {

                            if((min_slice_index + k) >= HII_DIM) {
                                slice_index_reducedLC = (min_slice_index + k) - HII_DIM;
                            }
                            else {
                                slice_index_reducedLC = (min_slice_index + k);
                            }
                            delta_T[HII_R_INDEX(i,j,slice_index_reducedLC)] = delta_T_RSD_LOS[slice_index_reducedLC];

                            ave += delta_T_RSD_LOS[slice_index_reducedLC];
                        }
                    }
                }
            }

            if(!USE_FCOLL_IONISATION_TABLE) {
                ave /= (float)HII_TOT_NUM_PIXELS;
            }
            else {
                ave /= ((float)HII_DIM*(float)HII_DIM*(float)total_in_z);
            }
        }
        else {

            // now add the velocity correction to the delta_T maps
            max_v_deriv = fabs(MAX_DVDR*H);

            for (i=0; i<HII_DIM; i++){
                for (j=0; j<HII_DIM; j++){
                    for (k=0; k<HII_DIM; k++){

                        dvdx = vel_gradient[HII_R_FFT_INDEX(i,j,k)];

                        // set maximum allowed gradient for this linear approximation
                        if (fabs(dvdx) > max_v_deriv){
                            if (dvdx < 0) dvdx = -max_v_deriv;
                            else dvdx = max_v_deriv;
                            //                               nonlin_ct++;
                        }

                        delta_T[HII_R_INDEX(i,j,k)] /= (dvdx/H + 1.0);

                        ave += delta_T[HII_R_INDEX(i,j,k)];
                    }
                }
            }
            ave /= (HII_TOT_NUM_PIXELS+0.0);
        }
    }

    // Note that the average brightness temperature will be less reliable if USE_FCOLL_IONISATION_TABLE is set as it uses only a sub-volume (i.e. larger sample variance)
    if(STORE_DATA) {
        aveTb[sample_index] = ave;
    }


//JBM: We can output the pow sp of velocities just to see it:
    if(PRINT_DEBUG_REL_VEL){
        smoothed_box_vcb = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
        if (!smoothed_box_vcb){
          fprintf(stderr, "drive_21.c: Error allocating memory for low-res box_vcb.\nAborting...\n");
        }

        fprintf(stderr, "Reading real-space vcb to get its power spectrum:\n");
        sprintf(filename, "%s/smoothed_vcb_x_z0.00_%i_%.0fMpc", BOXES_INPUT_FOLDER, HII_DIM, BOX_LEN);
        F = fopen(filename, "rb");
        if (!F){
          fprintf(stderr, "Couldn't open file %s for reading\nAborting...\n", filename);
        }
        if (mod_fread(smoothed_box_vcb, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
          fprintf(stderr, "drive_21.c: Read error occured!\n");
        }

        GeneratePS_vel(0,0.0);

        sprintf(filename, "%s/delta_vcb_spectrum.txt",OUTPUT_FOLDER);
        F=fopen(filename, "wt");
        for (ct=1; ct<NUM_BINS; ct++){
            if (in_bin_ct[ct]>0)
                fprintf(F, "%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
        }
        fclose(F);

        SAFE_FREE(smoothed_box_vcb);

    }//JBM:continue the rest of the code:




    /////////////////////////////  PRINT OUT THE POWERSPECTRUM  ///////////////////////////////

    if(USE_LIGHTCONE) {

        if(sample_index==0) {

            // NOTE: This only works under the assumption that the hightest redshift is >> larger than reionisation
            // i.e assumes that Z_HEAT_MAX >~ 20.

            // Set the highest redshift box for creating the light-cone
            memcpy(box_z2,delta_T,sizeof(float)*HII_TOT_NUM_PIXELS);
            z2_LC = redshifts[sample_index];
            t_z2_LC = gettime(z2_LC);
        }
        else {

            // LOS_direction is used in two separate locations (iterated), so we store the state for each individucal usage
            LOS_direction = Stored_LOS_direction_state_2;

            memcpy(box_z1,delta_T,sizeof(float)*HII_TOT_NUM_PIXELS);
            z1_LC = redshifts[sample_index];
            t_z1_LC = gettime(z1_LC);

            // Box will fill up, therefore, finish writing box, write out full box, continue onto the next box
            if(start_index_LC[N_USER_REDSHIFT-sample_index-1] > end_index_LC[N_USER_REDSHIFT-sample_index-1]) {

                for(ii=0;ii<end_index_LC[N_USER_REDSHIFT-sample_index-1];ii++) {
                    total_slice_ct = num_boxes_interp*HII_DIM + ii;

                    z_slice = slice_redshifts[full_index_LC[total_slice_ct]];
                    t_z_slice = gettime(z_slice);

                    for (j=0;j<HII_DIM; j++){
                        for (k=0; k<HII_DIM; k++){
                            fz1 = box_z1[coeval_box_pos(LOS_direction,j,k,ii)];
                            fz2 = box_z2[coeval_box_pos(LOS_direction,j,k,ii)];
                            fz = (fz2 - fz1) / (t_z2_LC - t_z1_LC) * (t_z_slice - t_z1_LC) + fz1; // linearly interpolate in z (time actually)
                            if(num_boxes_interp==total_num_boxes) {
                                box_interpolate_remainder[ii + remainder_LC*(k+HII_DIM*j)] = fz;
                            }
                            else {
                                box_interpolate[HII_R_INDEX(j, k, ii)] = fz;
                            }
                        }
                    }
                }

                if(num_boxes_interp==total_num_boxes) {

                    if(PRINT_LIGHTCONE_21cmBoxes) {
                        sprintf(filename, "%s/delta_T_%f_%f__zstart%09.5f_zend%09.5f_%i_%.0fMpc_lighttravel",OUTPUT_FOLDER,INDIVIDUAL_ID,INDIVIDUAL_ID_2, slice_redshifts[full_index_LC[num_boxes_interp*HII_DIM]], redshifts[0], HII_DIM, BOX_LEN);
                        F = fopen(filename, "wb");
                        mod_fwrite(box_interpolate_remainder, sizeof(float)*(remainder_LC*HII_DIM*HII_DIM), 1, F);
                        fclose(F);
                    }
                    ave = 0.0;
                }
                else {

                    if(PRINT_LIGHTCONE_21cmBoxes) {
                        sprintf(filename, "%s/delta_T_%f_%f__zstart%09.5f_zend%09.5f_%i_%.0fMpc_lighttravel",OUTPUT_FOLDER, INDIVIDUAL_ID,INDIVIDUAL_ID_2, slice_redshifts[full_index_LC[num_boxes_interp*HII_DIM]], slice_redshifts[full_index_LC[(num_boxes_interp+1)*HII_DIM]], HII_DIM, BOX_LEN);
                        F = fopen(filename, "wb");
                        mod_fwrite(box_interpolate, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F);
                        fclose(F);

                        //JBM:we also save in regular text

                        sprintf(filename, "%s/delta_T_%f_%f__zstart%09.5f_zend%09.5f_%i_%.0fMpc_lighttravel.txt",OUTPUT_FOLDER, INDIVIDUAL_ID,INDIVIDUAL_ID_2, slice_redshifts[full_index_LC[num_boxes_interp*HII_DIM]], slice_redshifts[full_index_LC[(num_boxes_interp+1)*HII_DIM]], HII_DIM, BOX_LEN);
                        F = fopen(filename, "w");

                        for(ii=0;ii<end_index_LC[N_USER_REDSHIFT-sample_index-1];ii++) {
                            for (j=0;j<HII_DIM; j++){
                                for (k=0; k<HII_DIM; k++){
                                    fprintf(F, "%le %le %le \n", j* BOX_LEN/(HII_DIM-1.0), k* BOX_LEN/(HII_DIM-1.0), box_interpolate[HII_R_INDEX(j, k, ii)]);
                                    }
                                }
                            }

                        fclose(F);
                    }




                    // We have filled the light-cone cube box, now generate the 21cm PS.
                    GeneratePS(0,0.0);

                    // now lets print out the k bins
                    if(PRINT_FILES) {

                        sprintf(filename, "%s/delTps_estimate_%f_%f_zstart%09.5f_zend%09.5f_%i_%.0fMpc_lighttravel.txt",OUTPUT_FOLDER,INDIVIDUAL_ID,INDIVIDUAL_ID_2,slice_redshifts[full_index_LC[num_boxes_interp*HII_DIM]], slice_redshifts[full_index_LC[(num_boxes_interp+1)*HII_DIM]], HII_DIM, BOX_LEN);
                        if(OVERWRITE_OUTPUT || !file_exists(filename)){
                            fprintf(stderr, "Writing output(1) to %s\n", filename);
                            F=fopen(filename, "wt");
                            for (ct=1; ct<NUM_BINS; ct++){
                                if (in_bin_ct[ct]>0)
                                    fprintf(F, "%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
                            }
                            fclose(F);
                            fprintf(stderr, "Finished writing\n");
                        }
                        strcpy(lightcone_box_names[num_boxes_interp], filename);
                    }
                }

                num_boxes_interp--;

                // Storing the LOS direction for creating the light-cones (reduced calculation flag uses the same code, so double counts changes to LOS_direction. Hence, we store the state
                Stored_LOS_direction_state_2 = LOS_direction;

                for(ii=start_index_LC[N_USER_REDSHIFT-sample_index-1];ii<HII_DIM;ii++) {
                    // (num_boxes_interp-1) as we are keeping track of the boxes in descending order (num_boxes_interp is decreased by one, further down)
                    // Yes, I know this is the counter intuitive way of doing it, but it renders Andrei's original code mostly intact
                    total_slice_ct = num_boxes_interp*HII_DIM + ii;

                    z_slice = slice_redshifts[full_index_LC[total_slice_ct]];
                    t_z_slice = gettime(z_slice);

                    for (j=0;j<HII_DIM; j++){
                        for (k=0; k<HII_DIM; k++){
                            fz1 = box_z1[coeval_box_pos(LOS_direction,j,k,ii)];
                            fz2 = box_z2[coeval_box_pos(LOS_direction,j,k,ii)];
                            fz = (fz2 - fz1) / (t_z2_LC - t_z1_LC) * (t_z_slice - t_z1_LC) + fz1; // linearly interpolate in z (time actually)
                            if(num_boxes_interp==total_num_boxes) {
                                box_interpolate_remainder[ii + remainder_LC*(k+HII_DIM*j)] = fz;
                            }
                            else {
                                box_interpolate[HII_R_INDEX(j, k, ii)] = fz;
                            }
                        }
                    }
                }
            }
            else{

                for(ii=start_index_LC[N_USER_REDSHIFT-sample_index-1];ii<end_index_LC[N_USER_REDSHIFT-sample_index-1];ii++) {
                    total_slice_ct = num_boxes_interp*HII_DIM + ii;

                    z_slice = slice_redshifts[full_index_LC[total_slice_ct]];
                    t_z_slice = gettime(z_slice);

                    for (j=0;j<HII_DIM; j++){
                        for (k=0; k<HII_DIM; k++){
                            fz1 = box_z1[coeval_box_pos(LOS_direction,j,k,ii)];
                            fz2 = box_z2[coeval_box_pos(LOS_direction,j,k,ii)];
                            fz = (fz2 - fz1) / (t_z2_LC - t_z1_LC) * (t_z_slice - t_z1_LC) + fz1; // linearly interpolate in z (time actually)
                            if(num_boxes_interp==total_num_boxes) {
                                box_interpolate_remainder[ii + remainder_LC*(k+HII_DIM*j)] = fz;
                            }
                            else {
                                box_interpolate[HII_R_INDEX(j, k, ii)] = fz;
                            }
                        }
                    }
                }
            }
            // ok, we are done with this box, move onto the next box
            memcpy(box_z2,delta_T,sizeof(float)*HII_TOT_NUM_PIXELS);
            z2_LC = z1_LC;
            t_z2_LC = gettime(z2_LC);

            if((num_boxes_interp==0)&&(sample_index==(N_USER_REDSHIFT-1))) {

                if(PRINT_LIGHTCONE_21cmBoxes) {
                    sprintf(filename, "%s/delta_T_%f_%f__zstart%09.5f_zend%09.5f_%i_%.0fMpc_lighttravel",OUTPUT_FOLDER,INDIVIDUAL_ID,INDIVIDUAL_ID_2, slice_redshifts[full_index_LC[num_boxes_interp*HII_DIM]], slice_redshifts[full_index_LC[(num_boxes_interp+1)*HII_DIM]], HII_DIM, BOX_LEN);
                    F = fopen(filename, "wb");
                    mod_fwrite(box_interpolate, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F);
                    fclose(F);
                }

                // We have filled the light-cone cube box, now generate the 21cm PS.
                GeneratePS(0,0.0);

                // now lets print out the k bins
                if(PRINT_FILES) {
                    sprintf(filename, "%s/delTps_estimate_%f_%f_zstart%09.5f_zend%09.5f_%i_%.0fMpc_lighttravel.txt",OUTPUT_FOLDER,INDIVIDUAL_ID,INDIVIDUAL_ID_2,slice_redshifts[full_index_LC[num_boxes_interp*HII_DIM]], slice_redshifts[full_index_LC[(num_boxes_interp+1)*HII_DIM]], HII_DIM, BOX_LEN);
                    if(OVERWRITE_OUTPUT || !file_exists(filename)){
                        fprintf(stderr, "Writing output(2) to %s\n", filename);
                        F=fopen(filename, "wt");
                        for (ct=1; ct<NUM_BINS; ct++){
                            if (in_bin_ct[ct]>0)
                                fprintf(F, "%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
                        }
                        fclose(F);
                        fprintf(stderr, "Finished writing\n");
                    }
                    strcpy(lightcone_box_names[num_boxes_interp], filename);
                }
            }
        }
    }

    if(!USE_LIGHTCONE) {

        if(PRINT_COEVAL_21cmBoxes) {
            sprintf(filename, "%s/delta_T_%f_%f_z%f1.6_%i_%.0fMpc",OUTPUT_FOLDER,INDIVIDUAL_ID,INDIVIDUAL_ID_2, REDSHIFT_SAMPLE, HII_DIM, BOX_LEN);
            F = fopen(filename, "wb");
            mod_fwrite(delta_T, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F);
            fclose(F);

        	//  //JBM:we can also save in regular text (readable)
            // sprintf(filename, "%s/readable_delta_T_%f_%f_z%f1.6_%i_%.0fMpc.txt",OUTPUT_FOLDER,INDIVIDUAL_ID,INDIVIDUAL_ID_2, REDSHIFT_SAMPLE, HII_DIM, BOX_LEN);
            // F = fopen(filename, "w");
            //
            // for(i=0;i<HII_DIM;i++) {
            //     for (j=0;j<HII_DIM; j++){
            //         for (k=0; k<HII_DIM; k++){
            //             fprintf(F, "%le %le %le %le \n", i* BOX_LEN/(HII_DIM-1.0), j* BOX_LEN/(HII_DIM-1.0), k* BOX_LEN/(HII_DIM-1.0), delta_T[HII_R_INDEX(i,j, k)]);
            //             }
            //         }
            //     }
            //
            // fclose(F);

        }

        GeneratePS(1,ave);

        // now lets print out the k bins
        if(PRINT_FILES) {
            sprintf(filename, "%s/delTps_estimate_%f_%f_%f.txt",OUTPUT_FOLDER,INDIVIDUAL_ID,INDIVIDUAL_ID_2,REDSHIFT_SAMPLE);
            if(OVERWRITE_OUTPUT || !file_exists(filename)){
                fprintf(stderr, "Writing output(3) to %s\n", filename);
                F=fopen(filename, "wt");
                for (ct=1; ct<NUM_BINS; ct++){
                    if (in_bin_ct[ct]>0)
                        fprintf(F, "%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
                }
                fclose(F);
                fprintf(stderr, "Finished writing\n");
            }
        }


        // now lets print out the k bins
        if(PRINT_FILES && DO_ANISOTROPY) {
            sprintf(filename, "%s/delTps_estimate_ani_%f_%f_%f.txt",OUTPUT_FOLDER,INDIVIDUAL_ID,INDIVIDUAL_ID_2,REDSHIFT_SAMPLE);
            GeneratePS_aniso(1, ave, filename);
//             F=fopen(filename, "wt");
//             for (ct=1; ct<NUM_BINS; ct++){
//                 if (in_bin_ct[ct]>0)
//                     fprintf(F, "%e \t %e \t %e \n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
//             }
//             fclose(F);
        }


    ///////////////////////////// END POWER SPECTRUM STUFF   ////////////////////////////////////

    //////////////////////////// End of perform 'delta_T.c' /////////////////////////////////////
    }

    SAFE_FREE(LOS_index);
    SAFE_FREE(slice_index);

    destroy_21cmMC_HII_arrays(skip_deallocate);
}

void ComputeInitialConditions() {

    /*
     Generates the initial conditions: gaussian random density field (DIM^3) as well as the equal or lower resolution velocity fields, and smoothed density field (HII_DIM^3).
     See INIT_PARAMS.H and ANAL_PARAMS.H to set the appropriate parameters.
     Output is written to ../Boxes

     Author: Andrei Mesinger
     Date: 9/29/06
     */

     /*
     Output is now written to the BOXES_INPUT_FOLDER defined in COSMOLOGY.H

     Thomas Smeman
     Date: 13/05/2024
     */

    fftwf_plan plan;

    unsigned long long ct;
    int n_x, n_y, n_z, i, j, k, ii;
    float k_x, k_y, k_z, k_mag, p, a, b, k_sq;
    double pixel_deltax;

    float f_pixel_factor;

    gsl_rng * r;

    /************  INITIALIZATION **********************/

    // Removed all references to threads as 21CMMC is always a single core implementation

    // seed the random number generators
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, RANDOM_SEED);

    // allocate array for the k-space and real-space boxes
    HIRES_box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*KSPACE_NUM_PIXELS);
    HIRES_box_saved = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*KSPACE_NUM_PIXELS);

    // now allocate memory for the lower-resolution box
    // use HII_DIM from ANAL_PARAMS
    LOWRES_density = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    LOWRES_vx = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    LOWRES_vy= (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    LOWRES_vz = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);

    if(SECOND_ORDER_LPT_CORRECTIONS){
        LOWRES_vx_2LPT = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
        LOWRES_vy_2LPT = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
        LOWRES_vz_2LPT = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
    }

    // find factor of HII pixel size / deltax pixel size
    f_pixel_factor = DIM/(float)HII_DIM;
    /************  END INITIALIZATION ******************/

    /************ CREATE K-SPACE GAUSSIAN RANDOM FIELD ***********/
    for (n_x=0; n_x<DIM; n_x++){
        // convert index to numerical value for this component of the k-mode: k = (2*pi/L) * n
        if (n_x>MIDDLE)
            k_x =(n_x-DIM) * DELTA_K;  // wrap around for FFT convention
        else
            k_x = n_x * DELTA_K;

        for (n_y=0; n_y<DIM; n_y++){
            // convert index to numerical value for this component of the k-mode: k = (2*pi/L) * n
            if (n_y>MIDDLE)
                k_y =(n_y-DIM) * DELTA_K;
            else
                k_y = n_y * DELTA_K;

            // since physical space field is real, only half contains independent modes
            for (n_z=0; n_z<=MIDDLE; n_z++){
                // convert index to numerical value for this component of the k-mode: k = (2*pi/L) * n
                k_z = n_z * DELTA_K;

                // now get the power spectrum; remember, only the magnitude of k counts (due to issotropy)
                // this could be used to speed-up later maybe
                k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);
                p = power_in_k(k_mag);

                // ok, now we can draw the values of the real and imaginary part
                // of our k entry from a Gaussian distribution
                a = gsl_ran_ugaussian(r);
                b = gsl_ran_ugaussian(r);
                HIRES_box[C_INDEX(n_x, n_y, n_z)] = sqrt(VOLUME*p/2.0) * (a + b*I);
            }
        }
    }

    /*****  Adjust the complex conjugate relations for a real array  *****/
    adj_complex_conj(HIRES_box);

    /*** Let's also create a lower-resolution version of the density field  ***/

    memcpy(HIRES_box_saved, HIRES_box, sizeof(fftwf_complex)*KSPACE_NUM_PIXELS);

    if (DIM != HII_DIM)
        filter(HIRES_box, 0, L_FACTOR*BOX_LEN/(HII_DIM+0.0));
    // FFT back to real space
    plan = fftwf_plan_dft_c2r_3d(DIM, DIM, DIM, (fftwf_complex *)HIRES_box, (float *)HIRES_box, FFTW_ESTIMATE);
    fftwf_execute(plan);
    // now sample the filtered box
    for (i=0; i<HII_DIM; i++){
        for (j=0; j<HII_DIM; j++){
            for (k=0; k<HII_DIM; k++){
                LOWRES_density[HII_R_INDEX(i,j,k)] =
                *((float *)HIRES_box + R_FFT_INDEX((unsigned long long)(i*f_pixel_factor+0.5),
                                             (unsigned long long)(j*f_pixel_factor+0.5),
                                             (unsigned long long)(k*f_pixel_factor+0.5)))/VOLUME;
            }
        }
    }
    /******* PERFORM INVERSE FOURIER TRANSFORM *****************/
    // add the 1/VOLUME factor when converting from k space to real space

    memcpy(HIRES_box, HIRES_box_saved, sizeof(fftwf_complex)*KSPACE_NUM_PIXELS);

    for (ct=0; ct<KSPACE_NUM_PIXELS; ct++){
        HIRES_box[ct] /= VOLUME;
    }
    plan = fftwf_plan_dft_c2r_3d(DIM, DIM, DIM, (fftwf_complex *)HIRES_box, (float *)HIRES_box, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();

    for (i=0; i<DIM; i++){
        for (j=0; j<DIM; j++){
            for (k=0; k<DIM; k++){
                *((float *)HIRES_density + R_FFT_INDEX(i,j,k)) = *((float *)HIRES_box + R_FFT_INDEX(i,j,k));
            }
        }
    }

    for(ii=0;ii<3;ii++) {

        memcpy(HIRES_box, HIRES_box_saved, sizeof(fftwf_complex)*KSPACE_NUM_PIXELS);
        /*** Now let's set the velocity field/dD/dt (in comoving Mpc) ***/

        for (n_x=0; n_x<DIM; n_x++){
            if (n_x>MIDDLE)
                k_x =(n_x-DIM) * DELTA_K;  // wrap around for FFT convention
            else
                k_x = n_x * DELTA_K;

            for (n_y=0; n_y<DIM; n_y++){
                if (n_y>MIDDLE)
                    k_y =(n_y-DIM) * DELTA_K;
                else
                    k_y = n_y * DELTA_K;

                for (n_z=0; n_z<=MIDDLE; n_z++){
                    k_z = n_z * DELTA_K;

                    k_sq = k_x*k_x + k_y*k_y + k_z*k_z;

                    // now set the velocities
                    if ((n_x==0) && (n_y==0) && (n_z==0)){ // DC mode
                        HIRES_box[0] = 0;
                    }
                    else{
                        if(ii==0) {
                            HIRES_box[C_INDEX(n_x,n_y,n_z)] *= k_x*I/k_sq/VOLUME;
                        }
                        if(ii==1) {
                            HIRES_box[C_INDEX(n_x,n_y,n_z)] *= k_y*I/k_sq/VOLUME;
                        }
                        if(ii==2) {
                            HIRES_box[C_INDEX(n_x,n_y,n_z)] *= k_z*I/k_sq/VOLUME;
                        }
                        // note the last factor of 1/VOLUME accounts for the scaling in real-space, following the FFT
                    }
                }
            }
        }

        if (DIM != HII_DIM)
            filter(HIRES_box, 0, L_FACTOR*BOX_LEN/(HII_DIM+0.0));

        plan = fftwf_plan_dft_c2r_3d(DIM, DIM, DIM, (fftwf_complex *)HIRES_box, (float *)HIRES_box, FFTW_ESTIMATE);
        fftwf_execute(plan);
        // now sample to lower res
        // now sample the filtered box
        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    if(ii==0) {
                        LOWRES_vx[HII_R_INDEX(i,j,k)] =
                        *((float *)HIRES_box + R_FFT_INDEX((unsigned long long)(i*f_pixel_factor+0.5),
                                                           (unsigned long long)(j*f_pixel_factor+0.5),
                                                           (unsigned long long)(k*f_pixel_factor+0.5)));
                    }
                    if(ii==1) {
                        LOWRES_vy[HII_R_INDEX(i,j,k)] =
                        *((float *)HIRES_box + R_FFT_INDEX((unsigned long long)(i*f_pixel_factor+0.5),
                                                           (unsigned long long)(j*f_pixel_factor+0.5),
                                                           (unsigned long long)(k*f_pixel_factor+0.5)));
                    }
                    if(ii==2) {
                        LOWRES_vz[HII_R_INDEX(i,j,k)] =
                        *((float *)HIRES_box + R_FFT_INDEX((unsigned long long)(i*f_pixel_factor+0.5),
                                                           (unsigned long long)(j*f_pixel_factor+0.5),
                                                           (unsigned long long)(k*f_pixel_factor+0.5)));
                    }
                }
            }
        }
    }
    // write out file

    /* *************************************************** *
     *              BEGIN 2LPT PART                        *
     * *************************************************** */

    // Generation of the second order Lagrangian perturbation theory (2LPT) corrections to the ZA
    // reference: Scoccimarro R., 1998, MNRAS, 299, 1097-1118 Appendix D

    // Parameter set in ANAL_PARAMS.H
    if(SECOND_ORDER_LPT_CORRECTIONS){
        // use six supplementary boxes to store the gradients of phi_1 (eq. D13b)
        // Allocating the boxes
#define PHI_INDEX(i, j) ((int) ((i) - (j)) + 3*((j)) - ((int)(j))/2  )
        // ij -> INDEX
        // 00 -> 0
        // 11 -> 3
        // 22 -> 5
        // 10 -> 1
        // 20 -> 2
        // 21 -> 4

        fftwf_complex *phi_1[6];

        for(i = 0; i < 3; ++i){
            for(j = 0; j <= i; ++j){
                phi_1[PHI_INDEX(i, j)] = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*KSPACE_NUM_PIXELS);
            }
        }

        for(i = 0; i < 3; ++i){
            for(j = 0; j <= i; ++j){

                // read in the box
                memcpy(HIRES_box, HIRES_box_saved, sizeof(fftwf_complex)*KSPACE_NUM_PIXELS);

                // generate the phi_1 boxes in Fourier transform
                for (n_x=0; n_x<DIM; n_x++){
                    if (n_x>MIDDLE)
                        k_x =(n_x-DIM) * DELTA_K;  // wrap around for FFT convention
                    else
                        k_x = n_x * DELTA_K;

                    for (n_y=0; n_y<DIM; n_y++){
                        if (n_y>MIDDLE)
                            k_y =(n_y-DIM) * DELTA_K;
                        else
                            k_y = n_y * DELTA_K;

                        for (n_z=0; n_z<=MIDDLE; n_z++){
                            k_z = n_z * DELTA_K;

                            k_sq = k_x*k_x + k_y*k_y + k_z*k_z;

                            float k[] = {k_x, k_y, k_z};
                            // now set the velocities
                            if ((n_x==0) && (n_y==0) && (n_z==0)){ // DC mode
                                phi_1[PHI_INDEX(i, j)][0] = 0;
                            }
                            else{
                                phi_1[PHI_INDEX(i, j)][C_INDEX(n_x,n_y,n_z)] = -k[i]*k[j]*HIRES_box[C_INDEX(n_x, n_y, n_z)]/k_sq/VOLUME;
                                // note the last factor of 1/VOLUME accounts for the scaling in real-space, following the FFT
                            }
                        }
                    }
                }
                // Now we can generate the real phi_1[i,j]
                plan = fftwf_plan_dft_c2r_3d(DIM, DIM, DIM, (fftwf_complex *)phi_1[PHI_INDEX(i, j)], (float *)phi_1[PHI_INDEX(i, j)], FFTW_ESTIMATE);
                fftwf_execute(plan);
            }
        }

        // Then we will have the laplacian of phi_2 (eq. D13b)
        // After that we have to return in Fourier space and generate the Fourier transform of phi_2
        int m, l;
        for (i=0; i<DIM; i++){
            for (j=0; j<DIM; j++){
                for (k=0; k<DIM; k++){
                    *( (float *)HIRES_box + R_FFT_INDEX((unsigned long long)(i), (unsigned long long)(j), (unsigned long long)(k) )) = 0.0;
                    for(m = 0; m < 3; ++m){
                        for(l = m+1; l < 3; ++l){
                            *((float *)HIRES_box + R_FFT_INDEX((unsigned long long)(i),(unsigned long long)(j),(unsigned long long)(k)) ) += ( *((float *)(phi_1[PHI_INDEX(l, l)]) + R_FFT_INDEX((unsigned long long) (i),(unsigned long long) (j),(unsigned long long) (k)))  ) * (  *((float *)(phi_1[PHI_INDEX(m, m)]) + R_FFT_INDEX((unsigned long long)(i),(unsigned long long)(j),(unsigned long long)(k)))  );
                            *((float *)HIRES_box + R_FFT_INDEX((unsigned long long)(i),(unsigned long long)(j),(unsigned long long)(k)) ) -= ( *((float *)(phi_1[PHI_INDEX(l, m)]) + R_FFT_INDEX((unsigned long long)(i),(unsigned long long) (j),(unsigned long long)(k) ) )  ) * (  *((float *)(phi_1[PHI_INDEX(l, m)]) + R_FFT_INDEX((unsigned long long)(i),(unsigned long long)(j),(unsigned long long)(k) ))  );
                            *((float *)HIRES_box + R_FFT_INDEX((unsigned long long)(i),(unsigned long long)(j),(unsigned long long)(k)) ) /= TOT_NUM_PIXELS;
                        }
                    }
                }
            }
        }

        plan = fftwf_plan_dft_r2c_3d(DIM, DIM, DIM, (float *)HIRES_box, (fftwf_complex *)HIRES_box, FFTW_ESTIMATE);
        fftwf_execute(plan);

        memcpy(HIRES_box_saved, HIRES_box, sizeof(fftwf_complex)*KSPACE_NUM_PIXELS);

        // Now we can store the content of box in a back-up file
        // Then we can generate the gradients of phi_2 (eq. D13b and D9)

        /***** Write out back-up k-box RHS eq. D13b *****/

        // For each component, we generate the velocity field (same as the ZA part)

        /*** Now let's set the velocity field/dD/dt (in comoving Mpc) ***/

        // read in the box
        // TODO correct free of phi_1

        for(ii=0;ii<3;ii++) {

            if(ii>0) {
                memcpy(HIRES_box, HIRES_box_saved, sizeof(fftwf_complex)*KSPACE_NUM_PIXELS);
            }
            // set velocities/dD/dt
            for (n_x=0; n_x<DIM; n_x++){
                if (n_x>MIDDLE)
                    k_x =(n_x-DIM) * DELTA_K;  // wrap around for FFT convention
                else
                    k_x = n_x * DELTA_K;

                for (n_y=0; n_y<DIM; n_y++){
                    if (n_y>MIDDLE)
                        k_y =(n_y-DIM) * DELTA_K;
                    else
                        k_y = n_y * DELTA_K;

                    for (n_z=0; n_z<=MIDDLE; n_z++){
                        k_z = n_z * DELTA_K;

                        k_sq = k_x*k_x + k_y*k_y + k_z*k_z;

                        // now set the velocities
                        if ((n_x==0) && (n_y==0) && (n_z==0)){ // DC mode
                            HIRES_box[0] = 0;
                        }
                        else{
                            if(ii==0) {
                                HIRES_box[C_INDEX(n_x,n_y,n_z)] *= k_x*I/k_sq;
                            }
                            if(ii==1) {
                                HIRES_box[C_INDEX(n_x,n_y,n_z)] *= k_y*I/k_sq;
                            }
                            if(ii==2) {
                                HIRES_box[C_INDEX(n_x,n_y,n_z)] *= k_z*I/k_sq;
                            }
                            // note the last factor of 1/VOLUME accounts for the scaling in real-space, following the FFT
                        }
                    }
                }
            }

            if (DIM != HII_DIM)
                filter(HIRES_box, 0, L_FACTOR*BOX_LEN/(HII_DIM+0.0));

            plan = fftwf_plan_dft_c2r_3d(DIM, DIM, DIM, (fftwf_complex *)HIRES_box, (float *)HIRES_box, FFTW_ESTIMATE);
            fftwf_execute(plan);
            // now sample to lower res
            // now sample the filtered box
            for (i=0; i<HII_DIM; i++){
                for (j=0; j<HII_DIM; j++){
                    for (k=0; k<HII_DIM; k++){
                        if(ii==0) {
                            LOWRES_vx_2LPT[HII_R_INDEX(i,j,k)] =
                            *((float *)HIRES_box + R_FFT_INDEX((unsigned long long)(i*f_pixel_factor+0.5),
                                                               (unsigned long long)(j*f_pixel_factor+0.5),
                                                               (unsigned long long)(k*f_pixel_factor+0.5)));
                        }
                        if(ii==1) {
                            LOWRES_vy_2LPT[HII_R_INDEX(i,j,k)] =
                            *((float *)HIRES_box + R_FFT_INDEX((unsigned long long)(i*f_pixel_factor+0.5),
                                                               (unsigned long long)(j*f_pixel_factor+0.5),
                                                               (unsigned long long)(k*f_pixel_factor+0.5)));
                        }
                        if(ii==2) {
                            LOWRES_vz_2LPT[HII_R_INDEX(i,j,k)] =
                            *((float *)HIRES_box + R_FFT_INDEX((unsigned long long)(i*f_pixel_factor+0.5),
                                                               (unsigned long long)(j*f_pixel_factor+0.5),
                                                               (unsigned long long)(k*f_pixel_factor+0.5)));
                        }
                    }
                }
            }
        }
        // deallocate the supplementary boxes
        for(i = 0; i < 3; ++i){
            for(j = 0; j <= i; ++j){
                SAFE_FFTWF_FREE(phi_1[PHI_INDEX(i,j)]);
            }
        }
    }
    /* *********************************************** *
     *               END 2LPT PART                     *
     * *********************************************** */

    // deallocate
    SAFE_FFTWF_FREE(HIRES_box);
    SAFE_FFTWF_FREE(HIRES_box_saved);
}

/*****  Adjust the complex conjugate relations for a real array  *****/
void adj_complex_conj(){
    int i, j, k;

    // corners
    HIRES_box[C_INDEX(0,0,0)] = 0;
    HIRES_box[C_INDEX(0,0,MIDDLE)] = crealf(HIRES_box[C_INDEX(0,0,MIDDLE)]);
    HIRES_box[C_INDEX(0,MIDDLE,0)] = crealf(HIRES_box[C_INDEX(0,MIDDLE,0)]);
    HIRES_box[C_INDEX(0,MIDDLE,MIDDLE)] = crealf(HIRES_box[C_INDEX(0,MIDDLE,MIDDLE)]);
    HIRES_box[C_INDEX(MIDDLE,0,0)] = crealf(HIRES_box[C_INDEX(MIDDLE,0,0)]);
    HIRES_box[C_INDEX(MIDDLE,0,MIDDLE)] = crealf(HIRES_box[C_INDEX(MIDDLE,0,MIDDLE)]);
    HIRES_box[C_INDEX(MIDDLE,MIDDLE,0)] = crealf(HIRES_box[C_INDEX(MIDDLE,MIDDLE,0)]);
    HIRES_box[C_INDEX(MIDDLE,MIDDLE,MIDDLE)] = crealf(HIRES_box[C_INDEX(MIDDLE,MIDDLE,MIDDLE)]);

    // do entire i except corners
    for (i=1; i<MIDDLE; i++){
        // just j corners
        for (j=0; j<=MIDDLE; j+=MIDDLE){
            for (k=0; k<=MIDDLE; k+=MIDDLE){
                HIRES_box[C_INDEX(i,j,k)] = conjf(HIRES_box[C_INDEX(DIM-i,j,k)]);
            }
        }

        // all of j
        for (j=1; j<MIDDLE; j++){
            for (k=0; k<=MIDDLE; k+=MIDDLE){
                HIRES_box[C_INDEX(i,j,k)] = conjf(HIRES_box[C_INDEX(DIM-i,DIM-j,k)]);
                HIRES_box[C_INDEX(i,DIM-j,k)] = conjf(HIRES_box[C_INDEX(DIM-i,j,k)]);
            }
        }
    } // end loop over i

    // now the i corners
    for (i=0; i<=MIDDLE; i+=MIDDLE){
        for (j=1; j<MIDDLE; j++){
            for (k=0; k<=MIDDLE; k+=MIDDLE){
                HIRES_box[C_INDEX(i,j,k)] = conjf(HIRES_box[C_INDEX(i,DIM-j,k)]);
            }
        }
    } // end loop over remaining j
}

void ComputePerturbField(float REDSHIFT_SAMPLE) {

    /*
     USAGE: perturb_field <REDSHIFT>

     PROGRAM PERTURB_FIELD uses the first-order Langragian displacement field to move the masses in the cells of the density field.
     The high-res density field is extrapolated to some high-redshift (INITIAL_REDSHIFT in ANAL_PARAMS.H), then uses the zeldovich approximation
     to move the grid "particles" onto the lower-res grid we use for the maps.  Then we recalculate the velocity fields on the perturbed grid.
     */

    fftwf_complex *LOWRES_density_perturb, *LOWRES_density_perturb_saved;
    fftwf_plan plan;

    float REDSHIFT, growth_factor, displacement_factor_2LPT, init_growth_factor, init_displacement_factor_2LPT, xf, yf, zf;
    float mass_factor, dDdt, f_pixel_factor, velocity_displacement_factor, velocity_displacement_factor_2LPT;
    unsigned long long ct, HII_i, HII_j, HII_k;
    int i,j,k, xi, yi, zi;
    double ave_delta, new_ave_delta;
    /***************   BEGIN INITIALIZATION   **************************/

    // perform a very rudimentary check to see if we are underresolved and not using the linear approx
    if ((BOX_LEN > DIM) && !EVOLVE_DENSITY_LINEARLY){
        fprintf(stderr, "perturb_field.c: WARNING: Resolution is likely too low for accurate evolved density fields\n It Is recommended that you either increase the resolution (DIM/Box_LEN) or set the EVOLVE_DENSITY_LINEARLY flag to 1\n");
    }

    growth_factor = dicke(REDSHIFT_SAMPLE);
    displacement_factor_2LPT = -(3.0/7.0) * growth_factor*growth_factor; // 2LPT eq. D8

    dDdt = ddickedt(REDSHIFT_SAMPLE); // time derivative of the growth factor (1/s)
    init_growth_factor = dicke(INITIAL_REDSHIFT);
    init_displacement_factor_2LPT = -(3.0/7.0) * init_growth_factor*init_growth_factor; // 2LPT eq. D8

    // allocate memory for the updated density, and initialize
    LOWRES_density_perturb = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    LOWRES_density_perturb_saved = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

    // check if the linear evolution flag was set
    if (EVOLVE_DENSITY_LINEARLY){
        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    *((float *)LOWRES_density_perturb + HII_R_FFT_INDEX(i,j,k)) = growth_factor*LOWRES_density[HII_R_INDEX(i,j,k)];
                }
            }
        }
    }
    // first order Zel'Dovich perturbation
    else{

        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    *((float *)LOWRES_density_perturb + HII_R_FFT_INDEX(i,j,k)) = 0.;
                }
            }
        }

        velocity_displacement_factor = (growth_factor-init_growth_factor) / BOX_LEN;

        // now add the missing factor of D
        for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
            LOWRES_vx[ct] *= velocity_displacement_factor; // this is now comoving displacement in units of box size
            LOWRES_vy[ct] *= velocity_displacement_factor; // this is now comoving displacement in units of box size
            LOWRES_vz[ct] *= velocity_displacement_factor; // this is now comoving displacement in units of box size
        }

        // find factor of HII pixel size / deltax pixel size
        f_pixel_factor = DIM/(float)HII_DIM;
        mass_factor = pow(f_pixel_factor, 3);

        /* ************************************************************************* *
         *                           BEGIN 2LPT PART                                 *
         * ************************************************************************* */
        // reference: reference: Scoccimarro R., 1998, MNRAS, 299, 1097-1118 Appendix D
        if(SECOND_ORDER_LPT_CORRECTIONS){

            // allocate memory for the velocity boxes and read them in
            velocity_displacement_factor_2LPT = (displacement_factor_2LPT - init_displacement_factor_2LPT) / BOX_LEN;

            // now add the missing factor in eq. D9
            for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
                LOWRES_vx_2LPT[ct] *= velocity_displacement_factor_2LPT; // this is now comoving displacement in units of box size
                LOWRES_vy_2LPT[ct] *= velocity_displacement_factor_2LPT; // this is now comoving displacement in units of box size
                LOWRES_vz_2LPT[ct] *= velocity_displacement_factor_2LPT; // this is now comoving displacement in units of box size
            }
        }

        /* ************************************************************************* *
         *                            END 2LPT PART                                  *
         * ************************************************************************* */

        /************  END INITIALIZATION ****************************/

        // go through the high-res box, mapping the mass onto the low-res (updated) box
        for (i=0; i<DIM;i++){
            for (j=0; j<DIM;j++){
                for (k=0; k<DIM;k++){

                    // map indeces to locations in units of box size
                    xf = (i+0.5)/(DIM+0.0);
                    yf = (j+0.5)/(DIM+0.0);
                    zf = (k+0.5)/(DIM+0.0);

                    // update locations
                    HII_i = (unsigned long long)(i/f_pixel_factor);
                    HII_j = (unsigned long long)(j/f_pixel_factor);
                    HII_k = (unsigned long long)(k/f_pixel_factor);
                    xf += LOWRES_vx[HII_R_INDEX(HII_i, HII_j, HII_k)];
                    yf += LOWRES_vy[HII_R_INDEX(HII_i, HII_j, HII_k)];
                    zf += LOWRES_vz[HII_R_INDEX(HII_i, HII_j, HII_k)];

                    // 2LPT PART
                    // add second order corrections
                    if(SECOND_ORDER_LPT_CORRECTIONS){
                        xf -= LOWRES_vx_2LPT[HII_R_INDEX(HII_i,HII_j,HII_k)];
                        yf -= LOWRES_vy_2LPT[HII_R_INDEX(HII_i,HII_j,HII_k)];
                        zf -= LOWRES_vz_2LPT[HII_R_INDEX(HII_i,HII_j,HII_k)];
                    }

                    xf *= HII_DIM;
                    yf *= HII_DIM;
                    zf *= HII_DIM;
                    while (xf >= (float)HII_DIM){ xf -= HII_DIM;}
                    while (xf < 0){ xf += HII_DIM;}
                    while (yf >= (float)HII_DIM){ yf -= HII_DIM;}
                    while (yf < 0){ yf += HII_DIM;}
                    while (zf >= (float)HII_DIM){ zf -= HII_DIM;}
                    while (zf < 0){ zf += HII_DIM;}
                    xi = xf;
                    yi = yf;
                    zi = zf;
                    if (xi >= HII_DIM){ xi -= HII_DIM;}
                    if (xi < 0) {xi += HII_DIM;}
                    if (yi >= HII_DIM){ yi -= HII_DIM;}
                    if (yi < 0) {yi += HII_DIM;}
                    if (zi >= HII_DIM){ zi -= HII_DIM;}
                    if (zi < 0) {zi += HII_DIM;}

                    // now move the mass
                    *( (float *)LOWRES_density_perturb + HII_R_FFT_INDEX(xi, yi, zi) ) +=
                    (1 + init_growth_factor*HIRES_density[R_FFT_INDEX(i,j,k)]);
                }
            }
        }

        // renormalize to the new pixel size, and make into delta
        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    *((float *)LOWRES_density_perturb + HII_R_FFT_INDEX(i,j,k) ) /= mass_factor;
                    *((float *)LOWRES_density_perturb + HII_R_FFT_INDEX(i,j,k) ) -= 1;
                }
            }
        }

        // deallocate
        for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
            LOWRES_vx[ct] /= velocity_displacement_factor; // convert back to z = 0 quantity
            LOWRES_vy[ct] /= velocity_displacement_factor;
            LOWRES_vz[ct] /= velocity_displacement_factor;
        }

        if(SECOND_ORDER_LPT_CORRECTIONS){
            for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
                LOWRES_vx_2LPT[ct] /= velocity_displacement_factor_2LPT; // convert back to z = 0 quantity
                LOWRES_vy_2LPT[ct] /= velocity_displacement_factor_2LPT;
                LOWRES_vz_2LPT[ct] /= velocity_displacement_factor_2LPT;
            }
        }
    }

    /****  Print and convert to velocities *****/
    if (EVOLVE_DENSITY_LINEARLY){
        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    *((float *)LOWRES_density_REDSHIFT + HII_R_INDEX(i,j,k)) = *((float *)LOWRES_density_perturb + HII_R_FFT_INDEX(i,j,k));
                }
            }
        }

        // transform to k-space
        plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)LOWRES_density_perturb, (fftwf_complex *)LOWRES_density_perturb, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
        fftwf_cleanup();

        // save a copy of the k-space density field
    }
    else{
        // transform to k-space
        plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)LOWRES_density_perturb, (fftwf_complex *)LOWRES_density_perturb, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
        fftwf_cleanup();

        //smooth the field
        if (!EVOLVE_DENSITY_LINEARLY && SMOOTH_EVOLVED_DENSITY_FIELD){
            HII_filter(LOWRES_density_perturb, 2, R_smooth_density*BOX_LEN/(float)HII_DIM);
        }

        // save a copy of the k-space density field
        memcpy(LOWRES_density_perturb_saved, LOWRES_density_perturb, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

        plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)LOWRES_density_perturb, (float *)LOWRES_density_perturb, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
        fftwf_cleanup();

        // normalize after FFT
        for(i=0; i<HII_DIM; i++){
            for(j=0; j<HII_DIM; j++){
                for(k=0; k<HII_DIM; k++){
                    *((float *)LOWRES_density_perturb + HII_R_FFT_INDEX(i,j,k)) /= (float)HII_TOT_NUM_PIXELS;
                    if (*((float *)LOWRES_density_perturb + HII_R_FFT_INDEX(i,j,k)) < -1) // shouldn't happen
                        *((float *)LOWRES_density_perturb + HII_R_FFT_INDEX(i,j,k)) = -1+FRACT_FLOAT_ERR;
                }
            }
        }

        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    *((float *)LOWRES_density_REDSHIFT + HII_R_INDEX(i,j,k)) = *((float *)LOWRES_density_perturb + HII_R_FFT_INDEX(i,j,k));
                }
            }
        }
        memcpy(LOWRES_density_perturb, LOWRES_density_perturb_saved, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    }

    float k_x, k_y, k_z, k_sq, dDdt_over_D;
    int n_x, n_y, n_z;

    dDdt_over_D = dDdt/growth_factor;

    for (n_x=0; n_x<HII_DIM; n_x++){
        if (n_x>HII_MIDDLE)
            k_x =(n_x-HII_DIM) * DELTA_K;  // wrap around for FFT convention
        else
            k_x = n_x * DELTA_K;

        for (n_y=0; n_y<HII_DIM; n_y++){
            if (n_y>HII_MIDDLE)
                k_y =(n_y-HII_DIM) * DELTA_K;
            else
                k_y = n_y * DELTA_K;

            for (n_z=0; n_z<=HII_MIDDLE; n_z++){
                k_z = n_z * DELTA_K;

                k_sq = k_x*k_x + k_y*k_y + k_z*k_z;

                // now set the velocities
                if ((n_x==0) && (n_y==0) && (n_z==0)) // DC mode
                    LOWRES_density_perturb[0] = 0;
                else{
                    if(VELOCITY_COMPONENT==1) // x-component
                        LOWRES_density_perturb[HII_C_INDEX(n_x,n_y,n_z)] *= dDdt_over_D*k_x*I/k_sq/(HII_TOT_NUM_PIXELS+0.0);
                    else if (VELOCITY_COMPONENT == 2)
                        LOWRES_density_perturb[HII_C_INDEX(n_x,n_y,n_z)] *= dDdt_over_D*k_y*I/k_sq/(HII_TOT_NUM_PIXELS+0.0);
                    else
                        LOWRES_density_perturb[HII_C_INDEX(n_x,n_y,n_z)] *= dDdt_over_D*k_z*I/k_sq/(HII_TOT_NUM_PIXELS+0.0);
                }
            }
        }
    }
    plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)LOWRES_density_perturb, (float *)LOWRES_density_perturb, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();

    for (i=0; i<HII_DIM; i++){
        for (j=0; j<HII_DIM; j++){
            for (k=0; k<HII_DIM; k++){
                *((float *)LOWRES_velocity_REDSHIFT + HII_R_INDEX(i,j,k)) = *((float *)LOWRES_density_perturb + HII_R_FFT_INDEX(i,j,k));
            }
        }
    }

    // deallocate
    SAFE_FFTWF_FREE(LOWRES_density_perturb);
    SAFE_FFTWF_FREE(LOWRES_density_perturb_saved);
}


void ReadFcollTable() {

    char filename[500];
    char dummy_string[500];
    FILE *F;

    double *PARAM_VALS_FCOLL = calloc(TOTAL_AVAILABLE_PARAMS_FCOLL_TABLE,sizeof(double));

    int i;

    // Note the below code does not yet work for ALPHA != 0

    sprintf(filename,"f_coll_lightcone_data_%d_%.0fMpc.txt", HII_DIM, BOX_LEN);
    F = fopen(filename,"rt");

    for(i=0;i<TOTAL_AVAILABLE_PARAMS_FCOLL_TABLE;i++) {
        fscanf(F,"%s\t%lf\n",&dummy_string,&PARAM_VALS_FCOLL[i]);
    }
    fclose(F);

    R_MFP_UB = PARAM_VALS_FCOLL[0];
    TVIR_LB_FCOLL = PARAM_VALS_FCOLL[1];
    TVIR_UB_FCOLL = PARAM_VALS_FCOLL[2];
    ZETA_PL_LB = PARAM_VALS_FCOLL[3];
    ZETA_PL_UB = PARAM_VALS_FCOLL[4];
    R_MFP_STEPS = (int)PARAM_VALS_FCOLL[5];
    TVIR_STEPS = (int)PARAM_VALS_FCOLL[6];
    PL_STEPS = (int)PARAM_VALS_FCOLL[7];

    R_MFP_MIN = fmax(R_BUBBLE_MIN, (L_FACTOR*BOX_LEN/(float)HII_DIM));

    R_MFP_BINWIDTH = ( R_MFP_UB - R_MFP_MIN )/((float)R_MFP_STEPS - 1.);
    TVIR_BINWIDTH = ( TVIR_UB_FCOLL - TVIR_LB_FCOLL )/((float)TVIR_STEPS - 1.);

    PL_BINWIDTH = ( ZETA_PL_UB - ZETA_PL_LB )/((float)PL_STEPS - 1.);

    TVIR_INT_1 = (int)floor((log10(ION_Tvir_MIN) - TVIR_LB_FCOLL)/TVIR_BINWIDTH);
    TVIR_INT_2 = TVIR_INT_1 + 1;

    TVIR_VAL_1 = TVIR_LB_FCOLL + ( TVIR_UB_FCOLL - TVIR_LB_FCOLL )*(float)TVIR_INT_1/((float)TVIR_STEPS - 1.);
    TVIR_VAL_2 = TVIR_LB_FCOLL + ( TVIR_UB_FCOLL - TVIR_LB_FCOLL )*(float)TVIR_INT_2/((float)TVIR_STEPS - 1.);

    ZETA_PL_INT_1 = (int)floor((EFF_FACTOR_PL_INDEX - ZETA_PL_LB)/PL_BINWIDTH);
    ZETA_PL_INT_2 = ZETA_PL_INT_1 + 1;

    ZETA_PL_VAL_1 = ZETA_PL_LB + ( ZETA_PL_UB - ZETA_PL_LB )*(float)ZETA_PL_INT_1/((float)PL_STEPS - 1.);
    ZETA_PL_VAL_2 = ZETA_PL_LB + ( ZETA_PL_UB - ZETA_PL_LB )*(float)ZETA_PL_INT_2/((float)PL_STEPS - 1.);

    if(INCLUDE_ZETA_PL) {
        SIZE_FIRST = N_USER_REDSHIFT*R_MFP_STEPS*TVIR_STEPS*PL_STEPS;
        SIZE_FINAL = N_USER_REDSHIFT*TVIR_STEPS*PL_STEPS;
    }
    else {
        SIZE_FIRST = N_USER_REDSHIFT*R_MFP_STEPS*TVIR_STEPS;
        SIZE_FINAL = N_USER_REDSHIFT*TVIR_STEPS;
    }

    Ionisation_fcoll_table = calloc(SIZE_FIRST,sizeof(double));
    Ionisation_fcoll_table_final = calloc(SIZE_FINAL,sizeof(double));

    // NOTE: No support for the power law index for zeta at the present time

    if(INCLUDE_ZETA_PL) {
        sprintf(filename, "Ionisation_fcoll_table_final_Rmax%f_Tmin%f_Tmax%f_PLmin%1.6f_PLmax%1.6f_%d_%.0fMpc", R_MFP_UB, TVIR_LB_FCOLL, TVIR_UB_FCOLL, ZETA_PL_LB, ZETA_PL_UB, HII_DIM, BOX_LEN);
    }
    else {
        sprintf(filename, "Ionisation_fcoll_table_final_Rmax%f_Tmin%f_Tmax%f_PL%1.6f_%d_%.0fMpc", R_MFP_UB, TVIR_LB_FCOLL, TVIR_UB_FCOLL, EFF_FACTOR_PL_INDEX, HII_DIM, BOX_LEN);
    }
    F = fopen(filename, "rb");
    fread(Ionisation_fcoll_table_final,N_USER_REDSHIFT*TVIR_STEPS,sizeof(double),F);
    fclose(F);

    if(INCLUDE_ZETA_PL) {
        sprintf(filename, "Ionisation_fcoll_table_Rmax%f_Tmin%f_Tmax%f_PLmin%1.6f_PLmax%1.6f_%d_%.0fMpc", R_MFP_UB, TVIR_LB_FCOLL, TVIR_UB_FCOLL, ZETA_PL_LB, ZETA_PL_UB, HII_DIM, BOX_LEN);
    }
    else {
        sprintf(filename, "Ionisation_fcoll_table_Rmax%f_Tmin%f_Tmax%f_PL%1.6f_%d_%.0fMpc", R_MFP_UB, TVIR_LB_FCOLL, TVIR_UB_FCOLL, EFF_FACTOR_PL_INDEX, HII_DIM, BOX_LEN);
    }
    F = fopen(filename, "rb");
    fread(Ionisation_fcoll_table,N_USER_REDSHIFT*R_MFP_STEPS*TVIR_STEPS,sizeof(double),F);
    fclose(F);
}



void GeneratePS(int CO_EVAL, double AverageTb) {

    fftwf_plan plan;

    int i,j,k,n_x, n_y, n_z,skip_zero_mode;
    float k_x, k_y, k_z, k_mag;
    double ave;
    unsigned long long ct;

    if(!CO_EVAL) {

        ave = 0.0;
        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    ave += box_interpolate[HII_R_INDEX(i,j,k)];
                }
            }
        }
        ave /= (HII_TOT_NUM_PIXELS+0.0);
    }
    else {
        ave = AverageTb;
    }

    for (ct=0; ct<NUM_BINS; ct++){
        p_box[ct] = k_ave[ct] = 0;
        in_bin_ct[ct] = 0;
    }

    // fill-up the real-space of the deldel box
    for (i=0; i<HII_DIM; i++){
        for (j=0; j<HII_DIM; j++){
            for (k=0; k<HII_DIM; k++){
                if(!CO_EVAL) {
                    *((float *)deldel_T_LC + HII_R_FFT_INDEX(i,j,k)) = (box_interpolate[HII_R_INDEX(i,j,k)]/ave - 1)*VOLUME/(HII_TOT_NUM_PIXELS+0.0);
                }
                else {
                    *((float *)deldel_T_LC + HII_R_FFT_INDEX(i,j,k)) = (delta_T[HII_R_INDEX(i,j,k)]/ave - 1)*VOLUME/(HII_TOT_NUM_PIXELS+0.0);
                }
                if (DIMENSIONAL_T_POWER_SPEC){
                    *((float *)deldel_T_LC + HII_R_FFT_INDEX(i,j,k)) *= ave;
                }
                // Note: we include the V/N factor for the scaling after the fft
            }
        }
    }

    // transform to k-space
    plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deldel_T_LC, (fftwf_complex *)deldel_T_LC, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // If the light-cone 21cm PS is to be calculated, one should avoid the k(k_x = 0, k_y = 0, k_z) modes (see Datta et al. 2012).
    if(!CO_EVAL) {

        // now construct the power spectrum file
        for (n_x=0; n_x<HII_DIM; n_x++){
            if (n_x>HII_MIDDLE)
                k_x =(n_x-HII_DIM) * DELTA_K;  // wrap around for FFT convention
            else
                k_x = n_x * DELTA_K;

            for (n_y=0; n_y<HII_DIM; n_y++){

                // avoid the k(k_x = 0, k_y = 0, k_z) modes
                if(n_x != 0 && n_y != 0) {

                    if (n_y>HII_MIDDLE)
                        k_y =(n_y-HII_DIM) * DELTA_K;
                    else
                        k_y = n_y * DELTA_K;

                    for (n_z=0; n_z<=HII_MIDDLE; n_z++){
                        k_z = n_z * DELTA_K;

                        k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                        // now go through the k bins and update
                        ct = 0;
                        k_floor = 0;
                        k_ceil = k_first_bin_ceil;
                        while (k_ceil < k_max){
                            // check if we fal in this bin
                            if ((k_mag>=k_floor) && (k_mag < k_ceil)){
                                in_bin_ct[ct]++;
                                p_box[ct] += pow(k_mag,3)*pow(cabs(deldel_T_LC[HII_C_INDEX(n_x, n_y, n_z)]), 2)/(2.0*PI*PI*VOLUME);
                                // note the 1/VOLUME factor, which turns this into a power density in k-space

                                k_ave[ct] += k_mag;
                                break;
                            }

                            ct++;
                            k_floor=k_ceil;
                            k_ceil*=k_factor;
                        }
                    }
                }
            }
        } // end looping through k box

    }
    else {

        // Co-eval box, so should sample the entire cube

        // now construct the power spectrum file
        for (n_x=0; n_x<HII_DIM; n_x++){
            if (n_x>HII_MIDDLE)
                k_x =(n_x-HII_DIM) * DELTA_K;  // wrap around for FFT convention
            else
                k_x = n_x * DELTA_K;

            for (n_y=0; n_y<HII_DIM; n_y++){

                if (n_y>HII_MIDDLE)
                    k_y =(n_y-HII_DIM) * DELTA_K;
                else
                    k_y = n_y * DELTA_K;

                for (n_z=0; n_z<=HII_MIDDLE; n_z++){
                    k_z = n_z * DELTA_K;

                    k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                    // now go through the k bins and update
                    ct = 0;
                    k_floor = 0;
                    k_ceil = k_first_bin_ceil;
                    while (k_ceil < k_max){
                        // check if we fal in this bin
                        if ((k_mag>=k_floor) && (k_mag < k_ceil)){
                            in_bin_ct[ct]++;
                            p_box[ct] += pow(k_mag,3)*pow(cabs(deldel_T_LC[HII_C_INDEX(n_x, n_y, n_z)]), 2)/(2.0*PI*PI*VOLUME);
                            // note the 1/VOLUME factor, which turns this into a power density in k-space

                            k_ave[ct] += k_mag;
                            break;
                        }

                        ct++;
                        k_floor=k_ceil;
                        k_ceil*=k_factor;
                    }
                }
            }
        } // end looping through k box
    }
}



//JBM:anisotropic version of the power spectrum
void GeneratePS_aniso(int CO_EVAL, double AverageTb, char* filename) {


	float mu_k; //JBM: mu (LoS)
  float mu_floor, mu_ceil; //the floor and ceiling of each bin
  int Nmu_bins = 10; //number of mu bins, from 0 to 1.
  int i;


  float **mu_k_avg; //the actual mu average inside each bin
  mu_k_avg = (float **)calloc(NUM_BINS,sizeof(float *));
  for(i=0;i<NUM_BINS;i++) {
      mu_k_avg[i] = (float *)calloc(Nmu_bins,sizeof(float));
  }

  float **k_avg; //the actual mu average inside each bin
  k_avg = (float **)calloc(NUM_BINS,sizeof(float *));
  for(i=0;i<NUM_BINS;i++) {
      k_avg[i] = (float *)calloc(Nmu_bins,sizeof(float));
  }


  float **pow21; //power spectrum and
  pow21 = (float **)calloc(NUM_BINS,sizeof(float *));
  for(i=0;i<NUM_BINS;i++) {
      pow21[i] = (float *)calloc(Nmu_bins,sizeof(float));
  }



  float **counter_kmu; //how many ks in each bin
  counter_kmu = (float **)calloc(NUM_BINS,sizeof(float *));
  for(i=0;i<NUM_BINS;i++) {
      counter_kmu[i] = (float *)calloc(Nmu_bins,sizeof(float));
  }














    fftwf_plan plan;

    int j,k,n_x, n_y, n_z,skip_zero_mode;
    float k_x, k_y, k_z, k_mag;
    double ave;
    unsigned long long ct;
    unsigned long long ct_mu;

    if(!CO_EVAL) {

        ave = 0.0;
        for (i=0; i<HII_DIM; i++){
            for (j=0; j<HII_DIM; j++){
                for (k=0; k<HII_DIM; k++){
                    ave += box_interpolate[HII_R_INDEX(i,j,k)];
                }
            }
        }
        ave /= (HII_TOT_NUM_PIXELS+0.0);
    }
    else {
        ave = AverageTb;
    }

    for (ct=0; ct<NUM_BINS; ct++){
        p_box[ct] = k_ave[ct] = 0;
        in_bin_ct[ct] = 0;
    }

    // fill-up the real-space of the deldel box
    for (i=0; i<HII_DIM; i++){
        for (j=0; j<HII_DIM; j++){
            for (k=0; k<HII_DIM; k++){
                if(!CO_EVAL) {
                    *((float *)deldel_T_LC + HII_R_FFT_INDEX(i,j,k)) = (box_interpolate[HII_R_INDEX(i,j,k)]/ave - 1)*VOLUME/(HII_TOT_NUM_PIXELS+0.0);
                }
                else {
                    *((float *)deldel_T_LC + HII_R_FFT_INDEX(i,j,k)) = (delta_T[HII_R_INDEX(i,j,k)]/ave - 1)*VOLUME/(HII_TOT_NUM_PIXELS+0.0);
                }
                if (DIMENSIONAL_T_POWER_SPEC){
                    *((float *)deldel_T_LC + HII_R_FFT_INDEX(i,j,k)) *= ave;
                }
                // Note: we include the V/N factor for the scaling after the fft
            }
        }
    }

    // transform to k-space
    plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deldel_T_LC, (fftwf_complex *)deldel_T_LC, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // If the light-cone 21cm PS is to be calculated, one should avoid the k(k_x = 0, k_y = 0, k_z) modes (see Datta et al. 2012).
    if(!CO_EVAL) {

        // now construct the power spectrum file
        for (n_x=0; n_x<HII_DIM; n_x++){
            if (n_x>HII_MIDDLE)
                k_x =(n_x-HII_DIM) * DELTA_K;  // wrap around for FFT convention
            else
                k_x = n_x * DELTA_K;

            for (n_y=0; n_y<HII_DIM; n_y++){

                // avoid the k(k_x = 0, k_y = 0, k_z) modes
                if(n_x != 0 && n_y != 0) {

                    if (n_y>HII_MIDDLE)
                        k_y =(n_y-HII_DIM) * DELTA_K;
                    else
                        k_y = n_y * DELTA_K;

                    for (n_z=0; n_z<=HII_MIDDLE; n_z++){
                        k_z = n_z * DELTA_K;

                        k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                        // now go through the k bins and update
                        ct = 0;
                        k_floor = 0;
                        k_ceil = k_first_bin_ceil;




                        while (k_ceil < k_max){
                            // check if we fal in this bin
                            if ((k_mag>=k_floor) && (k_mag < k_ceil)){
                                in_bin_ct[ct]++;
                                p_box[ct] += pow(k_mag,3)*pow(cabs(deldel_T_LC[HII_C_INDEX(n_x, n_y, n_z)]), 2)/(2.0*PI*PI*VOLUME);
                                // note the 1/VOLUME factor, which turns this into a power density in k-space

                                k_ave[ct] += k_mag;
                                break;
                            }

                            ct++;
                            k_floor=k_ceil;
                            k_ceil*=k_factor;
                        }
                    }
                }
            }
        } // end looping through k box

    }
    else {

        // Co-eval box, so should sample the entire cube

        // now construct the power spectrum file
        for (n_x=0; n_x<HII_DIM; n_x++){
            if (n_x>HII_MIDDLE)
                k_x =(n_x-HII_DIM) * DELTA_K;  // wrap around for FFT convention
            else
                k_x = n_x * DELTA_K;

            for (n_y=0; n_y<HII_DIM; n_y++){

                if (n_y>HII_MIDDLE)
                    k_y =(n_y-HII_DIM) * DELTA_K;
                else
                    k_y = n_y * DELTA_K;

                for (n_z=0; n_z<=HII_MIDDLE; n_z++){
                    k_z = n_z * DELTA_K;

                    k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);
                    if(k_mag==0){
                    	mu_k = 0.0;
                    }
                    else{
                    	mu_k = fabs(k_z/k_mag);
                    }



                    // now go through the k and mu bins and update
                    ct = 0;
                    k_floor = 0;
                    k_ceil = k_first_bin_ceil;

                    while (k_ceil < k_max){
                      //JBM also for mu
                      ct_mu = 0;
                      mu_floor = 0.0;
                      mu_ceil = 1.0/(Nmu_bins-1.0);

                      while (mu_ceil <= 1.0){


                        // printf("k=%.1e mu=%.1e \n", k_mag, mu_k);
                        // printf("kfl=%.1e mufl=%.1e \n", k_floor, mu_floor);
                        // printf("kcel=%.1e muceil=%.1e \n", k_ceil, mu_ceil);


                        // check if we fall in each bin
                        if ((k_mag>=k_floor) && (k_mag < k_ceil)
                              && (mu_k>=mu_floor) && (mu_k < mu_ceil)){

                            counter_kmu[ct][ct_mu]++;

                            pow21[ct][ct_mu] += pow(k_mag,3)*pow(cabs(deldel_T_LC[HII_C_INDEX(n_x, n_y, n_z)]), 2)/(2.0*PI*PI*VOLUME);
                            // note the 1/VOLUME factor, which turns this into a power density in k-space

                            k_avg[ct][ct_mu] += k_mag;
                            mu_k_avg[ct][ct_mu] += mu_k;



                            break;



                        }

                        ct_mu++;
                        mu_floor=mu_ceil;
                        mu_ceil+=1.0/(Nmu_bins-1.0);

                      }

                      ct++;
                      k_floor=k_ceil;
                      k_ceil*=k_factor;


                    }



                }
            }
        } // end looping through k box
    }


  	int lengthname=3000;
  	FILE *fp;
  	fp = fopen(filename,"w");

    for(ct=0;ct<NUM_BINS;ct++) {
      for(ct_mu=0;ct_mu<Nmu_bins;ct_mu++){
        if(counter_kmu[ct][ct_mu] > 0){
          fprintf(fp, "%.3e %.3e %.3e %.3e \n", k_avg[ct][ct_mu]/(counter_kmu[ct][ct_mu]+0.0) ,
          mu_k_avg[ct][ct_mu]/(counter_kmu[ct][ct_mu]+0.0), pow21[ct][ct_mu]/(counter_kmu[ct][ct_mu]+0.0)
              , pow21[ct][ct_mu]/(counter_kmu[ct][ct_mu]+0.0)/sqrt(counter_kmu[ct][ct_mu]+0.0));
        }
      }
    }

    fclose(fp);


    for(i=0;i<NUM_BINS;i++) {
        SAFE_FREE(pow21[i]);
    }
    SAFE_FREE(pow21);
    for(i=0;i<NUM_BINS;i++) {
        SAFE_FREE(counter_kmu[i]);
    }
    SAFE_FREE(counter_kmu);
    for(i=0;i<NUM_BINS;i++) {
        SAFE_FREE(k_avg[i]);
    }
    SAFE_FREE(k_avg);
    for(i=0;i<NUM_BINS;i++) {
        SAFE_FREE(mu_k_avg[i]);
    }
    SAFE_FREE(mu_k_avg);




}

//JBM: output velocity power spectrum, just a diagnostic
//note that it's P(v) and not P(v^2).


void GeneratePS_vel() {

    fftwf_plan plan;

    int i,j,k,n_x, n_y, n_z,skip_zero_mode;
    float k_x, k_y, k_z, k_mag;
    double ave;
    unsigned long long ct;

    double ave_vel = 26; //km/s, average velocity

    for (ct=0; ct<NUM_BINS; ct++){
        p_box[ct] = k_ave[ct] = 0;
        in_bin_ct[ct] = 0;
    }

    // fill-up the real-space of the deldel box
    for (i=0; i<HII_DIM; i++){
        for (j=0; j<HII_DIM; j++){
            for (k=0; k<HII_DIM; k++){
                  *((float *)deldel_T_LC + HII_R_FFT_INDEX(i,j,k)) = (smoothed_box_vcb[HII_R_INDEX(i,j,k)]/ave_vel - 1)*VOLUME/(HII_TOT_NUM_PIXELS+0.0);

                if (DIMENSIONAL_T_POWER_SPEC){
                    *((float *)deldel_T_LC + HII_R_FFT_INDEX(i,j,k)) *= ave_vel;
                }
                // Note: we include the V/N factor for the scaling after the fft
            }
        }
    }

    // transform to k-space
    plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deldel_T_LC, (fftwf_complex *)deldel_T_LC, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // If the light-cone 21cm PS is to be calculated, one should avoid the k(k_x = 0, k_y = 0, k_z) modes (see Datta et al. 2012).
        // Co-eval box, so should sample the entire cube

        // now construct the power spectrum file
        for (n_x=0; n_x<HII_DIM; n_x++){
            if (n_x>HII_MIDDLE)
                k_x =(n_x-HII_DIM) * DELTA_K;  // wrap around for FFT convention
            else
                k_x = n_x * DELTA_K;

            for (n_y=0; n_y<HII_DIM; n_y++){

                if (n_y>HII_MIDDLE)
                    k_y =(n_y-HII_DIM) * DELTA_K;
                else
                    k_y = n_y * DELTA_K;

                for (n_z=0; n_z<=HII_MIDDLE; n_z++){
                    k_z = n_z * DELTA_K;

                    k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

                    // now go through the k bins and update
                    ct = 0;
                    k_floor = 0;
                    k_ceil = k_first_bin_ceil;
                    while (k_ceil < k_max){
                        // check if we fal in this bin
                        if ((k_mag>=k_floor) && (k_mag < k_ceil)){
                            in_bin_ct[ct]++;
                            p_box[ct] += pow(k_mag,3)*pow(cabs(deldel_T_LC[HII_C_INDEX(n_x, n_y, n_z)]), 2)/(2.0*PI*PI*VOLUME);
                            // note the 1/VOLUME factor, which turns this into a power density in k-space

                            k_ave[ct] += k_mag;
                            break;
                        }

                        ct++;
                        k_floor=k_ceil;
                        k_ceil*=k_factor;
                    }
                }
            }
        } // end looping through k box

}

/**** Arrays declared and used *****/

void init_21cmMC_HII_arrays() {

    Overdense_spline_GL_low = calloc(Nlow,sizeof(float));
    Fcoll_spline_GL_low = calloc(Nlow,sizeof(float));
    second_derivs_low_GL = calloc(Nlow,sizeof(float));
    Overdense_spline_GL_high = calloc(Nhigh,sizeof(float));
    Fcoll_spline_GL_high = calloc(Nhigh,sizeof(float));
    second_derivs_high_GL = calloc(Nhigh,sizeof(float));

    deltax_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    deltax_unfiltered_original = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    deltax_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    xe_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    xe_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    deldel_T = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    deldel_T_LC = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

    deltax = (float *) calloc(HII_TOT_FFT_NUM_PIXELS,sizeof(float));
    Fcoll = (float *) calloc(HII_TOT_NUM_PIXELS,sizeof(float));
    xH = (float *)calloc(HII_TOT_NUM_PIXELS,sizeof(float));
    v = (float *) calloc(HII_TOT_FFT_NUM_PIXELS,sizeof(float));
    vel_gradient = (float *) calloc(HII_TOT_FFT_NUM_PIXELS,sizeof(float));
    delta_T = (float *) calloc(HII_TOT_NUM_PIXELS,sizeof(float));

    x_pos = calloc(N_RSD_STEPS,sizeof(float));
    x_pos_offset = calloc(N_RSD_STEPS,sizeof(float));
    delta_T_RSD_LOS = calloc(HII_DIM,sizeof(float));

    xi_low = calloc((NGLlow+1),sizeof(float));
    wi_low = calloc((NGLlow+1),sizeof(float));

    xi_high = calloc((NGLhigh+1),sizeof(float));
    wi_high = calloc((NGLhigh+1),sizeof(float));


    //JBM we want finer binning return to k_factor = 1.35; if want older results.
    // k_factor = 1.1;
    k_factor = 1.05;
    k_first_bin_ceil = DELTA_K;
    k_max = DELTA_K*HII_DIM;
    // initialize arrays
    // ghetto counting (lookup how to do logs of arbitrary bases in c...)
    NUM_BINS = 0;
    k_floor = 0;
    k_ceil = k_first_bin_ceil;
    while (k_ceil < k_max){
        NUM_BINS++;
        k_floor=k_ceil;
        k_ceil*=k_factor;
    }

    p_box = calloc(NUM_BINS,sizeof(double));
    k_ave = calloc(NUM_BINS,sizeof(double));
    in_bin_ct = (unsigned long long *)calloc(NUM_BINS,sizeof(unsigned long long));
}

// void init_21cmMC_Ts_arrays() {

//     int i,j;

//     fprintf(stderr, "Allocating box with HII_KSPACE_NUM_PIXELS = %d\n", HII_KSPACE_NUM_PIXELS);
//     fprintf(stderr, "Allocating box with sizeof(fftwf_complex) = %d\n", sizeof(fftwf_complex));

//     box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
//     if (box == NULL){
//         fprintf(stderr, "Box is not properly allocated");
//     }
//     unfiltered_box = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

//     //JBM:
//     box_vcb = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
//     unfiltered_vcb_box = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);

//     Tk_box = (float *) calloc(HII_TOT_NUM_PIXELS,sizeof(float));
//     x_e_box = (float *) calloc(HII_TOT_NUM_PIXELS,sizeof(float));
//     Ts = (float *) calloc(HII_TOT_NUM_PIXELS,sizeof(float));

//     inverse_diff = calloc(x_int_NXHII,sizeof(float));

//     zpp_growth = (float *)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));

//     fcoll_R_grid = (double ***)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double **));
//     dfcoll_dz_grid = (double ***)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double **));
//     for(i=0;i<NUM_FILTER_STEPS_FOR_Ts;i++) {
//         fcoll_R_grid[i] = (double **)calloc(zpp_interp_points,sizeof(double *));
//         dfcoll_dz_grid[i] = (double **)calloc(zpp_interp_points,sizeof(double *));
//         for(j=0;j<zpp_interp_points;j++) {
//             fcoll_R_grid[i][j] = (double *)calloc(dens_Ninterp,sizeof(double));
//             dfcoll_dz_grid[i][j] = (double *)calloc(dens_Ninterp,sizeof(double));
//         }
//     }

//     fcoll_R_array = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//     Sigma_Tmin_grid = (double *)calloc(zpp_interp_points,sizeof(double));

//     grid_dens = (double **)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double *));
//     for(i=0;i<NUM_FILTER_STEPS_FOR_Ts;i++) {
//         grid_dens[i] = (double *)calloc(dens_Ninterp,sizeof(double));
//     }

//     density_gridpoints = (double **)calloc(dens_Ninterp,sizeof(double *));
//     for(i=0;i<dens_Ninterp;i++) {
//         density_gridpoints[i] = (double *)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//     }
//     ST_over_PS_arg_grid = (double *)calloc(zpp_interp_points,sizeof(double));

//     //JBM: we define the F_ST(v) 2D array.
//     logFcoll_vcb = (double **)calloc(NZINT,sizeof(double *));
//     for(i=0;i<NZINT;i++) {
//         logFcoll_vcb[i] = (double *)calloc(NVINT,sizeof(double));
//     }

//     //JBM: and sigma(z,v,M_cool(z,v)).
//     sigmacool_vcb = (double **)calloc(NZINT,sizeof(double *));
//     for(i=0;i<NZINT;i++) {
//         sigmacool_vcb[i] = (double *)calloc(NVINT,sizeof(double));
//     }




//     dens_grid_int_vals = (short **)calloc(HII_TOT_NUM_PIXELS,sizeof(short *));
//     delNL0_rev = (float **)calloc(HII_TOT_NUM_PIXELS,sizeof(float *));
//     for(i=0;i<HII_TOT_NUM_PIXELS;i++) {
//         dens_grid_int_vals[i] = (short *)calloc((float)NUM_FILTER_STEPS_FOR_Ts,sizeof(short));
//         delNL0_rev[i] = (float *)calloc((float)NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
//     }

//     //JBM:values of velocity smoothed over some scale R
//     vcb_rev = (float **)calloc(HII_TOT_NUM_PIXELS,sizeof(float *));
//     for(i=0;i<HII_TOT_NUM_PIXELS;i++) {
//         vcb_rev[i] = (float *)calloc((float)NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
//     }

//     fcoll_interp1 = (double **)calloc(dens_Ninterp,sizeof(double *));
//     fcoll_interp2 = (double **)calloc(dens_Ninterp,sizeof(double *));
//     dfcoll_interp1 = (double **)calloc(dens_Ninterp,sizeof(double *));
//     dfcoll_interp2 = (double **)calloc(dens_Ninterp,sizeof(double *));
//     for(i=0;i<dens_Ninterp;i++) {
//         fcoll_interp1[i] = (double *)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//         fcoll_interp2[i] = (double *)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//         dfcoll_interp1[i] = (double *)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//         dfcoll_interp2[i] = (double *)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//     }

//     zpp_edge = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//     sigma_atR = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//     sigma_Tmin = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//     ST_over_PS = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//     sum_lyn = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));

//     zpp_for_evolve_list = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
//     R_values = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
//     SingleVal_float = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));

//     delNL0_bw = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
//     delNL0_Offset = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
//     delNL0_LL = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
//     delNL0_UL = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
//     delNL0_ibw = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
//     log10delNL0_diff = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));
//     log10delNL0_diff_UL = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(float));

//     freq_int_heat_tbl = (double **)calloc(x_int_NXHII,sizeof(double *));
//     freq_int_ion_tbl = (double **)calloc(x_int_NXHII,sizeof(double *));
//     freq_int_lya_tbl = (double **)calloc(x_int_NXHII,sizeof(double *));
//     freq_int_heat_tbl_diff = (double **)calloc(x_int_NXHII,sizeof(double *));
//     freq_int_ion_tbl_diff = (double **)calloc(x_int_NXHII,sizeof(double *));
//     freq_int_lya_tbl_diff = (double **)calloc(x_int_NXHII,sizeof(double *));
//     for(i=0;i<x_int_NXHII;i++) {
//         freq_int_heat_tbl[i] = (double *)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//         freq_int_ion_tbl[i] = (double *)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//         freq_int_lya_tbl[i] = (double *)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//         freq_int_heat_tbl_diff[i] = (double *)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//         freq_int_ion_tbl_diff[i] = (double *)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//         freq_int_lya_tbl_diff[i] = (double *)calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));
//     }

//     dstarlya_dt_prefactor = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(double));

//     SingleVal_int = calloc(NUM_FILTER_STEPS_FOR_Ts,sizeof(short));
// }

void destroy_21cmMC_HII_arrays(int skip_deallocate) {
    SAFE_FFTWF_FREE(deltax_unfiltered);
    SAFE_FFTWF_FREE(deltax_unfiltered_original);
    SAFE_FFTWF_FREE(deltax_filtered);
    SAFE_FFTWF_FREE(deldel_T);
    SAFE_FFTWF_FREE(deldel_T_LC);
    SAFE_FFTWF_FREE(xe_unfiltered);
    SAFE_FFTWF_FREE(xe_filtered);

    SAFE_FREE(xH);
    SAFE_FREE(deltax);
    SAFE_FREE(Fcoll);
    SAFE_FREE(delta_T);
    SAFE_FREE(v);
    SAFE_FREE(vel_gradient);
    SAFE_FREE(p_box);
    SAFE_FREE(k_ave);
    SAFE_FREE(in_bin_ct);

    SAFE_FREE(x_pos);
    SAFE_FREE(x_pos_offset);
    SAFE_FREE(delta_T_RSD_LOS);

    SAFE_FREE(Overdense_spline_GL_low);
    SAFE_FREE(Fcoll_spline_GL_low);
    SAFE_FREE(second_derivs_low_GL);
    SAFE_FREE(Overdense_spline_GL_high);
    SAFE_FREE(Fcoll_spline_GL_high);
    SAFE_FREE(second_derivs_high_GL);

    SAFE_FREE(xi_low);
    SAFE_FREE(wi_low);

    SAFE_FREE(xi_high);
    SAFE_FREE(wi_high);

    if(skip_deallocate!=1) {
        SAFE_FREE(Mass_Spline);
        SAFE_FREE(Sigma_Spline);
        SAFE_FREE(dSigmadm_Spline);
        SAFE_FREE(second_derivs_sigma);
        SAFE_FREE(second_derivs_dsigma);
    }
}

void destroy_21cmMC_Ts_arrays() {
    int i, j;

    printf("Freeing box and related pointers\n");

    SAFE_FFTWF_FREE(box);
    SAFE_FFTWF_FREE(unfiltered_box);
    SAFE_FFTWF_FREE(box_vcb);
    SAFE_FFTWF_FREE(unfiltered_vcb_box);

    SAFE_FREE(Tk_box);
    SAFE_FREE(x_e_box);
    SAFE_FREE(Ts);

    printf("Freeing fcoll_R_grid and dfcoll_dz_grid\n");

    if(fcoll_R_grid) {
        for(i = 0; i < NUM_FILTER_STEPS_FOR_Ts; i++) {
            if(fcoll_R_grid[i]) {
                for(j = 0; j < zpp_interp_points; j++) {
                    SAFE_FREE(fcoll_R_grid[i][j]);
                    SAFE_FREE(dfcoll_dz_grid[i][j]);
                }
                SAFE_FREE(fcoll_R_grid[i]);
                SAFE_FREE(dfcoll_dz_grid[i]);
            }
        }
        SAFE_FREE(fcoll_R_grid);
        SAFE_FREE(dfcoll_dz_grid);
    }

    printf("Freeing grid_dens\n");

    if(grid_dens) {
        for(i = 0; i < NUM_FILTER_STEPS_FOR_Ts; i++) {
            SAFE_FREE(grid_dens[i]);
        }
        SAFE_FREE(grid_dens);
    }

    printf("Freeing density_gridpoints\n");

    if(density_gridpoints) {
        for(i = 0; i < dens_Ninterp; i++) {
            SAFE_FREE(density_gridpoints[i]);
        }
        SAFE_FREE(density_gridpoints);
    }

    printf("Freeing logFcoll_vcb\n");

    if(logFcoll_vcb) {
        for(i = 0; i < NZINT; i++) {
            SAFE_FREE(logFcoll_vcb[i]);
        }
        SAFE_FREE(logFcoll_vcb);
    }

    printf("Freeing sigmacool_vcb\n");

    if(sigmacool_vcb) {
        for(i = 0; i < NZINT; i++) {
            SAFE_FREE(sigmacool_vcb[i]);
        }
        SAFE_FREE(sigmacool_vcb);
    }

    SAFE_FREE(ST_over_PS_arg_grid);
    SAFE_FREE(Sigma_Tmin_grid);
    SAFE_FREE(fcoll_R_array);
    SAFE_FREE(zpp_growth);
    SAFE_FREE(inverse_diff);

    printf("Freeing fcoll_interp and dfcoll_interp\n");

    if(fcoll_interp1) {
        for(i = 0; i < dens_Ninterp; i++) {
            SAFE_FREE(fcoll_interp1[i]);
            SAFE_FREE(fcoll_interp2[i]);
            SAFE_FREE(dfcoll_interp1[i]);
            SAFE_FREE(dfcoll_interp2[i]);
        }
        SAFE_FREE(fcoll_interp1);
        SAFE_FREE(fcoll_interp2);
        SAFE_FREE(dfcoll_interp1);
        SAFE_FREE(dfcoll_interp2);
    }

    printf("Freeing dens_grid_int_vals and delNL0_rev\n");

    if(dens_grid_int_vals) {
        for(i = 0; i < HII_TOT_NUM_PIXELS; i++) {
            SAFE_FREE(dens_grid_int_vals[i]);
            SAFE_FREE(delNL0_rev[i]);
        }
        SAFE_FREE(dens_grid_int_vals);
        SAFE_FREE(delNL0_rev);
    }

    printf("Freeing vcb_rev\n");

    if(vcb_rev) {
        for(i = 0; i < HII_TOT_NUM_PIXELS; i++) {
            SAFE_FREE(vcb_rev[i]);
        }
        SAFE_FREE(vcb_rev);
    }

    SAFE_FREE(delNL0_bw);
    SAFE_FREE(zpp_for_evolve_list);
    SAFE_FREE(R_values);
    SAFE_FREE(delNL0_Offset);
    SAFE_FREE(delNL0_LL);
    SAFE_FREE(delNL0_UL);
    SAFE_FREE(SingleVal_int);
    SAFE_FREE(SingleVal_float);
    SAFE_FREE(dstarlya_dt_prefactor);
    SAFE_FREE(delNL0_ibw);
    SAFE_FREE(log10delNL0_diff);
    SAFE_FREE(log10delNL0_diff_UL);

    SAFE_FREE(zpp_edge);
    SAFE_FREE(sigma_atR);
    SAFE_FREE(sigma_Tmin);
    SAFE_FREE(ST_over_PS);
    SAFE_FREE(sum_lyn);

    printf("Freeing freq_int tables\n");

    if(freq_int_heat_tbl) {
        for(i = 0; i < x_int_NXHII; i++) {
            SAFE_FREE(freq_int_heat_tbl[i]);
            SAFE_FREE(freq_int_ion_tbl[i]);
            SAFE_FREE(freq_int_lya_tbl[i]);
            SAFE_FREE(freq_int_heat_tbl_diff[i]);
            SAFE_FREE(freq_int_ion_tbl_diff[i]);
            SAFE_FREE(freq_int_lya_tbl_diff[i]);
        }
        SAFE_FREE(freq_int_heat_tbl);
        SAFE_FREE(freq_int_ion_tbl);
        SAFE_FREE(freq_int_lya_tbl);
        SAFE_FREE(freq_int_heat_tbl_diff);
        SAFE_FREE(freq_int_ion_tbl_diff);
        SAFE_FREE(freq_int_lya_tbl_diff);
    }

    printf("Calling free_FcollTable\n");

    free_FcollTable();

    printf("Finished destroying arrays\n");
}

