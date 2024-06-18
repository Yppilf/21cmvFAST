#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/Variables.h"

// Helper macro to safely free memory and set pointer to NULL
#define SAFE_FREE(ptr) if(ptr) { free(ptr); ptr = NULL; }
#define SAFE_FFTWF_FREE(ptr) if(ptr) { fftwf_free(ptr); ptr = NULL; }

int main(int argc, char ** argv){

    char filename[500];
    char dummy_string[500];
    char tau_idx[10];
    FILE *F;
    
    int i,N_REDSHIFTS;
    
    float taue;
    
    INDIVIDUAL_ID = atof(argv[1]);
    INDIVIDUAL_ID_2 = atof(argv[2]);
    tau_idx = atof(argv[3]);
    
    double *PARAM_COSMOLOGY_VALS = calloc(TOTAL_COSMOLOGY_FILEPARAMS,sizeof(double));
    
    /////////////////   Read in the cosmological parameter data     /////////////////
    
    sprintf(filename,"%sWalkerCosmology_%1.6lf_%1.6lf.txt",WALKER_FOLDER,INDIVIDUAL_ID,INDIVIDUAL_ID_2);
    // fprintf(stderr, "\nReading file %s\n", filename);
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

    // fprintf(stderr, "Assigned values from cosmology\n");
    
    // Minus 3 for argv[0]=0, Random ID (argv[1]) and Zeta (argv[2])
    // Use Zeta just incase the Random ID ends up the same (shouldn't happen)
    N_REDSHIFTS = (argc-3)/2;
    
    float *Redshifts = calloc(N_REDSHIFTS,sizeof(float));
    float *xH = calloc(N_REDSHIFTS,sizeof(float));
    // fprintf(stderr, "Memory allocated\n");
    
    for(i=0;i<N_REDSHIFTS;i++) {
        Redshifts[i] = (float)(atof(argv[2*i+3]));
        xH[i] = (float)(atof(argv[2*(i+2)]));
    }

    // fprintf(stderr, "Read data into arrays\n");
    
    taue = tau_e(0, Redshifts[N_REDSHIFTS-1], Redshifts, xH, N_REDSHIFTS);

//    printf("Tau = %lf\n",taue);
    
    sprintf(filename, "%s/Tau_e_%s_%s_%s.txt",OUTPUT_FOLDER,argv[1],argv[2],tau_idx);
    fprintf(stderr, "Printing Tau data to %s\n",filename);
    F=fopen(filename, "wt");
    fprintf(F, "%lf\n",taue);
    fclose(F);
    
    SAFE_FREE(Redshifts);
    SAFE_FREE(xH);
    
    return 0;
}
