#include <fftw3.h>
#include <stdio.h>

fftwf_complex *box;
int HII_KSPACE_NUM_PIXELS = 13590000;    // Taken from drive_21cmMC_streamlined.c
int size = sizeof(fftwf_complex);

int main(){
    fprintf(stderr, "Size = %d\n", size);

    box = (fftwf_complex *) fftwf_malloc(size*HII_KSPACE_NUM_PIXELS);

    if (box == NULL){
        fprintf(stderr, "Memory allocation for box failed\n");
        return 1;
    }

    box[1359001][0] = 3.0;
    box[1359000][0] *= 3.0;

    fftwf_free(box);

    return 0;
}