#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "thresholds.h"

int main(int argc, char **argv)
{
    threshold = 1e-11;
    epsilon = 1e-6;

    struct matrix *mtrx_a;
    struct matrix *mtrx_b;

    FILE *mtrx_a_file = fopen("matrices/big_A.txt", "r");
    FILE *mtrx_b_file = fopen("matrices/big_B.txt", "r");
    
    if ((!mtrx_a_file) || (!mtrx_b_file)) {
	    printf("Couldn't open file.\n");
	    return -1;
    }
    
    if ((matrix_fread(mtrx_a_file, &mtrx_a)) ||
            (matrix_fread(mtrx_b_file, &mtrx_b))) {
	    printf("Couldn't read matrix.\n");
	    return -1;
    }
    
    fclose(mtrx_a_file);
    fclose(mtrx_b_file);

    qz_process(&mtrx_a, &mtrx_b);

    exit(0);
}
