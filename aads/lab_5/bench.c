
#include <stdio.h>
#include "matrix.h"

int main(int argc, char **argv)
{
    struct matrix *mtrx_a;
    struct matrix *mtrx_b;

    FILE *mtrx_a_file = fopen("matrix_a", "r");
    FILE *mtrx_b_file = fopen("matrix_b", "r");
    
    if ((!mtrx_a_file) || (!mtrx_b_file)) {
	printf("Couldn't open file.\n");
	return -1;
    }
    
    if ((matrix_fread(mtrx_a_file, &mtrx_a)) || (matrix_fread(mtrx_b_file, &mtrx_b))) {
	printf("Couldn't read matrix.\n");
	return -1;
    }

    fclose(mtrx_a_file);
    fclose(mtrx_b_file);
    
    return 0;
}
