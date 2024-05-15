#include <stdio.h>
#include "matrix.h"

int main(int argc, char **argv)
{
    struct matrix *mtrx_a;
    struct matrix *mtrx_b;

    FILE *mtrx_a_file = fopen("D:/c-programs/aads/lab_5/matrices/matrix_a", "r");
    FILE *mtrx_b_file = fopen("D:/c-programs/aads/lab_5/matrices/matrix_b", "r");
    
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
    
    struct matrix *mtrx_c = matrix_mul(mtrx_a, mtrx_b);

    matrix_print(mtrx_a);
    matrix_print(mtrx_b);

    matrix_print(mtrx_c);
    
    matrix_free(mtrx_a);
    matrix_free(mtrx_b);
    matrix_free(mtrx_c);

    return 0;
}
