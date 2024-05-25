#include <stdio.h>
#include "matrix.h"

int main(int argc, char **argv)
{
    struct matrix *mtrx_a;
    struct matrix *mtrx_b;

    FILE *mtrx_a_file = fopen("matrices/matrix_a", "r");
    FILE *mtrx_b_file = fopen("matrices/matrix_b", "r");
    
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

    //struct matrix *mtrx_c = matrix_get_column_vector(mtrx_a, 0);
    //struct matrix *housholder_c = matrix_housholder_vector(mtrx_c);

    //struct matrix *cut_vector = matrix_get_cut_column_vector(mtrx_a, 1, 1);

    printf("I am before housholder transformation\n");
    matrix_housholder_transformation(mtrx_a);
    matrix_print(mtrx_a);
    //matrix_print(mtrx_b);

    //matrix_print(cut_vector);
    
    //matrix_print(mtrx_c);
    //matrix_print(housholder_c);
    
    //matrix_free(mtrx_a);
    //matrix_free(mtrx_b);
    //matrix_free(mtrx_c);
    //matrix_free(housholder_c);

    return 0;
}
