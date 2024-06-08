#include <stdio.h>
#include <stdlib.h>
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


    printf("Old matrix B\n");

    matrix_print(mtrx_b);

    printf("\n");

    //struct matrix *houshold_matrix;

    //struct matrix *mtrx_c = matrix_cut_row_matrix(mtrx_a, 1, 2, 2, 3);
    
    printf("\n");

    matrix_hessenberg_triangular_transformation(&mtrx_a, mtrx_b);

    printf("\n"); 

    //struct matrix *housholder_c = matrix_housholder_vector(mtrx_c);

    //struct matrix *cut_vector = matrix_get_cut_column_vector(mtrx_a, 1, 1);

    //matrix_housholder_transformation(mtrx_a, &houshold_matrix);

    //matrix_print(mtrx_c);

    printf("New matrix B\n");

    matrix_print(mtrx_b);

    printf("\n");

    //matrix_print(houshold_matrix);
    
    //printf("\n");

    struct matrix *mtrx_b_inverse = matrix_get_triangular_inverse(mtrx_b);

    if (!mtrx_b_inverse) {
        printf("Matrix B is degenerate.\n");
        exit(1);
    }

    printf("Matrix B inverse\n");
    matrix_print(mtrx_b_inverse);
    printf("\n");   

    struct matrix *b_b_i = matrix_mul(mtrx_b, mtrx_b_inverse);
    printf("Is it really inverse:\n");
    matrix_print(b_b_i);
    printf("\n");
    //struct matrix *check = matrix_mul(houshold_matrix, mtrx_a);
    
    //printf("\n");

    //matrix_print(check);
    //matrix_print(cut_vector);
    
    //matrix_print(mtrx_c);
    //matrix_print(housholder_c);
    
    matrix_free(mtrx_a);
    matrix_free(mtrx_b);
    matrix_free(mtrx_b_inverse);
    matrix_free(b_b_i);
    //matrix_free(check);
    //matrix_free(houshold_matrix);

    exit(0);
}
