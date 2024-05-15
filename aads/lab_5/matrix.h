#include <complex.h>

/*
 * Numbers stored column-wise to easily get
 * vectors for Hausholder method
 */
struct matrix {
    size_t rows;
    size_t columns;
    double _Complex *values;
};

int matrix_fread(FILE *mtrx_file, struct matrix **mtrx);
struct matrix *matrix_sum_seq(struct matrix *mtrx_1, struct matrix *mtrx_2);
struct matrix *matrix_sum_conc(struct matrix *mtrx_1, struct matrix *mtrx_2);
