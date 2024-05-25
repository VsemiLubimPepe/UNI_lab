#include <complex.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"


struct matrix *matrix_init(size_t rows, size_t columns)
{
    struct matrix *mtrx = (struct matrix *)malloc(sizeof(struct matrix));
    if (!mtrx)
        return NULL;

    mtrx->values =
        (double _Complex *)malloc(sizeof(double _Complex) * rows * columns);
    if (!(mtrx->values)) {
        free(mtrx);
        return NULL;
    }

    mtrx->rows = rows;
    mtrx->columns = columns;

    return mtrx;
}

struct matrix *matrix_get_column_vector(struct matrix *mtrx, size_t column)
{
    size_t rows = mtrx->rows;

    struct matrix *vector = matrix_init(rows, 1);
    if (!vector)
	return NULL;

    memcpy(vector->values, mtrx->values + column * rows, sizeof(double _Complex) * rows);

    return vector;
}
/*
struct matrix *matrix_get_identity_matrix(size_t order)
{
    struct matrix *identity_matrix = matrix_init(order, order);
    if (!identity_matrix)
	return NULL;

    for (size_t i = 0; i < order; i++) {
	
    }
}
*/
struct matrix *matrix_get_cut_column_vector(struct matrix *mtrx, size_t column, size_t first_row)
{
    size_t rows = mtrx->rows;

    struct matrix *vector = matrix_init(rows - first_row, 1);
    if (!vector)
	return NULL;

    printf("I am before copying in geting cut vectro\n");
    memcpy(vector->values, mtrx->values + column * rows + first_row, sizeof(double _Complex) * (rows - first_row));

    return vector;
}

struct matrix *matrix_cut_matrix(struct matrix *mtrx, size_t row, size_t column)
{
    size_t rows = mtrx->rows;
    size_t columns = mtrx->columns;

    struct matrix *cut_mtrx = matrix_init(rows - row, columns - column);
    if (!cut_mtrx)
	return NULL;

    for (size_t i = 0; i < columns - column; i++) {
	memcpy(cut_mtrx->values + i * (rows - row), mtrx->values + (column + i) * rows + row, sizeof(double _Complex) * (rows - row));
    }

    return cut_mtrx;
}

int matrix_paste_matrix(struct matrix *mtrx_1, struct matrix *mtrx_2, size_t row, size_t column)
{
    size_t mtrx_1_rows = mtrx_1->rows;
    size_t mtrx_1_columns = mtrx_1->columns;
    size_t mtrx_2_rows = mtrx_2->rows;
    size_t mtrx_2_columns = mtrx_2->columns;

    if ((mtrx_1_rows - row != mtrx_2_rows) || (mtrx_1_columns - column != mtrx_2_columns))
	return -1;

    for (size_t i = 0; i < mtrx_2_columns; i++) {
	memcpy(mtrx_1->values + (column + i) * mtrx_1_rows + row, mtrx_2->values + i * mtrx_2_rows,
	       sizeof(double _Complex) * mtrx_2_rows);
    }

    return 0;
}

int matrix_fread(FILE *mtrx_file, struct matrix **mtrx)
{
    size_t rows;
    size_t columns;
    fscanf(mtrx_file, "%zu %zu ", &rows, &columns);
    
    *mtrx = matrix_init(rows, columns);
    if (!*mtrx)
	return -1;
    
    for (size_t i = 0; i < rows; i++) {
	for (size_t j = 0; j < columns; j++) {
	    double real;
	    double imag;
	    fscanf(mtrx_file, "%lf %lf*i", &real, &imag);
	    double _Complex read_number = real + imag * I;
	    ((*mtrx)->values)[i + j * rows] = read_number;
	}
    }
    
    return 0;
}

void matrix_print(struct matrix *mtrx)
{
    size_t rows = mtrx->rows;
    size_t columns = mtrx->columns;

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            double _Complex element = (mtrx->values)[i + j * rows];
            printf("%f%+f*i ", creal(element), cimag(element));
        }
        printf("\n");
    }

    return;
}

void matrix_free(struct matrix *mtrx)
{
    free(mtrx->values);
    free(mtrx);

    return;
}

struct matrix *matrix_sum(struct matrix *mtrx_1, struct matrix *mtrx_2)
{
    if ((mtrx_1->rows != mtrx_2->rows) || (mtrx_1->columns != mtrx_2->columns))
        return NULL;

    size_t rows = mtrx_1->rows;
    size_t columns = mtrx_1->columns;
 
    struct matrix *mtrx = matrix_init(rows, columns);
    if(!mtrx)
        return NULL;

    for (size_t j = 0; j < columns; j++) {
        for(size_t i = 0; i < rows; i++) {
            (mtrx->values)[i + j * columns] = 
                (mtrx_1->values)[i + j * columns] + (mtrx_2->values)[i + j * columns];
        }
    }

    return mtrx;
}

struct matrix *matrix_dif(struct matrix *mtrx_1, struct matrix *mtrx_2)
{
    if ((mtrx_1->rows != mtrx_2->rows) || (mtrx_1->columns != mtrx_2->columns))
        return NULL;

    size_t rows = mtrx_1->rows;
    size_t columns = mtrx_1->columns;

    struct matrix *mtrx = matrix_init(rows, columns);
    if(!mtrx)
        return NULL;

    for (size_t j = 0; j < columns; j++) {
        for(size_t i = 0; i < rows; i++) {
            (mtrx->values)[i + j * columns] =
                (mtrx_1->values)[i + j * columns] - (mtrx_2->values)[i + j * columns];
        }
    }

    return mtrx;
}

struct matrix *matrix_const_mul(struct matrix *mtrx, double _Complex number)
{
    size_t rows = mtrx->rows;
    size_t columns = mtrx->columns;

    struct matrix *new_mtrx = matrix_init(rows, columns);
    if (!new_mtrx)
	return NULL;
    
    for (size_t i = 0; i < rows; i++) {
	for (size_t j = 0; j < columns; j++) {
	    (new_mtrx->values)[j + i * columns] = number * (mtrx->values)[j + i * columns];
	}
    }

    return new_mtrx;
}

struct matrix *matrix_mul(struct matrix *mtrx_1, struct matrix *mtrx_2)
{
    if (mtrx_1->columns != mtrx_2->rows)
        return NULL;

    size_t rows = mtrx_1->rows;
    size_t columns = mtrx_2->columns;
    size_t common_dimension = mtrx_1->columns;
    struct matrix* mtrx = matrix_init(rows, columns);

    if (!mtrx)
        return NULL;

    for (size_t j = 0; j < columns; j++) {
        for (size_t i = 0; i < rows; i++) {
            (mtrx->values)[i + j * columns] = 0.0 + 0.0*I;
            for (size_t n = 0; n < common_dimension; n++) {
                (mtrx->values)[i + j * columns] +=
                    (mtrx_1->values)[i + n * rows] * 
                    (mtrx_2->values)[n + j * common_dimension];
            }
        }
    }

    return mtrx;
}

struct matrix *matrix_transpose(struct matrix *mtrx)
{
    size_t rows = mtrx->rows;
    size_t columns = mtrx->columns;

    struct matrix *trans_mtrx = matrix_init(columns, rows);
    if (!trans_mtrx)
        return NULL;
    
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            (trans_mtrx->values)[i * columns + j] = 
            (mtrx->values)[i + j * rows];
        }
    }

    return trans_mtrx;
}

struct givens_pair *givens_pair_init(double _Complex c, double _Complex s)
{
    return (struct givens_pair *)malloc(sizeof(struct givens_pair));
}

struct givens_pair *givens_pair_find(double _Complex a, double _Complex b)
{
    if ((creal(b) == 0.0) && (cimag(b) == 0.0))
	return givens_pair_init(1 + 0 * I, 0 + 0 * I);

    double _Complex tau;
    double _Complex s;
    double _Complex c;
    
    if (cabs(b) > cabs(a)) {
	tau = -a / b;
	s = 1 / csqrt(1 + tau * tau);
	c = s * tau;
    }
    else {
	tau = -b / a;
	c = 1 / csqrt(1 + tau * tau);
	s = c * tau;
    }

    return givens_pair_init(c, s);
}

struct matrix *matrix_givens_left_mul(struct matrix *mtrx, struct givens_pair *givens_pair)
{
    size_t rows = mtrx->rows;
    size_t columns = mtrx->columns;

    struct matrix *rotated_mtrx = matrix_init(rows, columns);
    if (!rotated_mtrx)
	return NULL;

    memcpy(rotated_mtrx->values, mtrx->values,
	   sizeof(double _Complex) * rows * columns);

    for (size_t j = 0; j < columns; j++) {
	double _Complex tau_1 = (rotated_mtrx->values)[j * columns];
	double _Complex tau_2 = (rotated_mtrx->values)[j * columns + 1];
	(rotated_mtrx->values)[j * columns] =
	    givens_pair->c * tau_1 - givens_pair->s * tau_2;
	(rotated_mtrx->values)[j * columns + 1] =
	    givens_pair->s * tau_1 + givens_pair->c * tau_2;
    }

    return rotated_mtrx;
}

struct matrix *matrix_givens_right_mul(struct matrix *mtrx, struct givens_pair *givens_pair)
{
    size_t rows = mtrx->rows;
    size_t columns = mtrx->columns;

    struct matrix *rotated_mtrx = matrix_init(rows, columns);
    if (!rotated_mtrx)
	return NULL;

    memcpy(rotated_mtrx->values, mtrx->values,
	   sizeof(double _Complex) * rows * columns);

    for (size_t i = 0; i < rows; i++) {
	double _Complex tau_1 = (rotated_mtrx->values)[i];
	double _Complex tau_2 = (rotated_mtrx->values)[i + rows];
	(rotated_mtrx->values)[i] =
	    givens_pair->c * tau_1 - givens_pair->s * tau_2;
	(rotated_mtrx->values)[i + rows] =
	    givens_pair->s * tau_1 + givens_pair->c * tau_2;
    }

    return rotated_mtrx;
}

int  matrix_vector_norm(struct matrix *vector, double _Complex *norm)
{
    struct matrix *transposed_vector = matrix_transpose(vector);
    if (!transposed_vector)
	return -1;

    struct matrix *norm_matrix = matrix_mul(transposed_vector, vector);
    if (!norm_matrix) {
	matrix_free(transposed_vector);
	return -1;
    }

    *norm = csqrt((norm_matrix->values)[0]);
    
    matrix_free(transposed_vector);
    matrix_free(norm_matrix);
    
    return 0;
}


struct matrix *matrix_housholder_vector(struct matrix *vector)
{
    struct matrix *housholder_vector = matrix_init(vector->rows, vector->columns);
    
    size_t vector_length;
    if (vector->rows > vector->columns)
	vector_length = vector->rows;
    else
	vector_length = vector->columns;

    double _Complex mu;

    if (matrix_vector_norm(vector, &mu))
	return NULL;

    if (mu != 0.0 + 0.0 * I) {
	double _Complex beta = (vector->values)[0];
	beta += (vector->values)[0] / cabs((vector->values)[0]) * mu;
	
	for (size_t i = 1; i < vector_length; i++)
	    (housholder_vector->values)[i] = (vector->values)[i] / beta; 
    }
    (housholder_vector->values)[0] = 1;

    return housholder_vector;
}

// row.house
struct matrix *matrix_housholder_left_mul(struct matrix *mtrx, struct matrix *housholder_vector)
{
    struct matrix *transposed_housholder_vector = matrix_transpose(housholder_vector);
    if (!transposed_housholder_vector)
	return NULL;
    
    struct matrix *hvt_hv_matrix = matrix_mul(transposed_housholder_vector, housholder_vector);
    if (!hvt_hv_matrix) {
	matrix_free(transposed_housholder_vector);
	return NULL;
    }

    double _Complex hvt_hv = (hvt_hv_matrix->values)[0];
    matrix_free(hvt_hv_matrix);
    matrix_free(transposed_housholder_vector);

    double _Complex beta = -2 / hvt_hv;
    
    struct matrix *transposed_matrix = matrix_transpose(mtrx);
    if (!transposed_matrix)
	return NULL;
    
    struct matrix *omega = matrix_mul(transposed_matrix, housholder_vector);
    omega = matrix_const_mul(omega, beta);
    
    matrix_free(transposed_matrix);

    struct matrix *transposed_omega = matrix_transpose(omega);
    if (!transposed_omega) {
	matrix_free(omega);
	return NULL;
    }

    matrix_free(omega);
    
    struct matrix *v_omegat = matrix_mul(housholder_vector, transposed_omega);
    if (!v_omegat) {
	matrix_free(transposed_omega);
	return NULL;
    }

    struct matrix *result_matrix = matrix_sum(mtrx, v_omegat);
    if (!result_matrix) {
	matrix_free(v_omegat);
	return NULL;
    }

    matrix_free(v_omegat);

    return result_matrix;
}

// col.house
struct matrix *matrix_housholder_right_mul(struct matrix *mtrx, struct matrix *housholder_vector)
{
    struct matrix *transposed_housholder_vector = matrix_transpose(housholder_vector);
    if (!transposed_housholder_vector)
	return NULL;
    
    struct matrix *hvt_hv_matrix = matrix_mul(transposed_housholder_vector, housholder_vector);
    if (!hvt_hv_matrix) {
	matrix_free(transposed_housholder_vector);
	return NULL;
    }

    double _Complex hvt_hv = (hvt_hv_matrix->values)[0];
    matrix_free(hvt_hv_matrix);

    double _Complex beta = -2 / hvt_hv;
    
    struct matrix *omega = matrix_mul(mtrx, housholder_vector);
    omega = matrix_const_mul(omega, beta);
    
    struct matrix *omega_vt = matrix_mul(omega, transposed_housholder_vector);
    if (!omega_vt) {
	matrix_free(omega);
	matrix_free(transposed_housholder_vector);
	return NULL;
    }

    matrix_free(omega);
    matrix_free(transposed_housholder_vector);

    struct matrix *result_matrix = matrix_sum(mtrx, omega_vt);
    if (!result_matrix) {
	matrix_free(omega_vt);
	return NULL;
    }

    matrix_free(omega_vt);

    return result_matrix;
}


// Add creating copy of initial matrix
void matrix_housholder_transformation(struct matrix *mtrx)
{
    size_t columns = mtrx->columns;
    size_t rows = mtrx->rows;
    
    
    for (size_t i = 0; i < rows - 1; i++) {
	struct matrix *vector = matrix_get_cut_column_vector(mtrx, i, i);
	if (!vector)
	    printf("Vector is NULL.\n");
	struct matrix *housholder_vector = matrix_housholder_vector(vector);
	if (!housholder_vector)
	    printf("Housholder vector is NULL.\n");
	matrix_free(vector);
	struct matrix *matrix_cut = matrix_cut_matrix(mtrx, i, i);
	if (!matrix_cut)
	    printf("Matrix cut is NULL.\n");
	struct matrix *cut_hous_mul_left_matrix = matrix_housholder_left_mul(mtrx, housholder_vector);
	if (!cut_hous_mul_left_matrix)
	    printf("cut_house_mul_left_matrix is NULL.\n");
	matrix_free(matrix_cut);
	matrix_paste_matrix(mtrx, cut_hous_mul_left_matrix, i, i);
	matrix_free(cut_hous_mul_left_matrix);
    }
}

