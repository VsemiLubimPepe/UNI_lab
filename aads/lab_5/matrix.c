#include <complex.h>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"
#include "thresholds.h"

/*
 * Zero threshold for complex number parts.
 */
double threshold;
/*
 * Threshold multiplier for values of matrix A subdiagonal elements.
 */
double epsilon;

/*
 * Value of zero in complex form.
 */
const double _Complex zero = 0.0 + 0.0 * I;

/*
 * Prints complex number NUM.
 */
void complex_print(double _Complex num)
{
    printf("%f%+f*i\n", creal(num), cimag(num));
}

/*
 * Prepares empty matrix of size ROWS x COLUMNS.
 */
struct matrix *matrix_init(size_t rows, size_t columns)
{
    // Allocating memory for matrix struct.
    struct matrix *mtrx = malloc(sizeof(struct matrix));
    if (!mtrx)
        return NULL;

    // Allocating memory for matrix elements.
    mtrx->values = malloc(sizeof(double _Complex) * rows * columns);
    if (!(mtrx->values)) {
        free(mtrx);
        return NULL;
    }

    // Storing info about matrix size.
    mtrx->rows = rows;
    mtrx->columns = columns;

    return mtrx;
}

/*
 * Returns identity matrix of ORDER order (allocates new matrix).
 */
struct matrix *matrix_get_identity_matrix(size_t order)
{
    // Preparing empty matrix.
    struct matrix *identity_matrix = matrix_init(order, order);
    if (!identity_matrix)
	    return NULL;

    // Making identity matrix.
    for (size_t j = 0; j < order; j++)
        for (size_t i = 0; i < order; i++) {
            if (i == j)
                (identity_matrix->values)[j + i * order] = 1.0 + 0.0 * I;
            else
                (identity_matrix->values)[j + i * order] = zero;
        }

    return identity_matrix;
}

/*
 * Return zero matrix of order ORDER.
 */
struct matrix *matrix_get_zeros(size_t order)
{
    // Preparing empty matrix.
    struct matrix *zeros_matrix = matrix_init(order, order);
    if (!zeros_matrix)
        return NULL;

    // Making zero matrix.
    for (size_t i = 0; i < order; i++)
        for (size_t j = 0; j < order; j++)
            (zeros_matrix->values)[j + i * order] = zero;
    
    return zeros_matrix;
}

/*
 * Returns submatrix starting from element (ROW, COLUMN) 
 * and ending on element (END_ROW, END_COLUMN).
 */
struct matrix *matrix_cut_matrix(struct matrix *mtrx, size_t row, 
    size_t end_row, size_t column, size_t end_column)
{
    // Rows of MTRX
    size_t rows = mtrx->rows;

    // Rows of cut matrix.
    size_t cut_rows = end_row - row + 1;
    // Columns of cut matrix.
    size_t cut_columns = end_column - column + 1;

    // Preparing empty cut matrix.
    struct matrix *cut_mtrx = matrix_init(cut_rows, cut_columns);  
    if(!cut_mtrx)
        return NULL;

    // Filling cut matrix with needed values of MTRX.
    for (size_t j = 0; j < cut_columns; j++)
        for (size_t i = 0; i < cut_rows; i++)
            (cut_mtrx->values)[j * (cut_rows) + i] = 
                (mtrx->values)[row + i + (column + j) * rows];

    return cut_mtrx;
}

/*
 * Inserts submatrix MTRX_2 into matrix MTRX_1 starting 
 * from element (ROW, COLUMN).
 */
int matrix_paste_matrix(struct matrix *mtrx_1, struct matrix *mtrx_2, 
    size_t row, size_t column)
{
    // Rows of target matrix.
    size_t mtrx_1_rows = mtrx_1->rows;
    // Columns of target matrix.
    size_t mtrx_1_columns = mtrx_1->columns;
    // Rows of paste matrix.
    size_t mtrx_2_rows = mtrx_2->rows;
    // Columns of paste matrix
    size_t mtrx_2_columns = mtrx_2->columns;

    // Checking if paste matrix would fit in designated place of target matrix.
    if ((mtrx_1_rows - row < mtrx_2_rows) || 
        (mtrx_1_columns - column < mtrx_2_columns))
	    return -1;

    // Filling values from paste matrix in target matrix.
    for (size_t j = 0; j < mtrx_2_columns; j++)
	    for (size_t i = 0; i < mtrx_2_rows; i++)
            (mtrx_1->values)[row + i + (column + j) * mtrx_1_rows] = 
                (mtrx_2->values)[j * mtrx_2_rows + i];

    return 0;
}

/*
 * Reads matrix from a file.
 * File format:
 * ROWS COLUMNS
 * MATRIX_ROW
 * ...
 * MATRIX_ROW
 */
int matrix_fread(FILE *mtrx_file, struct matrix **mtrx)
{
    // Rows of read matrix.
    size_t rows;
    // Columns of read matrix.
    size_t columns;

    // Checking if ROWS and COLUMNS were read successfully.
    if (fscanf(mtrx_file, "%zu %zu", &rows, &columns) != 2)
        return -1;
    
    // Preparing empty matrix.
    *mtrx = matrix_init(rows, columns);
    if (!*mtrx)
	    return -1;
    
    // Reading values from file into matrix.
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < columns; j++) {
            double real;
            double imag;
            if (fscanf(mtrx_file, "%lf %lf*i", &real, &imag) != 2)
                return -1;
            double _Complex read_number = real + imag * I;
            ((*mtrx)->values)[i + j * rows] = read_number;
        }
    
    return 0;
}

/*
 * Frees allocated memory of the matrix MTRX.
 */
void matrix_free(struct matrix *mtrx)
{
    // Freeing memory allocated for matrix elements.
    free(mtrx->values);
    // Freeing memory allocated for matrix struct.
    free(mtrx);
}

/*
 * Sets zeros in parts of complex number, if they are close to zero.
 */
double _Complex complex_threshold_cut(double _Complex number)
{
    // Real part of NUMBER.
    double new_real = creal(number);
    // Imaginary part of NUMBER.
    double new_imag = cimag(number);
    
    // Checking if real part of NUMBER should be zeroed.
    if ((new_real < threshold) && (new_real > -threshold))
        new_real = 0.0;
    // Checking if imaginary part of NUMBER should be zeroed.
    if ((new_imag < threshold) && (new_imag > -threshold))
        new_imag = 0.0;

    return new_real + new_imag * I; 
}

/*
 * Sets zeros in parts of complex numbers of matrix, if they are close to zero.
 */
void matrix_threshold_cut(struct matrix *mtrx)
{
    // Rows of MTRX.
    size_t rows = mtrx->rows;
    // Columns of MTRX.
    size_t columns = mtrx->columns;

    // Handling of all elements of MTRX.
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < columns; j++)
            mtrx->values[i + j * rows] = 
                complex_threshold_cut(mtrx->values[i + j * rows]);
}

/*
 * Returns the matrix reflected relatively to the first (or last) 
 * diagonal element of MTRX (allocates memory).
 */
struct matrix *matrix_dot_mirror(struct matrix *mtrx)
{
    // Rows of MTRX.
    size_t rows = mtrx->rows;
    // Columns of MTRX.
    size_t columns = mtrx->columns;

    // Preparing empty mirrored matrix.
    struct matrix *mirror = matrix_init(rows, columns);

    // Filling mirrored matrix.
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < columns; j++)
            (mirror->values)[i + j * rows] = 
                conj((mtrx->values)[rows - 1 - i + (columns - 1 - j) * rows]);

    return mirror;
}

/*
 * Returns the submatrix from matrix MTRX for finding minor (allocates memory).
 */
struct matrix *matrix_get_minor_matrix(struct matrix *mtrx, size_t row, 
    size_t column)
{
    // Order of MTRX.
    size_t order = mtrx->rows;
    // Row which should be omited.
    size_t row_counter = 0;
    // Column which should be omited.
    size_t column_counter = 0;

    // Preparing empty minor matrix.
    struct matrix *minor_matrix = matrix_init(order - 1, order - 1);

    // Filling minor matrix.
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < column; j++) {
            (minor_matrix->values)[row_counter + column_counter * (order - 1)] =
                (mtrx->values)[i + j * order];
            column_counter++;
        }
        for (size_t j = column + 1; j < order; j++) {
            (minor_matrix->values)[row_counter + column_counter * (order - 1)] =
                (mtrx->values)[i + j * order];
            column_counter++;
        }
        column_counter = 0;
        row_counter++;
    }
    for (size_t i = row + 1; i < order; i++) {
        for (size_t j = 0; j < column; j++) {
            (minor_matrix->values)[row_counter + column_counter * (order - 1)] =
                (mtrx->values)[i + j * order];
            column_counter++;
        }
        for (size_t j = column + 1; j < order; j++) {
            (minor_matrix->values)[row_counter + column_counter * (order - 1)] =
                (mtrx->values)[i + j * order];
            column_counter++;
        }
        column_counter = 0;
        row_counter++;
    }

    return minor_matrix;
}

/*
 * Returns the sum of argument matrices (allocates new matrix).
 */
struct matrix *matrix_sum(struct matrix *mtrx_1, struct matrix *mtrx_2)
{
    // Checking if sizes of matrices are equal.
    if ((mtrx_1->rows != mtrx_2->rows) || (mtrx_1->columns != mtrx_2->columns))
        return NULL;

    // Rows of matrices.
    size_t rows = mtrx_1->rows;
    // Columns of matrices.
    size_t columns = mtrx_1->columns;
    
    // Preparing empty sum matrix.
    struct matrix *mtrx = matrix_init(rows, columns);
    if(!mtrx)
        return NULL;

    // Filling sum matrix.
    for (size_t j = 0; j < columns; j++)
        for(size_t i = 0; i < rows; i++)
            (mtrx->values)[i + j * rows] = 
                (mtrx_1->values)[i + j * rows] + 
                (mtrx_2->values)[i + j * rows];

    // Handling all values which are close to zero.
    matrix_threshold_cut(mtrx);

    return mtrx;
}

/*
 * Returns the difference of argument matrices (allocates new matrix).
 */
struct matrix *matrix_dif(struct matrix *mtrx_1, struct matrix *mtrx_2)
{
    // Checking if sizes of matrices are equal.
    if ((mtrx_1->rows != mtrx_2->rows) || (mtrx_1->columns != mtrx_2->columns))
        return NULL;

    // Rows of matrices.
    size_t rows = mtrx_1->rows;
    // Columns of matrices.
    size_t columns = mtrx_1->columns;

    // Preparing empty difference matrix.
    struct matrix *mtrx = matrix_init(rows, columns);
    if(!mtrx)
        return NULL;

    // Filling difference matrix.
    for (size_t j = 0; j < columns; j++)
        for(size_t i = 0; i < rows; i++)
            (mtrx->values)[i + j * rows] =
                (mtrx_1->values)[i + j * rows] - 
                (mtrx_2->values)[i + j * rows];

    // Handling all values which are close to zero.
    matrix_threshold_cut(mtrx);

    return mtrx;
}

/*
 * Returns the product of MTRX and CONST value (allocates memory).
 */
struct matrix *matrix_const_mul(struct matrix *mtrx, double _Complex number)
{
    // Rows of matrix.
    size_t rows = mtrx->rows;
    // Columns of matrix.
    size_t columns = mtrx->columns;

    // Preparing empty product matrix.
    struct matrix *new_mtrx = matrix_init(rows, columns);
    if (!new_mtrx)
	    return NULL;
    
    // Filling product matrix.
    for (size_t j = 0; j < columns; j++)
        for (size_t i = 0; i < rows; i++)
            (new_mtrx->values)[j * rows + i] = 
                number * (mtrx->values)[j * rows + i];

    // Handling all values which are close to zero.
    matrix_threshold_cut(new_mtrx);

    return new_mtrx;
}

/*
 * Returns the product of argument matrices (allocates new matrix).
 */
struct matrix *matrix_mul(struct matrix *mtrx_1, struct matrix *mtrx_2)
{
    // Checking if matrices can be multiplied.
    if (mtrx_1->columns != mtrx_2->rows)
        return NULL;

    size_t rows = mtrx_1->rows;
    size_t columns = mtrx_2->columns;
    size_t common_dimension = mtrx_1->columns;
    
    // Preparing empty product matrix.
    struct matrix* mtrx = matrix_init(rows, columns);
    if (!mtrx)
        return NULL;
    
    // Filling product matrix.
    for (size_t j = 0; j < columns; j++)
        for (size_t i = 0; i < rows; i++) {
            (mtrx->values)[i + j * rows] = zero;
            for (size_t n = 0; n < common_dimension; n++) {
                (mtrx->values)[i + j * rows] +=
                    (mtrx_1->values)[i + n * rows] * 
                    (mtrx_2->values)[n + j * common_dimension];
            }
        }

    // Handling all values which are close to zero.
    matrix_threshold_cut(mtrx);

    return mtrx;
}

/*
 * Returns Hermitian transpose of MTRX.
 */
struct matrix *matrix_hermitian_transpose(struct matrix *mtrx)
{
    size_t rows = mtrx->rows;
    size_t columns = mtrx->columns;

    // Preparing empty Hermitian transpose of MTRX
    struct matrix *trans_mtrx = matrix_init(columns, rows);
    if (!trans_mtrx)
        return NULL;
    
    // Filling Hermitian transpose of MTRX
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < columns; j++)
            (trans_mtrx->values)[i * columns + j] = 
                conj((mtrx->values)[i + j * rows]);

    return trans_mtrx;
}

/*
 * Returns givens pair (allocates new givens pair).
 */
struct givens_pair *givens_pair_init(double _Complex c, double _Complex s)
{
    // Allocating memory for givens pair struct.
    struct givens_pair *givens_pair = malloc(sizeof(struct givens_pair));
    givens_pair->c = c;
    givens_pair->s = s;
    return givens_pair;
}

/*
 * Returns Givens pair for Givens rotation (allocates memory).
 */
struct givens_pair *givens_pair_find(double _Complex a, double _Complex b)
{
    if (b == zero)
	    return givens_pair_init(1 + 0 * I, 0 + 0 * I);

    double _Complex tau;
    // Elements of givens pair.
    double _Complex s;
    double _Complex c;
    
    // Choosing such givens pair elements that there would be no precision loss.
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

    // Handling values which are close to zero.
    c = complex_threshold_cut(c);
    s = complex_threshold_cut(s);

    return givens_pair_init(c, s);
}

/*
 * Returns the product of Givens matrix and MTRX (GM * MTRX) (allocates memory).
 */
struct matrix *matrix_givens_left_mul(struct matrix *mtrx, 
    struct givens_pair *givens_pair)
{
    size_t rows = mtrx->rows;
    size_t columns = mtrx->columns;

    // Copying MTRX.
    struct matrix *rotated_mtrx = 
        matrix_cut_matrix(mtrx, 0, rows - 1, 0, columns - 1);
    if (!rotated_mtrx)
	    return NULL;

    // Changing values of MTRX which were affected by rotation.
    for (size_t j = 0; j < columns; j++) {
        double _Complex tau_1 = (rotated_mtrx->values)[j * rows];
        double _Complex tau_2 = (rotated_mtrx->values)[j * rows + 1];
        (rotated_mtrx->values)[j * rows] =
            givens_pair->c * tau_1 - givens_pair->s * tau_2;
        (rotated_mtrx->values)[j * rows + 1] =
            givens_pair->s * tau_1 + givens_pair->c * tau_2;
    }

    // Handling elements which are close to zero.
    matrix_threshold_cut(rotated_mtrx);

    return rotated_mtrx;
}

/*
 * Returns the product of MTRX and Givens matrix (MTRX * GM) (allocates memory).
 */
struct matrix *matrix_givens_right_mul(struct matrix *mtrx, 
    struct givens_pair *givens_pair)
{
    size_t rows = mtrx->rows;
    size_t columns = mtrx->columns;

    // Copying MTRX.
    struct matrix *rotated_mtrx = 
        matrix_cut_matrix(mtrx, 0, rows - 1, 0, columns - 1);
    if (!rotated_mtrx)
	    return NULL;

    // Changing values of MTRX which were affected by rotation.
    for (size_t i = 0; i < rows; i++) {
        double _Complex tau_1 = (rotated_mtrx->values)[i];
        double _Complex tau_2 = (rotated_mtrx->values)[i + rows];
        (rotated_mtrx->values)[i + rows] =
            givens_pair->c * tau_1 - givens_pair->s * tau_2;
        (rotated_mtrx->values)[i] =
            givens_pair->s * tau_1 + givens_pair->c * tau_2;
    }

    // Handling elements which are close to zero.
    matrix_threshold_cut(rotated_mtrx);

    return rotated_mtrx;
}

/*
 * Returns the euclidean norm of VECTOR.
 */
void matrix_vector_norm(struct matrix *vector, double _Complex *norm)
{
    size_t rows = vector->rows;

    double _Complex found_norm = 0;

    // Calculating VECTOR norm.
    for (size_t i = 0; i < rows; i++)
        found_norm += (vector->values)[i] * conj((vector->values)[i]);

    *norm = csqrt(found_norm);
}

/*
 * Returns the Householder vector of VECTOR (allocates memory).
 */
struct matrix *matrix_housholder_vector(struct matrix *vector)
{
    // Preparing empty housholder vector.
    struct matrix *housholder_vector = 
        matrix_init(vector->rows, vector->columns);

    size_t vector_length = vector->rows;
    
    // Variable for VECTOR norm.
    double _Complex mu;

    // Calculating VECTOR norm.
    matrix_vector_norm(vector, &mu);

    // Calculating VECTOR norm.
    if (mu != zero) {
        double _Complex beta = (vector->values)[0];
        if (beta != zero)
            beta += (vector->values)[0] / cabs((vector->values)[0]) * mu;
        else 
            beta += mu;

        for (size_t i = 1; i < vector_length; i++)
            (housholder_vector->values)[i] = (vector->values)[i] / beta; 
    }
    (housholder_vector->values)[0] = 1;

    // Handling all elements close to zero.
    matrix_threshold_cut(housholder_vector);

    return housholder_vector;
}

/*
 * Returns Housholder reflection matrix (allocates memory).
 */
struct matrix *matrix_housholder_vector_matrix(struct matrix *vector)
{
    size_t length = vector->rows;
    
    struct matrix *housholder_vector = matrix_housholder_vector(vector);
    struct matrix *transposed_hous_vec = 
        matrix_hermitian_transpose(housholder_vector);
    struct matrix *hv_hvt = matrix_mul(housholder_vector, transposed_hous_vec);
    struct matrix *hvt_hv = matrix_mul(transposed_hous_vec, housholder_vector);
    
    matrix_free(housholder_vector);
    matrix_free(transposed_hous_vec);
    
    struct matrix *hv_hvt_2 = matrix_const_mul(hv_hvt, 2);
    struct matrix *fin = matrix_const_mul(hv_hvt_2, 1 / hvt_hv->values[0]);
    
    matrix_free(hv_hvt);
    matrix_free(hvt_hv);
    matrix_free(hv_hvt_2);

    struct matrix *identity = matrix_get_identity_matrix(length);
    struct matrix *hous_mtrx = matrix_dif(identity, fin);
    
    matrix_free(fin);
    matrix_free(identity);

    return hous_mtrx;
}

/*
 * Returns the product of Housholder matrix and MTRX 
 * (HM * MTRX) (allocates memory).
 */
struct matrix *matrix_housholder_left_mul(struct matrix *mtrx, 
    struct matrix *housholder_vector)
{
    struct matrix *transposed_housholder_vector = 
        matrix_hermitian_transpose(housholder_vector);
    if (!transposed_housholder_vector)
	    return NULL;

    struct matrix *hvt_hv_matrix = matrix_mul(transposed_housholder_vector, 
        housholder_vector);
    if (!hvt_hv_matrix) {
        matrix_free(transposed_housholder_vector);
        return NULL;
    }

    double _Complex hvt_hv = (hvt_hv_matrix->values)[0];
    double _Complex beta = -2 / hvt_hv;

    matrix_free(hvt_hv_matrix);
    matrix_free(transposed_housholder_vector);
    
    struct matrix *transposed_matrix = matrix_hermitian_transpose(mtrx);
    if (!transposed_matrix)
	    return NULL;
    struct matrix *pre_omega = matrix_mul(transposed_matrix, housholder_vector);
    struct matrix *omega = matrix_const_mul(pre_omega, beta);
    
    matrix_free(pre_omega);
    matrix_free(transposed_matrix);

    struct matrix *transposed_omega = matrix_hermitian_transpose(omega);
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

    matrix_free(transposed_omega);

    struct matrix *result_matrix = matrix_sum(mtrx, v_omegat);
    if (!result_matrix) {
        matrix_free(v_omegat);
        return NULL;
    }

    matrix_free(v_omegat);

    return result_matrix;
}

/*
 * Returns the product of MTRX and Householder matrix 
 * (MTRX * HM) (allocates memory).
 */
struct matrix *matrix_housholder_right_mul(struct matrix *mtrx, 
    struct matrix *housholder_vector)
{
    struct matrix *transposed_housholder_vector = 
        matrix_hermitian_transpose(housholder_vector);
    if (!transposed_housholder_vector)
	    return NULL;

    struct matrix *hvt_hv_matrix = 
        matrix_mul(transposed_housholder_vector, housholder_vector);
    if (!hvt_hv_matrix) {
        matrix_free(transposed_housholder_vector);
        return NULL;
    }

    double _Complex hvt_hv = (hvt_hv_matrix->values)[0];
    matrix_free(hvt_hv_matrix);

    double _Complex beta = -2 / hvt_hv;
    
    struct matrix *pre_omega = matrix_mul(mtrx, housholder_vector);
    struct matrix *omega = matrix_const_mul(pre_omega, beta);

    matrix_free(pre_omega);
    
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

/*
 * Returns upper triangular form of the MTRX (using Householder matrices).
 */
void matrix_housholder_transformation(struct matrix *mtrx, 
    struct matrix **houshold_matrix)
{
    size_t columns = mtrx->columns;
    size_t rows = mtrx->rows;
    
    // Array for all found householder vectors.
    struct matrix *housholder_vectors[columns];
    
    // Transforming MTRX.
    for (size_t i = 0; i < columns; i++) {
	    struct matrix *vector = matrix_cut_matrix(mtrx, i, rows - 1, i, i);
        struct matrix *housholder_vector = matrix_housholder_vector(vector);
        housholder_vectors[i] = housholder_vector;
        struct matrix *matrix_cut = 
            matrix_cut_matrix(mtrx, i, rows - 1, i, columns - 1);
        struct matrix *cut_hous_mul_left_matrix = 
            matrix_housholder_left_mul(matrix_cut, housholder_vector);
        matrix_free(matrix_cut);
        matrix_free(vector);
        matrix_paste_matrix(mtrx, cut_hous_mul_left_matrix, i, i);
        matrix_free(cut_hous_mul_left_matrix);
    }
    
    // Matrix - product of all found Householder matrices.
    *houshold_matrix = matrix_get_identity_matrix(rows);
    
    for (size_t i = 0; i < columns; i++) {
        struct matrix *matrix_cut = 
            matrix_cut_matrix(*houshold_matrix, rows - 1 - i, rows - 1,
                columns - 1 - i, columns - 1);
            
        struct matrix *cut_hous_mul_left_matrix = 
            matrix_housholder_left_mul(matrix_cut, 
                housholder_vectors[columns - 1 - i]);

        matrix_free(matrix_cut);

        matrix_paste_matrix(*houshold_matrix, cut_hous_mul_left_matrix, 
            columns - 1 - i, columns - 1 - i);
        
        matrix_free(cut_hous_mul_left_matrix);
        matrix_free(housholder_vectors[columns - 1 - i]);
    }
}

/*
 * Returns MTRX_A in upper Hessenberg form and MTRX_A in upper triangular 
 * form (using Givens rotations and Housholder transformation). 
 */
void matrix_hessenberg_triangular_transformation(struct matrix **mtrx_a, 
    struct matrix *mtrx_b)
{
    // Transforming MTRX_B to triangular form.
    struct matrix *housholder_matrix;
    matrix_housholder_transformation(mtrx_b, &housholder_matrix);

    struct matrix *transposed_housholder_matrix = 
        matrix_hermitian_transpose(housholder_matrix);
    matrix_free(housholder_matrix);

    // Applying MTRX_B transformation to MTRX_A.
    struct matrix *new_mtrx_a = 
        matrix_mul(transposed_housholder_matrix, *mtrx_a);
    matrix_free(transposed_housholder_matrix);

    size_t order = new_mtrx_a->rows;

    // Transforming MTRX_A to Hessenberg form.
    for (size_t j = 0; j < order - 2; j++)
        for (size_t i = order - 1; i > j + 1; i--) {
            struct givens_pair *givens_pair = 
                givens_pair_find((new_mtrx_a->values)[i - 1 + j * order], 
                    (new_mtrx_a->values)[i + j * order]);

            struct matrix *cut_mtrx_a = 
                matrix_cut_matrix(new_mtrx_a, i - 1, i, j, order - 1);
            struct matrix *rot_cut_mtrx_a = 
                matrix_givens_left_mul(cut_mtrx_a, givens_pair);
            matrix_paste_matrix(new_mtrx_a, rot_cut_mtrx_a, i-1, j);

            matrix_free(rot_cut_mtrx_a);
            matrix_free(cut_mtrx_a);

            struct matrix *cut_mtrx_b = 
                matrix_cut_matrix(mtrx_b, i - 1, i, i - 1, order - 1);
            struct matrix *rot_cut_mtrx_b = 
                matrix_givens_left_mul(cut_mtrx_b, givens_pair);
            matrix_paste_matrix(mtrx_b, rot_cut_mtrx_b, i - 1, i - 1);

            matrix_free(rot_cut_mtrx_b);
            matrix_free(cut_mtrx_b);
            free(givens_pair);

            givens_pair = 
                givens_pair_find((mtrx_b->values)[i + (i - 1) * order], 
                    (mtrx_b->values)[i + i * order]);

            cut_mtrx_a = 
                matrix_cut_matrix(new_mtrx_a, 0, order - 1, i - 1, i);
            rot_cut_mtrx_a = matrix_givens_right_mul(cut_mtrx_a, givens_pair);
            matrix_paste_matrix(new_mtrx_a, rot_cut_mtrx_a, 0, i - 1);

            matrix_free(rot_cut_mtrx_a);
            matrix_free(cut_mtrx_a);

            cut_mtrx_b = matrix_cut_matrix(mtrx_b, 0, i, i - 1, i);
            rot_cut_mtrx_b = matrix_givens_right_mul(cut_mtrx_b, givens_pair);
            matrix_paste_matrix(mtrx_b, rot_cut_mtrx_b, 0, i - 1);

            matrix_free(rot_cut_mtrx_b);
            matrix_free(cut_mtrx_b);
            free(givens_pair);
        }

    matrix_free(*mtrx_a);
    *mtrx_a = new_mtrx_a;
}

/*
 * Returns eigenvalues of second order MTRX (allocates memory).
 */
struct matrix *matrix_get_second_order_eigenvalues(struct matrix *mtrx)
{
    size_t order = mtrx->rows;
    struct matrix *eigenvalues = matrix_init(order, 1);

    double _Complex a = 1.0 + 0.0 * I;
    double _Complex b = - (mtrx->values[0] + mtrx->values[3]);
    double _Complex c = mtrx->values[0] * mtrx->values[3] - mtrx->values[1] * 
        mtrx->values[2];

    b = complex_threshold_cut(b);
    c = complex_threshold_cut(c);

    double _Complex D = b * b - 4 * a * c;
    D = complex_threshold_cut(D);

    if (D == 0) {
        double _Complex lambda = - b / (2 * a);
        eigenvalues->values[0] = lambda;
        eigenvalues->values[1] = lambda;
    }
    else {
        double _Complex lambda_1 = (- b + csqrt(D)) / (2 * a);
        double _Complex lambda_2 = (- b - csqrt(D)) / (2 * a);
        eigenvalues->values[0] = lambda_2;
        eigenvalues->values[1] = lambda_1; 
    }

    return eigenvalues;
}

/*
 * Returns diagonal block matrix of MTRX_1, 
 * MTRX_2 and MTRX_3 (allocates memory).
 */
struct matrix *matrix_diagonal_compose(struct matrix *mtrx_1, 
    struct matrix *mtrx_2, struct matrix *mtrx_3)
{
    size_t order_1 = 0;
    size_t order_2 = 0;
    size_t order_3 = 0;
    if (mtrx_1)
        order_1 = mtrx_1->rows;
    if (mtrx_2)
        order_2 = mtrx_2->rows;
    if (mtrx_3)
        order_3 = mtrx_3->rows;
    
    size_t order = 0;
    if (order_1)
        order += order_1;
    if (order_2)    
        order += order_2;
    if (order_3)
        order += order_3;

    struct matrix *diag_comp = matrix_get_zeros(order);

    size_t left_coord = 0;
    if (mtrx_1) {
        matrix_paste_matrix(diag_comp, mtrx_1, left_coord, left_coord);
        left_coord += order_1;
    }
    if (mtrx_2) {
        matrix_paste_matrix(diag_comp, mtrx_2, left_coord, left_coord);
        left_coord += order_2;
    }
    if (mtrx_3)
        matrix_paste_matrix(diag_comp, mtrx_3, left_coord, left_coord);
    
    return diag_comp;
}

/*
 * Returns the determinant of the MTRX.
 */
double _Complex matrix_determinant(struct matrix *mtrx)
{
    size_t order = mtrx->rows;

    double _Complex determinant = zero;

    if (order == 1)
        return (mtrx->values)[0];

    if (order == 2) {
        determinant = (mtrx->values)[0] * (mtrx->values)[3] - 
            (mtrx->values)[1] * (mtrx->values)[2];
        // Handling value close to zero
        complex_threshold_cut(determinant);
        return determinant;
    }

    // Calculating determinant by column 
    // expansion (first column contains more zeros)
    for (size_t i = 0; i < order; i++)
        if ((mtrx->values)[i] != zero) {
            double _Complex sign;
            if (i % 2 == 0)
                sign = 1.0 + 0.0 * I;
            else
                sign = -1.0 + 0.0 * I;
            struct matrix *min_mtrx = matrix_get_minor_matrix(mtrx, i, 0);
            determinant += 
                sign * (mtrx->values)[i] * matrix_determinant(min_mtrx);
            matrix_free(min_mtrx);
        }

    // Handling value close to zero
    complex_threshold_cut(determinant);
    return determinant;
}

/*
 * Retuns MTRX inverse (allocates memory).
 */
struct matrix *matrix_inverse(struct matrix *mtrx)
{
    size_t order = mtrx->rows;

    double _Complex determinant = matrix_determinant(mtrx);

    if (determinant == zero) 
        return NULL;

    struct matrix *inverse = matrix_init(order, order);

    for (size_t i = 0; i < order; i++)
        for (size_t j = 0; j < order; j++) {
            struct matrix *min_mtrx = matrix_get_minor_matrix(mtrx, i, j);
            double _Complex sign;
            if ((i + j) % 2 == 0)
                sign = 1.0 + 0.0 * I;
            else
                sign = -1.0 + 0.0 * I;
            double _Complex minor = sign * matrix_determinant(min_mtrx);
            matrix_free(min_mtrx);
            (inverse->values)[i + j * order] = minor;
        }

    struct matrix *transposed_inverse = matrix_hermitian_transpose(inverse);
    matrix_free(inverse);

    inverse = matrix_const_mul(transposed_inverse, 1 / determinant);
    matrix_free(transposed_inverse);

    return inverse;
}

/*
 * Probably transforms MTRX_A to the upper quasitriangular form.
 */
void qz_step(struct matrix **mtrx_a, struct matrix **mtrx_b)
{
    size_t order = (*mtrx_a)->rows;

    struct matrix *b_inverse = matrix_inverse(*mtrx_b);
    struct matrix *m = matrix_mul(*mtrx_a, b_inverse);
    matrix_free(b_inverse);

    struct matrix *m_lower = 
        matrix_cut_matrix(m, order - 2, order - 1, order - 2, order - 1);

    struct matrix *eigenvalues = 
        matrix_get_second_order_eigenvalues(m_lower);
    double _Complex a = eigenvalues->values[0];
    double _Complex b = eigenvalues->values[1];

    matrix_free(eigenvalues);
    matrix_free(m_lower);

    struct matrix *id_order = matrix_get_identity_matrix(order);
    struct matrix *id_first_col = 
        matrix_cut_matrix(id_order, 0, order - 1, 0, 0);
    struct matrix *a_id = matrix_const_mul(id_order, a);
    struct matrix *b_id = matrix_const_mul(id_order, b);

    matrix_free(id_order);

    struct matrix *m_a_id = matrix_dif(m, a_id);
    matrix_free(a_id);

    struct matrix *m_b_id = matrix_dif(m, b_id);
    matrix_free(m);
    matrix_free(b_id);

    struct matrix *mul = matrix_mul(m_a_id, m_b_id);
    matrix_free(m_a_id);
    matrix_free(m_b_id);

    struct matrix *xyz_vector = matrix_mul(mul, id_first_col);
    matrix_free(mul);
    matrix_free(id_first_col);

    struct matrix *xyz_cut = matrix_cut_matrix(xyz_vector, 0, 2, 0, 0);
    matrix_free(xyz_vector);

    for (size_t i = 0; i < order - 2; i++) {
        struct matrix *vector_houshold_mtrx = 
            matrix_housholder_vector_matrix(xyz_cut);

        struct matrix *id_1;
        if (i == 0)
            id_1 = NULL;
        else
            id_1 = matrix_get_identity_matrix(i);
        
        struct matrix *id_2;
        if (i == order - 3)
            id_2 = NULL;
        else
            id_2 = matrix_get_identity_matrix(order - 3 - i);

        struct matrix *diag = 
            matrix_diagonal_compose(id_1, vector_houshold_mtrx, id_2);

        matrix_free(vector_houshold_mtrx);

        struct matrix *new_mtrx_a = matrix_mul(diag, *mtrx_a);
        struct matrix *new_mtrx_b = matrix_mul(diag, *mtrx_b);
        matrix_free(diag);
        matrix_free(*mtrx_a);
        matrix_free(*mtrx_b);
        *mtrx_a = new_mtrx_a;
        *mtrx_b = new_mtrx_b;

        struct matrix *b_vector = matrix_init(3, 1);
        (b_vector->values)[0] = (new_mtrx_b->values)[i + 2 + (i + 2) * order];
        (b_vector->values)[1] = (new_mtrx_b->values)[i + 2 + (i + 1) * order];
        (b_vector->values)[2] = (new_mtrx_b->values)[i + 2 + i * order];

        struct matrix *b_hous = matrix_housholder_vector_matrix(b_vector);
        matrix_free(b_vector);
        struct matrix *b_hous_mirror = matrix_dot_mirror(b_hous);
        matrix_free(b_hous);

        diag = matrix_diagonal_compose(id_1, b_hous_mirror, id_2);
        matrix_free(b_hous_mirror);
        new_mtrx_a = matrix_mul(*mtrx_a, diag);
        new_mtrx_b = matrix_mul(*mtrx_b, diag);
        matrix_free(diag);
        matrix_free(*mtrx_a);
        matrix_free(*mtrx_b);
        *mtrx_a = new_mtrx_a;
        *mtrx_b = new_mtrx_b;

        b_vector = matrix_init(2, 1);
        (b_vector->values)[0] = (new_mtrx_b->values)[i + 1 + (i + 1) * order];
        (b_vector->values)[1] = (new_mtrx_b->values)[i + 1 + i * order];
        
        b_hous = matrix_housholder_vector_matrix(b_vector);
        matrix_free(b_vector);
        b_hous_mirror = matrix_dot_mirror(b_hous);
        matrix_free(b_hous);

        if (id_2)
            matrix_free(id_2);
        id_2 = matrix_get_identity_matrix(order - 2 - i);
        diag = matrix_diagonal_compose(id_1, b_hous_mirror, id_2);
        matrix_free(b_hous_mirror);
        if (id_1)
            matrix_free(id_1);
        matrix_free(id_2);
    
        new_mtrx_a = matrix_mul(*mtrx_a, diag);
        new_mtrx_b = matrix_mul(*mtrx_b, diag);
        matrix_free(diag);
        matrix_free(*mtrx_a);
        matrix_free(*mtrx_b);
        *mtrx_a = new_mtrx_a;
        *mtrx_b = new_mtrx_b;

        (xyz_cut->values)[0] = ((*mtrx_a)->values)[i + 1 + i * order];
        (xyz_cut->values)[1] = ((*mtrx_a)->values)[i + 2 + i * order];
        if (i != order - 3)
            (xyz_cut->values)[2] = ((*mtrx_a)->values)[i + 3 + i * order];
    }

    struct matrix *xy = matrix_init(2, 1);
    (xy->values)[0] = (xyz_cut->values)[0];
    (xy->values)[1] = (xyz_cut->values)[1];
    matrix_free(xyz_cut);

    struct matrix *xy_hous = matrix_housholder_vector_matrix(xy);
    struct matrix *id = matrix_get_identity_matrix(order - 2);
    struct matrix *diag = matrix_diagonal_compose(id, xy_hous, NULL);
    matrix_free(xy_hous);
    struct matrix *new_mtrx_a = matrix_mul(diag, *mtrx_a);
    struct matrix *new_mtrx_b = matrix_mul(diag, *mtrx_b);
    matrix_free(diag);
    matrix_free(*mtrx_a);
    matrix_free(*mtrx_b);
    *mtrx_a = new_mtrx_a;
    *mtrx_b = new_mtrx_b;

    (xy->values)[0] = 
        ((*mtrx_b)->values)[order - 1 + (order - 1) * order];
    (xy->values)[1] = 
        ((*mtrx_b)->values)[order - 1 + (order - 2) * order];

    xy_hous = matrix_housholder_vector_matrix(xy);
    matrix_free(xy);
    
    struct matrix *xy_hous_mirror = matrix_dot_mirror(xy_hous);
    matrix_free(xy_hous);

    diag = matrix_diagonal_compose(id, xy_hous_mirror, NULL);
    matrix_free(xy_hous_mirror);

    new_mtrx_a = matrix_mul(*mtrx_a, diag);
    new_mtrx_b = matrix_mul(*mtrx_b, diag);
    matrix_free(diag);
    matrix_free(*mtrx_a);
    matrix_free(*mtrx_b);
    *mtrx_a = new_mtrx_a;
    *mtrx_b = new_mtrx_b;
}  

/*
 * Prints general eigenvalue of first order matrices MTRX_A and MTRX_B.
 */
int matrix_general_1eigenvalue(struct matrix *mtrx_a, 
    struct matrix *mtrx_b)
{
    if ((mtrx_b->values)[0] == zero) {
        if ((mtrx_a->values)[0] == zero) {
            printf("Infinte number of eigenvalues (initial matrix B "
                "was singular)\n");
                
            matrix_free(mtrx_a);
            matrix_free(mtrx_b);

            return 1;
        }
        printf("One soulution is lost (initial matrix B was singular)\n");

        matrix_free(mtrx_a);
        matrix_free(mtrx_b);

        return 0;
    }

    double _Complex lambda = (mtrx_a->values)[0] / (mtrx_b->values)[0];
    complex_print(lambda);

    matrix_free(mtrx_a);
    matrix_free(mtrx_b);

    return 0;
}

/*
 * Prints general eigenvalue of first order matrices MTRX_A and MTRX_B.
 */
int matrix_general_2eigenvalues(struct matrix *mtrx_a, 
    struct matrix *mtrx_b)
{
    double _Complex a = (mtrx_b->values)[0] * (mtrx_b->values)[3];
    double _Complex b = (mtrx_b->values)[2] * (mtrx_a->values)[1] - 
        (mtrx_b->values)[3] * (mtrx_a->values)[0] - 
        (mtrx_b->values)[0] * (mtrx_a->values)[3];
    double _Complex c = (mtrx_a->values)[0] * (mtrx_a->values)[3] - 
        (mtrx_a->values)[1] * (mtrx_a->values)[2];

    if (a == zero) {
        if (b == zero) {
            if (c == zero)
                printf("Infinite number of eigenvalues (initial matrix B "
                    "was singular)\n");
            else 
                printf("There is no eigenvalues (initial matrix B "
                    "was singular)\n");
            
            matrix_free(mtrx_a);
            matrix_free(mtrx_b);
            
            return 1;
        }
        else {
            printf("One eigenvalue is lost (initial matrix B "
                "was singular)\n");
            double _Complex lambda = - c / b;
            complex_print(lambda);
            
            matrix_free(mtrx_a);
            matrix_free(mtrx_b);
            
            return 0;
        }
    }
    
    double _Complex d = b * b - 4 * a * c;

    double _Complex lambda_1 = (- b - csqrt(d)) / (2 * a);
    double _Complex lambda_2 = (- b + csqrt(d)) / (2 * a);

    complex_print(lambda_1);
    complex_print(lambda_2);

    matrix_free(mtrx_a);
    matrix_free(mtrx_b);

    return 0;
}

/*
 * Returns index of zero on the diagonal of MTRX.
 */
size_t matrix_diagonal_zero(struct matrix *mtrx)
{
    size_t order = mtrx->rows;

    for (size_t i = 0; i < order; i++)
        if ((mtrx->values)[i + i * order] == zero)
            return i;
    
    return order;
}

/*
 * Moves zero on MTRX_B diagonal to the (n, n) position and 
 * sets zero on position (n, n - 1) of MTRX_A.
 */
void zero_shift(struct matrix **mtrx_a, struct matrix **mtrx_b, size_t zero_pos)
{
    size_t order = (*mtrx_a)->rows;

    for (size_t i = zero_pos + 1; i < order; i++) {
        struct givens_pair *givens_pair = 
            givens_pair_find(((*mtrx_b)->values)[i - 1 + i * order],
                ((*mtrx_b)->values)[i + i * order]);
        
        struct matrix *cut_a = 
            matrix_cut_matrix(*mtrx_a, i - 1, i, 0, order - 1);
        struct matrix *cut_a_rot = matrix_givens_left_mul(cut_a, givens_pair);
        matrix_free(cut_a);
        matrix_paste_matrix(*mtrx_a, cut_a_rot, i - 1, 0);
        matrix_free(cut_a_rot);
        
        struct matrix *cut_b = 
            matrix_cut_matrix(*mtrx_b, i - 1, i, 0, order - 1);
        struct matrix *cut_b_rot = matrix_givens_left_mul(cut_b, givens_pair);
        matrix_free(cut_b);
        matrix_paste_matrix(*mtrx_b, cut_b_rot, i - 1, 0);
        matrix_free(cut_b_rot);

        free(givens_pair);

        if (i > 1) {
            givens_pair = 
                givens_pair_find(((*mtrx_a)->values)[i + (i - 2) * order],
                    ((*mtrx_a)->values)[i + (i - 1) * order]);
            cut_a = matrix_cut_matrix(*mtrx_a, 0, order - 1, i - 2, i - 1);
            cut_a_rot = matrix_givens_right_mul(cut_a, givens_pair);
            matrix_free(cut_a);
            matrix_paste_matrix(*mtrx_a, cut_a_rot, 0, i - 2);
            matrix_free(cut_a_rot);
            
            cut_b = matrix_cut_matrix(*mtrx_b, 0, order - 1, i - 2, i - 1);
            cut_b_rot = matrix_givens_right_mul(cut_b, givens_pair);
            matrix_free(cut_b);
            matrix_paste_matrix(*mtrx_b, cut_b_rot, 0, i - 2);
            matrix_free(cut_b_rot);

            free(givens_pair);
        }
    }

    struct givens_pair *givens_pair = 
        givens_pair_find(((*mtrx_a)->values)[order -1 + (order - 2) * order],
            ((*mtrx_a)->values)[order - 1 + (order - 1) * order]);
    
    struct matrix *cut_a = 
        matrix_cut_matrix(*mtrx_a, 0, order - 1, order - 2, order - 1);
    struct matrix *cut_a_rot = matrix_givens_right_mul(cut_a, givens_pair);
    matrix_free(cut_a);
    matrix_paste_matrix(*mtrx_a, cut_a_rot, 0, order - 2);
    matrix_free(cut_a_rot);

    struct matrix *cut_b = 
        matrix_cut_matrix(*mtrx_b, 0, order - 1, order - 2, order - 1);
    struct matrix *cut_b_rot = matrix_givens_right_mul(cut_b, givens_pair);
    matrix_free(cut_b);
    matrix_paste_matrix(*mtrx_b, cut_b_rot, 0, order - 2);
    matrix_free(cut_b_rot);

    free(givens_pair);
}

/*
 * Reduces sizes of matrices MTRX_A and MTRX_B (if MTRX_B is singular). 
 * Moreover, checks if there is infinite number of 
 * eigenvalues (when (n,n) MTRX_A value is also zero).
 */
int diagonal_zero_reduce(struct matrix **mtrx_a, struct matrix **mtrx_b)
{
    size_t order = (*mtrx_a)->rows;
    size_t zero_position = matrix_diagonal_zero(*mtrx_b);

    while (zero_position != order) {
        
        if (order == 1) {
            if (((*mtrx_a)->values)[0] == zero) {
                printf("Infinite number of eigenvalues (initial matrix B "
                    "was singular)\n");
                return 1;
            }
            printf("One eigenvalue is lost (initial matrix B "
                "was singular)\n");
            return 1;
        }
        
        zero_shift(mtrx_a, mtrx_b, zero_position);

        if (((*mtrx_a)->values)[order - 1 + (order - 1) * order] == zero) {
            printf("Infinite number of eigenvalues (initial matrix B "
                    "was singular)\n");
            return 1;
        }
        else
            printf("One eigenvalue is lost (initial matrix B "
                "was singular)\n");
        
        struct matrix *new_mtrx_a = 
            matrix_cut_matrix(*mtrx_a, 0, order - 2, 0, order - 2);
        matrix_free(*mtrx_a);
        *mtrx_a = new_mtrx_a;
        
        struct matrix *new_mtrx_b = 
            matrix_cut_matrix(*mtrx_b, 0, order - 2, 0, order - 2);
        matrix_free(*mtrx_b);
        *mtrx_b = new_mtrx_b;

        order = (*mtrx_a)->rows;
        zero_position = matrix_diagonal_zero(*mtrx_b);
    }

    return 0;
}

/*
 * Prints eigenvalues of matrix pair MTRX_A and MTRX_B. 
 * Solves generalised problem det(A - lambda * B) = 0.
 */
void qz_process(struct matrix **mtrx_a, struct matrix **mtrx_b)
{   
    size_t order = (*mtrx_a)->rows;

    if (order == 2) {
        matrix_general_2eigenvalues(*mtrx_a, *mtrx_b);
        return;
    }
    if (order == 1) {
        matrix_general_1eigenvalue(*mtrx_a, *mtrx_b);
        return;
    }
        
    matrix_hessenberg_triangular_transformation(mtrx_a, *mtrx_b);

    if (diagonal_zero_reduce(mtrx_a, mtrx_b))
        return;

    while (1) {
        order = (*mtrx_a)->rows;
        for (size_t j = 1; j < order; j++)
            if (cabs(((*mtrx_a)->values)[j + (j - 1) * order]) <= 
            epsilon * (cabs(((*mtrx_a)->values)[j - 1 + (j - 1) * order]) + 
                cabs(((*mtrx_a)->values)[j + j * order])))
                    ((*mtrx_a)->values)[j + (j - 1) * order] = zero;
        
        int found = 1;
        while (found) {
            order = (*mtrx_a)->rows;
        
            if (order == 2) {
                matrix_general_2eigenvalues(*mtrx_a, *mtrx_b);
                return;
            }

            if (order == 1) {
                matrix_general_1eigenvalue(*mtrx_a, *mtrx_b);
                return;
            }

            if (((*mtrx_a)->values)[2 + order] == zero) {
                struct matrix *small_a = matrix_cut_matrix(*mtrx_a, 0, 1, 0, 1);
                struct matrix *new_mtrx_a = 
                    matrix_cut_matrix(*mtrx_a, 2, order - 1, 2, order - 1);
                struct matrix *small_b = matrix_cut_matrix(*mtrx_b, 0, 1, 0, 1);
                struct matrix *new_mtrx_b = 
                    matrix_cut_matrix(*mtrx_b, 2, order - 1, 2, order - 1);
                matrix_general_2eigenvalues(small_a, small_b);
                matrix_free(*mtrx_a);
                matrix_free(*mtrx_b);
                *mtrx_a = new_mtrx_a;
                *mtrx_b = new_mtrx_b;
                continue;
            }

            if (((*mtrx_a)->values)[1] == zero) {
                struct matrix *small_a = matrix_cut_matrix(*mtrx_a, 0, 0, 0, 0);
                struct matrix *new_mtrx_a = 
                    matrix_cut_matrix(*mtrx_a, 1, order - 1, 1, order - 1);
                struct matrix *small_b = matrix_cut_matrix(*mtrx_b, 0, 0, 0, 0);
                struct matrix *new_mtrx_b = 
                    matrix_cut_matrix(*mtrx_b, 1, order - 1, 1, order - 1);
                matrix_general_1eigenvalue(small_a, small_b);
                matrix_free(*mtrx_a);
                matrix_free(*mtrx_b);
                *mtrx_a = new_mtrx_a;
                *mtrx_b = new_mtrx_b;
                continue;
            }

            if (((*mtrx_a)->values)[order - 2 + (order - 3) * order] == zero) {
                struct matrix *small_a = 
                    matrix_cut_matrix(*mtrx_a, order - 2, order - 1,
                        order - 2, order - 1);
                struct matrix *new_mtrx_a = 
                    matrix_cut_matrix(*mtrx_a, 0, order - 3, 0, order - 3);
                struct matrix *small_b = 
                    matrix_cut_matrix(*mtrx_b, order - 2, order - 1, 
                        order - 2, order - 1);
                struct matrix *new_mtrx_b = 
                    matrix_cut_matrix(*mtrx_b, 0, order - 3, 0, order - 3);
                matrix_general_2eigenvalues(small_a, small_b);
                matrix_free(*mtrx_a);
                matrix_free(*mtrx_b);
                *mtrx_a = new_mtrx_a;
                *mtrx_b = new_mtrx_b;
                continue;
            }

            if (((*mtrx_a)->values)[order - 1 + (order - 2) * order] == zero) {
                struct matrix *small_a = 
                    matrix_cut_matrix(*mtrx_a, order - 1, order - 1, 
                        order - 1, order - 1);
                struct matrix *new_mtrx_a = 
                    matrix_cut_matrix(*mtrx_a, 0, order - 2, 0, order - 2);
                struct matrix *small_b = 
                    matrix_cut_matrix(*mtrx_b, order - 1, order - 1, 
                        order - 1, order - 1);
                struct matrix *new_mtrx_b = 
                    matrix_cut_matrix(*mtrx_b, 0, order - 2, 0, order - 2);
                matrix_general_1eigenvalue(small_a, small_b);
                matrix_free(*mtrx_a);
                matrix_free(*mtrx_b);
                *mtrx_a = new_mtrx_a;
                *mtrx_b = new_mtrx_b;
                continue;
            }

            found = 0;
        }

        qz_step(mtrx_a, mtrx_b);
    }
}