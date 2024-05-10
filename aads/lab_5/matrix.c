#include <complex.h>
#include <malloc.h>

/*
 * Numbers stored column-wise to easily get
 * vectors for Hausholder method
 */
struct matrix {
    size_t rows;
    size_t columns;
    double _Complex *values;
};

struct matrix *matrix_init(size_t rows, size_t columns)
{
    struct matrix *mtrx = (struct matrix *)malloc(sizeof(struct matrix));
    if (!mtrx)
        return NULL;

    mtrx->values =
        (double _Complex *)malloc(sizeof(double _Complex) * columns * rows);
    if (!mtrx->values) {
        free(mtrx);
        return NULL;
    }

    mtrx->rows = rows;
    mtrx->columns = columns;

    return mtrx;
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
            for (size_t n = 0; n < common_dimension; n++) {
                (mtrx->values)[i + j * columns] +=
                    (mtrx_1->values)[i + n] * (mtrx_2->values)[n + j];
            }
        }
    }

    return mtrx;
}
