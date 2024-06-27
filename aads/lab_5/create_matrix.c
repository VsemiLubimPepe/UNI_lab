#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


double get_random_double(double min, double max)
{
    double range = max - min;
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        printf("Usage: ./create_matrix ORDER\n");
        return -1;
    }

    size_t order = strtoull(argv[1], NULL, 10);
    if ((errno == EINVAL) || (errno == ERANGE)) {
        printf("Argument is not a number or it doesn't fit in interval.\n");
        return -1;
    }

    FILE *matrix_a_file = fopen("matrices/matrix_a", "w");
    if (!matrix_a_file) {
        printf("Couldn't open file to write matrix values.\n");
        return -1;
    }

    fprintf(matrix_a_file, "%zu %zu \n", order, order);

    srand(time(NULL));
    
    for (size_t i = 0; i < order; i++) {
        for (size_t j = 0; j < order; j++) {
            fprintf(matrix_a_file, "%f%+f*i ",
                get_random_double(-1000.0, 1000.0), 
                get_random_double(-1000.0, 1000.0));
        }
        fprintf(matrix_a_file, "\n");
    }

    fclose(matrix_a_file);
    
    FILE *matrix_b_file = fopen("matrices/matrix_b", "w");
    if (!matrix_b_file) {
        printf("Couldn't open file to write matrix values.\n");
        return -1;
    }

    fprintf(matrix_b_file, "%zu %zu \n", order, order);
    
    for (size_t i = 0; i < order; i++) {
        for (size_t j = 0; j < order; j++) {
            fprintf(matrix_b_file, "%f%+f*i ",
                get_random_double(-1000.0, 1000.0), 
                get_random_double(-1000.0, 1000.0));
        }
        fprintf(matrix_b_file, "\n");
    }

    fclose(matrix_b_file);

    return 0;
}
