#include <stdio.h>
#include <stdlib.h>

int main() {
    // Define test matrix dimensions
    int test_rows = 4096*2;
    int test_cols = 300*2;
    
    // Allocate memory for the test matrix
    double *test_matrix = (double *)malloc(test_rows * test_cols * sizeof(double));

    // Check if memory allocation was successful
    if (test_matrix == NULL) {
        perror("Error allocating memory for test matrix");
        exit(EXIT_FAILURE);
    } else {
        printf("Successfully allocated memory for test matrix\n");
    }

    // Optionally: Fill matrix with some test data to verify access
    for (int i = 0; i < test_rows; i++) {
        for (int j = 0; j < test_cols; j++) {
            test_matrix[i * test_cols + j] = (double)(i * test_cols + j);
        }
    }

    // Optionally: Print part of the matrix to confirm allocation and access
    printf("Matrix[0][0] = %.2f\n", test_matrix[0]);
    printf("Matrix[99][99] = %.2f\n", test_matrix[test_rows * test_cols - 1]);

    // Free the allocated memory
    free(test_matrix);
    printf("Memory for test matrix freed successfully\n");

    return 0;
}
