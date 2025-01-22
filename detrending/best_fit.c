#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define N 3  // Degree of polynomial 

double legendre_polynomial(int k, double x) {
    if (k == 0) return 1.0;
    if (k == 1) return x;

    double Pk_minus_2 = 1.0;
    double Pk_minus_1 = x;
    double Pk = 0.0;

    for (int n = 2; n <= k; n++) {
        Pk = ((2 * n - 1) * x * Pk_minus_1 - (n - 1) * Pk_minus_2) / n;
        Pk_minus_2 = Pk_minus_1;
        Pk_minus_1 = Pk;
    }
    return Pk;
}

// Cholesky Decomposition
int cholesky_decomposition(double A[N][N], double L[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += L[i][k] * L[j][k];
            }
            if (i == j) {
                L[i][j] = sqrt(A[i][j] - sum);
            } else {
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
        for (int j = i + 1; j < N; j++) {
            L[i][j] = 0.0; // Fill upper triangle with zeros
        }
    }
    return 1;
}

// Forward substitution to solve Lz = b
void forward_substitution(double L[N][N], double b[N], double z[N]) {
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * z[j];
        }
        z[i] = (b[i] - sum) / L[i][i];
    }
}

// Backward substitution to solve L^T c = z
void backward_substitution(double L[N][N], double z[N], double c[N]) {
    for (int i = N - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < N; j++) {
            sum += L[j][i] * c[j];
        }
        c[i] = (z[i] - sum) / L[i][i];
    }
}

double evalute_polynomial(double x, double *c) {
    double result = 0;
    for (int i = 0; i < N; i++) {
        result += c[i] * legendre_polynomial(i, x);
    }
    return result;
}

void read_dat_to_matrix(const char *filename, double **matrix, int *rows, int *cols) {
    // Define the expected matrix shape
    int frequency_points = 4096;  // Number of rows
    int time_points = 300;        // Number of columns

    // Calculate the expected file size (4096 * 300 * sizeof(double) bytes)
    size_t expected_size = frequency_points * time_points * sizeof(double);

    // Open the binary file for reading
    FILE *fstream = fopen(filename, "rb");
    if (fstream == NULL) {
        printf("\nFile opening failed\n");
        exit(EXIT_FAILURE);
    }

    // Get the file size
    fseek(fstream, 0, SEEK_END);
    size_t file_size = ftell(fstream);
    fseek(fstream, 0, SEEK_SET);

    // Check if the file size matches the expected size
    if (file_size != expected_size) {
        printf("Error: The file size does not match the expected matrix size.\n");
        printf("File size: %zu bytes\n", file_size);
        printf("Expected matrix size: %zu bytes\n", expected_size);
        fclose(fstream);
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the matrix as a 1D array
    *matrix = (double *)malloc(frequency_points * time_points * sizeof(double));
    if (*matrix == NULL) {
        printf("Memory allocation failed for matrix\n");
        fclose(fstream);
        exit(EXIT_FAILURE);
    }

    // Read the data into the matrix (1D array)
    size_t read_size = fread(*matrix, sizeof(double), frequency_points * time_points, fstream);
    if (read_size != frequency_points * time_points) {
        printf("Error: Not enough data read from the file\n");
        fclose(fstream);
        free(*matrix);
        exit(EXIT_FAILURE);
    }

    *rows = frequency_points;
    *cols = time_points;
    printf("File read successfully. Matrix dimensions: %d x %d\n", *rows, *cols);
    fclose(fstream);
}
void write_matrix_to_csv(const char *filename, double *matrix, int rows, int cols) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Use appropriate formatting for CSV
            fprintf(file, "%.6f", matrix[i * cols + j]);

            if (j < cols - 1) {
                fprintf(file, ",");
            }
        }
        fprintf(file, "\n");
    }
    printf("Detrended data saved. Rows: %d, Columns: %d\n", rows, cols);

    fclose(file);
}

double *allocate_matrix(int rows, int cols) {
    printf("Trying to allocate %zu bytes for detrended matrix\n", rows * cols * sizeof(double));


    double *matrix = (double *)malloc(rows * cols * sizeof(double));
    if (matrix == NULL) {
        // Print a detailed error message and exit if memory allocation fails
        perror("Error allocating memory for detrended matrix");
        exit(EXIT_FAILURE);
    }
    printf("Memory allocated for detrended matrix\n");  // Print success message
    return matrix;
}

void free_matrix(double *matrix) {
    free(matrix);
}

int main() {
    const char *input_filename = "/home/ritavash/Documents/Simulated_Data/complex_gaussian_noise.dat";
    const char *output_filename = "detrended_data.csv";
    double *matrix = NULL;
    int rows = 0, cols = 0;

    // Read .dat file
    read_dat_to_matrix(input_filename, &matrix, &rows, &cols);
    printf("Matrix dimensions: %d x %d\n", rows, cols);

    // Allocate detrended matrix
    printf("Trying to allocate %zu bytes for detrended matrix\n", rows * cols * sizeof(double));

    double *detrended_matrix = allocate_matrix(rows, cols);
    printf("Memory allocated for detrended matrix");
    // Process data
    double A[N][N] = {0}, L[N][N] = {0}, b[N] = {0}, c[N] = {0};
    for (int freq = 0; freq < cols; freq++) {
        printf("Processing frequency %d\n", freq);
        memset(A, 0, sizeof(A));
        memset(b, 0, sizeof(b));
        memset(c, 0, sizeof(c));

        // Fill A and b
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                for (int i = 0; i < rows; i++) {
                    double Pj = legendre_polynomial(j, matrix[i * cols]);
                    double Pk = legendre_polynomial(k, matrix[i * cols]);
                    A[j][k] += Pj * Pk;
                }
            }
            for (int i = 0; i < rows; i++) {
                b[j] += matrix[i * cols + freq] * legendre_polynomial(j, matrix[i * cols]);
            }
        }

        cholesky_decomposition(A, L);

        double z[N];
        forward_substitution(L, b, z);
        backward_substitution(L, z, c);

        for (int i = 0; i < rows; i++) {
            double trend = evalute_polynomial(matrix[i * cols], c);
            detrended_matrix[i * cols + freq] = matrix[i * cols + freq] - trend;
        }
    }

    // Write detrended data to .dat file
    write_matrix_to_csv(output_filename, detrended_matrix, rows, cols);

   


    // Free memory
    free_matrix(matrix);
    free_matrix(detrended_matrix);

    return 0;
}
