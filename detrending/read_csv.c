#include <stdio.h>
#include <stdlib.h>

int main() {
    // Define the expected matrix shape
    int frequency_points = 4096;  // Number of rows
    int time_points = 300;        // Number of columns

    // Calculate the expected file size (4096 * 300 * sizeof(double) bytes)
    size_t expected_size = frequency_points * time_points * sizeof(double);

    // Open the binary file for reading
    FILE *fstream = fopen("/home/ritavash/Documents/Simulated_Data/complex_gaussian_noise.dat", "rb");
    if (fstream == NULL) {
        printf("\nFile opening failed\n");
        return -1;
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
        return -1;
    }

    // Allocate memory for the matrix as a contiguous block
    double *mat = (double *)malloc(frequency_points * time_points * sizeof(double));
    if (mat == NULL) {
        printf("Error: Memory allocation failed\n");
        fclose(fstream);
        return -1;
    }

    // Read the data into the matrix
    size_t read_size = fread(mat, sizeof(double), frequency_points * time_points, fstream);
    if (read_size != frequency_points * time_points) {
        printf("Error: Not enough data read from the file\n");
        free(mat);
        fclose(fstream);
        return -1;
    }

    // Print the shape of the matrix
    printf("Matrix shape: %d x %d\n", frequency_points, time_points);

    // Free the allocated memory
    free(mat);

    fclose(fstream);
    return 0;
}

