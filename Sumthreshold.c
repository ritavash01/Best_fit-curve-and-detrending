#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SMOOTHING_WINDOW 2  // Half-window size for smoothing
#define LAMBDA_UP 1.5       // Weight for upward deviation
#define LAMBDA_DOWN 1.0     // Weight for downward deviation

// Function for asymmetric smoothing
void asymmetric_smoothing(double *signal, double *baseline, int size, int window, double lambda_up, double lambda_down) {
    double *temp_baseline = (double *)malloc(size * sizeof(double));
    if (temp_baseline == NULL) {
        perror("Memory allocation failed for temporary baseline");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < size; i++) {
        temp_baseline[i] = signal[i];
    }

    for (int iter = 0; iter < 10; iter++) {
        for (int i = 0; i < size; i++) {
            double sum_weights = 0.0;
            double weighted_sum = 0.0;

            for (int j = -window; j <= window; j++) {
                int idx = i + j;
                if (idx >= 0 && idx < size) {
                    double weight = 1.0;
                    if (signal[idx] > temp_baseline[idx]) {
                        weight = lambda_up;
                    } else if (signal[idx] < temp_baseline[idx]) {
                        weight = lambda_down;
                    }
                    sum_weights += weight;
                    weighted_sum += weight * temp_baseline[idx];
                }
            }
            baseline[i] = weighted_sum / sum_weights;
        }

        // Update the temporary baseline
        for (int i = 0; i < size; i++) {
            temp_baseline[i] = baseline[i];
        }
    }

    free(temp_baseline);
}

// Function to read binary data into matrix
void read_dat_to_matrix(const char *filename, double **matrix, int *rows, int *cols) {
    int frequency_points = 4096;
    int time_points = 300;
    size_t expected_size = frequency_points * time_points * sizeof(double);

    FILE *fstream = fopen(filename, "rb");
    if (fstream == NULL) {
        printf("\nFile opening failed\n");
        exit(EXIT_FAILURE);
    }

    fseek(fstream, 0, SEEK_END);
    size_t file_size = ftell(fstream);
    fseek(fstream, 0, SEEK_SET);

    if (file_size != expected_size) {
        printf("Error: The file size does not match the expected matrix size.\n");
        fclose(fstream);
        exit(EXIT_FAILURE);
    }

    *matrix = (double *)malloc(frequency_points * time_points * sizeof(double));
    if (*matrix == NULL) {
        printf("Memory allocation failed for matrix\n");
        fclose(fstream);
        exit(EXIT_FAILURE);
    }

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

// Function to write matrix to CSV
void write_matrix_to_csv(const char *filename, double *matrix, int rows, int cols) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            fprintf(file, "%.6f", matrix[i * cols + j]);
            if (j < cols - 1) {
                fprintf(file, ",");
            }
        }
        fprintf(file, "\n");
    }
    printf("Processed data saved to %s. Rows: %d, Columns: %d\n", filename, rows, cols);
    fclose(file);
}

// Function to write marker matrix to CSV
void write_marker_matrix_to_csv(const char *filename, int *matrix, int rows, int cols) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            fprintf(file, "%d", matrix[i * cols + j]);
            if (j < cols - 1) {
                fprintf(file, ",");
            }
        }
        fprintf(file, "\n");
    }
    printf("Marker matrix saved to %s. Rows: %d, Columns: %d\n", filename, rows, cols);
    fclose(file);
}

// Function to allocate matrix memory
double *allocate_matrix(int rows, int cols) {
    double *matrix = (double *)malloc(rows * cols * sizeof(double));
    if (matrix == NULL) {
        perror("Error allocating memory for matrix");
        exit(EXIT_FAILURE);
    }
    return matrix;
}

// Function to free matrix memory
void free_matrix(double *matrix) {
    free(matrix);
}

// Function to calculate median
int compare(const void *a, const void *b) {
    double diff = (*(double *)a - *(double *)b);
    return (diff > 0) - (diff < 0); // Return 1, -1, or 0
}

// Function to calculate median
double calculate_median(double *array, int size) {
    qsort(array, size, sizeof(double), compare);
    if (size % 2 == 0) {
        return (array[size / 2 - 1] + array[size / 2]) / 2.0;
    } else {
        return array[size / 2];
    }
}

// Function to calculate MAD
double calculate_MAD(double *data, int size) {
    double median = calculate_median(data, size);
    double *abs_deviation = (double *)malloc(size * sizeof(double));
    if (!abs_deviation) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < size; i++) {
        abs_deviation[i] = fabs(data[i] - median);
    }

    double MAD = calculate_median(abs_deviation, size) * 1.4826;
    free(abs_deviation);
    return MAD;
}

// Function to apply MAD filter
double apply_MAD_filter(double *data, int size, double beta) {
    double median = calculate_median(data, size);
    double mad = calculate_MAD(data, size);
    double threshold_1 = median + beta * mad;
    return threshold_1;
}

// Function to calculate threshold list and write to CSV
double *threshold_list(double *data, int size, double beta, double rho, int bin_array_size, const char *output_filename) {
    double threshold_1 = apply_MAD_filter(data, size, beta);
    double *thresholds = (double *)malloc(bin_array_size * sizeof(double));
    if (thresholds == NULL) {
        printf("Memory allocation failed!\n");
    }

    // Open the CSV file for writing the thresholds
    FILE *file = fopen(output_filename, "w");
    if (file == NULL) {
        perror("Error opening file for writing thresholds");
        exit(EXIT_FAILURE);
    }

    // Write header for the CSV
    fprintf(file, "Index,Threshold\n");

    for (int i = 0; i < bin_array_size; i++) {
        thresholds[i] = threshold_1 / pow(rho, log2(i + 1));
        fprintf(file, "%d,%.6f\n", i, thresholds[i]);
        printf("Threshold[%d] = %.6f\n", i, thresholds[i]);
    }

    // Close the file after writing
    fclose(file);

    return thresholds;
}

// Function for sumthreshold algorithm
void sumthreshold(double *data, int size, double *bin_array, double *thresholds, int *markers) {
    int *current_markers = (int *)malloc(size * sizeof(int));
    if (!current_markers) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < size; i++) {
        markers[i] = 0;
    }

    int max_iterations = 0;
    while (bin_array[max_iterations] > 0) {
        max_iterations++;
    }

    for (int iteration = 0; iteration < max_iterations; iteration++) {
        int window_size = (int)bin_array[iteration];
        double threshold = thresholds[iteration];

        for (int i = 0; i < size; i++) {
            current_markers[i] = markers[i];
        }

        for (int i = 0; i < size; i++) {
            double sum = 0.0;
            int count = 0;
            // Accumulate values within the window centered around i
            for (int j = i - window_size; j <= i + window_size; j++) {
                if (j >= 0 && j < size) {
                    sum += data[j];
                    count++;
                }
            }

            double avg = (count > 0) ? (sum / count) : 0;

            if (avg > threshold) {
                markers[i] = 1; // Mark as flagged
                data[i] = threshold; // Replace flagged value with threshold of iteration
            }
        }
    }

    free(current_markers);
}

// Apply sumthreshold in time direction
void apply_sumthreshold_in_time_direction(double *matrix, int rows, int cols, double *bin_array, double beta, double rho, int bin_array_size, int *marker_matrix) {
    double *baseline = (double *)malloc(cols * sizeof(double));
    int *markers = (int *)malloc(cols * sizeof(int));

    for (int row = 0; row < rows; row++) {
        double *signal = &matrix[row * cols];

        // Perform asymmetric smoothing to calculate the baseline
        asymmetric_smoothing(signal, baseline, cols, SMOOTHING_WINDOW, LAMBDA_UP, LAMBDA_DOWN);

        // Apply sumthreshold for filtering
        double *thresholds = threshold_list(signal, cols, beta, rho, bin_array_size, "time_thresholds.csv");
        sumthreshold(signal, cols, bin_array, thresholds, markers);

        // Update the marker matrix and matrix itself
        for (int i = 0; i < cols; i++) {
            marker_matrix[row * cols + i] = markers[i];
            matrix[row * cols + i] = signal[i];  // Update the matrix with the filtered signal
        }

        // Free memory for thresholds
        free(thresholds);
    }

    free(baseline);
    free(markers);
}

// Apply sumthreshold in frequency direction
void apply_sumthreshold_in_frequency_direction(double *matrix, int rows, int cols, double *bin_array, double beta, double rho, int bin_array_size, int *marker_matrix) {
    double *baseline = (double *)malloc(rows * sizeof(double));
    int *markers = (int *)malloc(rows * sizeof(int));

    for (int col = 0; col < cols; col++) {
        double *signal = (double *)malloc(rows * sizeof(double));

        // Initialize the signal array and markers with the previous marker array
        for (int row = 0; row < rows; row++) {
            signal[row] = matrix[row * cols + col];
            markers[row] = marker_matrix[row * cols + col];
        }

        // Apply sumthreshold for filtering
        double *thresholds = threshold_list(signal, rows, beta, rho, bin_array_size, "freq_thresholds.csv");
        sumthreshold(signal, rows, bin_array, thresholds, markers);

        // Write the filtered column back to the matrix and update the marker matrix
        for (int row = 0; row < rows; row++) {
            matrix[row * cols + col] = signal[row];
            marker_matrix[row * cols + col] = markers[row]; // Update marker matrix
        }

        // Free memory for thresholds
        free(thresholds);
        free(signal);
    }

    free(baseline);
    free(markers);
}

int main() {
    const char *input_filename = "/home/ritavash/Documents/Simulated_Data/complex_gaussian_noise.dat";
    const char *output_filename_time = "final_data_time.csv";
    const char *output_filename_freq = "final_data_freq.csv";
    const char *marker_filename_time = "marker_matrix_time.csv";
    const char *marker_filename_freq = "marker_matrix_freq.csv";

    double *matrix = NULL;
    int rows = 0, cols = 0;

    // Read the binary .dat file into a matrix
    read_dat_to_matrix(input_filename, &matrix, &rows, &cols);

    // Allocate memory for the final matrix and marker matrix
    double *final_matrix_time = allocate_matrix(rows, cols);
    memcpy(final_matrix_time, matrix, rows * cols * sizeof(double));

    int *marker_matrix = (int *)malloc(rows * cols * sizeof(int));
    memset(marker_matrix, 0, rows * cols * sizeof(int)); // Initialize marker matrix to 0

    // Define parameters for sumthreshold
    double bin_array[] = {1}; // Example bin sizes
    double beta = 5.0;              // Threshold scaling factor
    double rho = 1.5;               // Exponential decay factor
    int bin_array_size = sizeof(bin_array) / sizeof(bin_array[0]);

    // Apply sumthreshold in time direction and save the result
    apply_sumthreshold_in_time_direction(final_matrix_time, rows, cols, bin_array, beta, rho, bin_array_size, marker_matrix);
    write_matrix_to_csv(output_filename_time, final_matrix_time, rows, cols);

    // Save the marker matrix after time run
    write_marker_matrix_to_csv(marker_filename_time, marker_matrix, rows, cols);

    // Apply sumthreshold in frequency direction and save the result
    apply_sumthreshold_in_frequency_direction(final_matrix_time, rows, cols, bin_array, beta, rho, bin_array_size, marker_matrix);
    write_matrix_to_csv(output_filename_freq, final_matrix_time, rows, cols);

    // Save the marker matrix after frequency run
    write_marker_matrix_to_csv(marker_filename_freq, marker_matrix, rows, cols);

    // Free memory for allocated matrices
    free_matrix(matrix);
    free_matrix(final_matrix_time);
    free(marker_matrix);

    return 0;
}
