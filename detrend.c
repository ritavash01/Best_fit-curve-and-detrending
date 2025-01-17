#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define M 10 //number of data points 
#define N 3 //degree of polynomial 


double legendre_polynomial(int k, double x){
    if ( k == 0) return 1.0;
    if ( k == 1) return x;

    double Pk_minus_2 = 1.0;
    double Pk_minus_1 = x; 
    double Pk = 0.0;

    for (int n = 2; n <=k; n++){
        Pk = ((2 * n - 1) * x * Pk_minus_1 - (n - 1) * Pk_minus_2) / n;
        Pk_minus_2 = Pk_minus_1;
        Pk_minus_1 = Pk;

    return Pk;
    }
}
// Cholesky Decomposition 
int cholesky_decomposition(double A[N][N], double L[N][N]){
    double sum = 0.0;
    for (int i = 0; i < N; i++){
        for (int j = 0; j <= i; j++){
            for (int k = 0; k < j; k++){
                sum += L[i][k]*L[j][k];
            }
            if (i == j){
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
double evalute_polynomial( double x, double *c){
    double result = 0; 
    for (int i = 0; i < N; i++) {
        result += c[i] * legendre_polynomial(i, x);
    }
    return result;
}

int main() {
    double x[M] = {-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8}; // If we are fitting the frequency values then it would be the time axis points 
    double y[M] = {1, 0.64, 0.36, 0.16, 0.04, 0, -0.04, -0.16, -0.36, -0.64}; // Observed intensity values

    //Initialise the different matrices 
    double A[N][N] = {0};  // Normal matrix
    double b[N] = {0};     // Right-hand side vector
    double L[N][N] = {0};  // Cholesky lower triangular matrix
    double c[N] = {0};     // Coefficients


    // Compute A and b just using the formula 
    for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {
            for (int i = 0; i < M; i++) {
                double Pj = legendre_polynomial(j, x[i]);
                double Pk = legendre_polynomial(k, x[i]);
                A[j][k] += Pj * Pk;
            }
        }
        for (int i = 0; i < M; i++) {
            b[j] += y[i] * legendre_polynomial(j, x[i]);
        }
    }



    cholesky_decomposition(A, L);

    // Solve for coefficients
    double z[N];
    forward_substitution(L, b, z);
    backward_substitution(L, z, c);

   

    //Find the trend and subtract to get the residuals
    double residuals[M];
    for (int i = 0; i < M; i ++){
        double trend = evalute_polynomial(x[i], c);
        residuals[i] = y[i] - trend;
    }

      printf("Original Data, Trend, Residuals:\n");
    for (int i = 0; i < M; i++) {
        double trend = evalute_polynomial(x[i], c);
        printf("x[%d] = %f, y[%d] = %f, Trend = %f, Residual = %f\n", i, x[i], i, y[i], trend, residuals[i]);
    }

 

    return 0;


}
