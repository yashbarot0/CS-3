#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * This code sets up the linear system Ay = b for the Poisson equation
 * -∆u = f on the unit square [0,1]² with u = 0 on the boundary.
 */

// Function f(x,y) = 2π² sin(πx) sin(πy)
double f(double x, double y) {
    return 2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
}

// Create the matrix A and right-hand side b for the Poisson problem
void setup_poisson_system(int N, double **A, double *b) {
    int n = (N-1) * (N-1); // Number of interior grid points
    double h = 1.0 / N;    // Grid spacing
    double h2 = h * h;     // h²

    // Initialize A to zeros
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = 0.0;
        }
    }

    // Fill matrix A and vector b
    for (int i = 0; i < N-1; i++) {
        for (int j = 0; j < N-1; j++) {
            int idx = i + j * (N-1); // Linear index of grid point (i+1, j+1)
            
            // Set diagonal element (coefficient of u_{i,j})
            A[idx][idx] = 4.0;
            
            // Set off-diagonal elements (for neighboring points)
            if (i > 0) {
                A[idx][idx-1] = -1.0; // left neighbor
            }
            if (i < N-2) {
                A[idx][idx+1] = -1.0; // right neighbor
            }
            if (j > 0) {
                A[idx][idx-(N-1)] = -1.0; // bottom neighbor
            }
            if (j < N-2) {
                A[idx][idx+(N-1)] = -1.0; // top neighbor
            }
            
            // Set right-hand side b
            double x = (i+1) * h;
            double y = (j+1) * h;
            b[idx] = -f(x, y) * h2;
        }
    }
}

// Function to display a matrix (for debugging)
void display_matrix(double **A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%6.2f ", A[i][j]);
        }
        printf("\n");
    }
}

// Function to display a vector (for debugging)
void display_vector(double *v, int n) {
    for (int i = 0; i < n; i++) {
        printf("%6.2f\n", v[i]);
    }
}

// Simple test to validate the system setup
void test_setup(int N) {
    int n = (N-1) * (N-1);
    
    // Allocate memory for A and b
    double **A = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        A[i] = (double *)malloc(n * sizeof(double));
    }
    double *b = (double *)malloc(n * sizeof(double));
    
    // Setup the system
    setup_poisson_system(N, A, b);
    
    // Display for a small test case
    if (N <= 5) {
        printf("Matrix A:\n");
        display_matrix(A, n);
        printf("\nVector b:\n");
        display_vector(b, n);
    } else {
        printf("System set up for N = %d (matrix size = %d x %d)\n", N, n, n);
    }
    
    // Free memory
    for (int i = 0; i < n; i++) {
        free(A[i]);
    }
    free(A);
    free(b);
}

int main() {
    // Test with a small grid
    test_setup(4); // N=4 gives a 9x9 matrix (3x3 interior grid points)
    
    return 0;
}