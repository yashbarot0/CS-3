#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>

// Function to calculate the norm of a vector
double vector_norm(double *v, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}

// Function to perform matrix-vector multiplication for dense matrix
void matrix_vector_mult(double **A, double *v, double *result, int n) {
    for (int i = 0; i < n; i++) {
        result[i] = 0.0;
        for (int j = 0; j < n; j++) {
            result[i] += A[i][j] * v[j];
        }
    }
}

// Conjugate Gradient method that saves residuals for convergence analysis
int conjugate_gradient_with_residuals(double **A, double *b, double *x, int n, 
                                      double reltol, double abstol,
                                      int max_iter, double *residuals) {
    double *r = (double *)malloc(n * sizeof(double));
    double *p = (double *)malloc(n * sizeof(double));
    double *Ap = (double *)malloc(n * sizeof(double));
    
    // Initialize x to zeros
    for (int i = 0; i < n; i++) {
        x[i] = 0.0;
    }
    
    // r = b - A*x (since x=0, r = b)
    for (int i = 0; i < n; i++) {
        r[i] = b[i];
        p[i] = r[i];
    }
    
    double r_norm = vector_norm(r, n);
    double initial_r_norm = r_norm;
    double tol = fmax(reltol * initial_r_norm, abstol);
    
    residuals[0] = r_norm; // Save initial residual
    
    int iter;
    for (iter = 0; iter < max_iter && r_norm > tol; iter++) {
        // Compute A*p
        matrix_vector_mult(A, p, Ap, n);
        
        // Calculate alpha = (r^T * r) / (p^T * A * p)
        double r_dot_r = 0.0;
        double p_dot_Ap = 0.0;
        for (int i = 0; i < n; i++) {
            r_dot_r += r[i] * r[i];
            p_dot_Ap += p[i] * Ap[i];
        }
        double alpha = r_dot_r / p_dot_Ap;
        
        // Update x and r
        double r_dot_r_new = 0.0;
        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        
        // Calculate new r_dot_r for beta
        for (int i = 0; i < n; i++) {
            r_dot_r_new += r[i] * r[i];
        }
        
        // Calculate beta = (r_{k+1}^T * r_{k+1}) / (r_k^T * r_k)
        double beta = r_dot_r_new / r_dot_r;
        
        // Update p
        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta * p[i];
        }
        
        // Update r_norm for next iteration check
        r_norm = sqrt(r_dot_r_new);
        
        // Save residual for this iteration
        residuals[iter+1] = r_norm;
    }
    
    // Free memory
    free(r);
    free(p);
    free(Ap);
    
    return iter;
}

// Function to generate the matrix A as defined in the problem
void generate_matrix_A(double **A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = (double)(n - abs(i - j)) / n;
        }
    }
}

// Function to generate the right-hand side vector b with all 1s
void generate_vector_b(double *b, int n) {
    for (int i = 0; i < n; i++) {
        b[i] = 1.0;
    }
}

// Function to calculate the convergence rate from residuals and compare to theory
void analyze_convergence(double *residuals, int iterations, double kappa) {
    printf("\nConvergence Analysis:\n");
    printf("Theoretical convergence factor: %.4f\n", 
           (sqrt(kappa) - 1) / (sqrt(kappa) + 1));
    
    // Calculate observed convergence rate from last few iterations
    if (iterations > 10) {
        double avg_rate = 0.0;
        int count = 0;
        
        for (int i = iterations - 10; i < iterations; i++) {
            if (residuals[i] > 0 && residuals[i+1] > 0) {
                double rate = residuals[i+1] / residuals[i];
                avg_rate += rate;
                count++;
            }
        }
        
        if (count > 0) {
            avg_rate /= count;
            printf("Observed convergence factor (last 10 iterations): %.4f\n", avg_rate);
        }
    }
}

// Function to save residuals to file for plotting
void save_residuals(double *residuals, int iterations, int n) {
    char filename[100];
    sprintf(filename, "residuals_n%d.txt", n);
    FILE *fp = fopen(filename, "w");
    
    for (int i = 0; i <= iterations; i++) {
        fprintf(fp, "%d %.15e\n", i, residuals[i]);
    }
    
    fclose(fp);
    printf("Residuals saved to %s\n", filename);
}

int main() {
    // Test with different matrix sizes
    int sizes[] = {100, 1000, 10000, 100000};
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);
    
    printf("CG Convergence Study for Dense Matrix\n");
    printf("------------------------------------\n");
    
    for (int s = 0; s < num_sizes; s++) {
        int n = sizes[s];
        printf("\nMatrix size: %d x %d\n", n, n);
        
        // Allocate memory
        double **A = (double **)malloc(n * sizeof(double *));
        for (int i = 0; i < n; i++) {
            A[i] = (double *)malloc(n * sizeof(double));
        }
        double *b = (double *)malloc(n * sizeof(double));
        double *x = (double *)malloc(n * sizeof(double));
        
        // Max iterations should be sufficient for convergence
        int max_iter = 5 * n;
        double *residuals = (double *)malloc((max_iter + 1) * sizeof(double));
        
        // Set up the system
        generate_matrix_A(A, n);
        generate_vector_b(b, n);
        
        // Use relative tolerance based on machine epsilon
        double reltol = sqrt(DBL_EPSILON);
        double abstol = 0.0;
        
        printf("Solving with reltol = %.15e, abstol = %.15e\n", reltol, abstol);
        
        // Solve using CG and measure time
        clock_t start = clock();
        int iterations = conjugate_gradient_with_residuals(A, b, x, n, 
                                                          reltol, abstol,
                                                          max_iter, residuals);
        clock_t end = clock();
        
        double cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
        
        printf("CG iterations: %d\n", iterations);
        printf("Time to solution: %.6f seconds\n", cpu_time);
        printf("Final residual: %.15e\n", residuals[iterations]);
        
        // For this matrix, theoretical condition number would depend on eigenvalues
        // For simplicity, we'll use approximations for different n values
        // In practice, you would compute this with a specialized library
        double kappa;
        if (n == 100) kappa = 36.8; // Example value, should be computed properly
        else if (n == 1000) kappa = 208.0;
        else if (n == 10000) kappa = 1124.0;
        else kappa = 3560.0;
        
        printf("Approximated condition number: %.1f\n", kappa);
        
        // Analyze convergence
        analyze_convergence(residuals, iterations, kappa);
        
        // Save residuals for plotting
        save_residuals(residuals, iterations, n);
        
        // Free memory
        for (int i = 0; i < n; i++) {
            free(A[i]);
        }
        free(A);
        free(b);
        free(x);
        free(residuals);
        
        printf("\n");
    }
    
    return 0;
}