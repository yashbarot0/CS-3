#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Function f(x,y) = 2π² sin(πx) sin(πy)
double f(double x, double y) {
    return 2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
}

// Exact solution u(x,y) = sin(πx) sin(πy)
double exact_solution(double x, double y) {
    return sin(M_PI * x) * sin(M_PI * y);
}

// Function to calculate the norm of a vector
double vector_norm(double *v, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}

// Function to perform matrix-vector multiplication for Poisson problem
// Instead of storing the full matrix, we use the stencil directly
void poisson_matrix_vector_mult(double *result, double *v, int N) {
    int n = (N-1) * (N-1);
    
    // Initialize result to zero
    for (int i = 0; i < n; i++) {
        result[i] = 0.0;
    }
    
    // Compute A*v without explicitly forming A
    for (int i = 0; i < N-1; i++) {
        for (int j = 0; j < N-1; j++) {
            int idx = i + j * (N-1);
            
            // Apply stencil: 4*center - neighbors
            result[idx] = 4.0 * v[idx];
            
            // Left neighbor
            if (i > 0) {
                result[idx] -= v[idx-1];
            }
            
            // Right neighbor
            if (i < N-2) {
                result[idx] -= v[idx+1];
            }
            
            // Bottom neighbor
            if (j > 0) {
                result[idx] -= v[idx-(N-1)];
            }
            
            // Top neighbor
            if (j < N-2) {
                result[idx] -= v[idx+(N-1)];
            }
        }
    }
}

// Conjugate Gradient method
int conjugate_gradient(int N, double *b, double *x, double tol, int max_iter) {
    int n = (N-1) * (N-1);
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
    
    int iter;
    for (iter = 0; iter < max_iter && r_norm > tol * initial_r_norm; iter++) {
        // Compute A*p
        poisson_matrix_vector_mult(Ap, p, N);
        
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
    }
    
    // Free memory
    free(r);
    free(p);
    free(Ap);
    
    return iter;
}

// Function to solve the Poisson problem and measure performance
void solve_poisson_problem(int N) {
    int n = (N-1) * (N-1);
    double h = 1.0 / N;
    double h2 = h * h;
    
    // Allocate memory
    double *b = (double *)malloc(n * sizeof(double));
    double *x = (double *)malloc(n * sizeof(double));
    
    // Set up right-hand side vector b
    for (int j = 0; j < N-1; j++) {
        for (int i = 0; i < N-1; i++) {
            int idx = i + j * (N-1);
            double x_pos = (i+1) * h;
            double y_pos = (j+1) * h;
            b[idx] = -f(x_pos, y_pos) * h2;
        }
    }
    
    // Solve using CG
    double tol = 1e-8;
    int max_iter = 10000;
    
    clock_t start = clock();
    int iterations = conjugate_gradient(N, b, x, tol, max_iter);
    clock_t end = clock();
    
    double cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    printf("Grid size N = %d, Matrix size = %d x %d\n", N, n, n);
    printf("CG iterations: %d\n", iterations);
    printf("Time to solution: %.6f seconds\n", cpu_time);
    
    // Calculate error compared to exact solution
    double max_error = 0.0;
    for (int j = 0; j < N-1; j++) {
        for (int i = 0; i < N-1; i++) {
            int idx = i + j * (N-1);
            double x_pos = (i+1) * h;
            double y_pos = (j+1) * h;
            double exact = exact_solution(x_pos, y_pos);
            double error = fabs(x[idx] - exact);
            if (error > max_error) {
                max_error = error;
            }
        }
    }
    printf("Maximum error: %.6e\n\n", max_error);
    
    // Free memory
    free(b);
    free(x);
}

// Function to save solution for plotting
void save_solution(int N, double *x) {
    char filename[100];
    sprintf(filename, "solution_N%d.txt", N);
    FILE *fp = fopen(filename, "w");
    
    double h = 1.0 / N;
    
    // Include boundary points (which are zero)
    for (int j = 0; j <= N; j++) {
        double y = j * h;
        for (int i = 0; i <= N; i++) {
            double x_pos = i * h;
            double u_val;
            
            if (i == 0 || i == N || j == 0 || j == N) {
                // Boundary point
                u_val = 0.0;
            } else {
                // Interior point
                int idx = (i-1) + (j-1) * (N-1);
                u_val = x[idx];
            }
            
            fprintf(fp, "%lf %lf %lf\n", x_pos, y, u_val);
        }
        fprintf(fp, "\n"); // Add blank line for gnuplot
    }
    
    fclose(fp);
    printf("Solution saved to %s\n", filename);
}

int main() {
    // Solve for increasing grid sizes
    int grid_sizes[] = {8, 16, 32, 64, 128, 256};
    int num_sizes = sizeof(grid_sizes) / sizeof(grid_sizes[0]);
    
    printf("Performance of CG solver for Poisson problem\n");
    printf("-------------------------------------------\n");
    
    for (int i = 0; i < num_sizes; i++) {
        solve_poisson_problem(grid_sizes[i]);
    }
    
    // For visualization, solve and save the N=64 solution
    int N = 64;
    int n = (N-1) * (N-1);
    double h = 1.0 / N;
    double h2 = h * h;
    
    double *b = (double *)malloc(n * sizeof(double));
    double *x = (double *)malloc(n * sizeof(double));
    
    // Set up right-hand side vector b
    for (int j = 0; j < N-1; j++) {
        for (int i = 0; i < N-1; i++) {
            int idx = i + j * (N-1);
            double x_pos = (i+1) * h;
            double y_pos = (j+1) * h;
            b[idx] = -f(x_pos, y_pos) * h2;
        }
    }
    
    // Solve using CG
    double tol = 1e-8;
    int max_iter = 10000;
    conjugate_gradient(N, b, x, tol, max_iter);
    
    // Save solution for plotting
    save_solution(N, x);
    
    // Free memory
    free(b);
    free(x);
    
    return 0;
}