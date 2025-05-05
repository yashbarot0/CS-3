# Conjugate Gradient (CG) Solver - Detailed Solution

## 3.1 Basics: The Poisson Problem

We need to solve the Poisson equation on a unit square:
- Domain: $\Omega = (0, 1)^2$
- Equation: $-\Delta u(x) = f(x)$ on the interior of $\Omega$
- Boundary condition: $u(x) = 0$ on $\partial\Omega$ (the boundary)
- Where $\Delta u = \frac{\partial^2 u}{\partial x_1^2} + \frac{\partial^2 u}{\partial x_2^2}$ is the Laplacian operator

### Discretization
To solve this numerically, we discretize the domain into a regular grid:
- We divide the domain into an $N \times N$ grid, giving $N+1$ points in each direction
- Grid spacing: $h = \frac{1}{N}$
- Grid points: $(x_1, x_2) = (ih, jh)$ where $i,j = 0,1,\ldots,N$

Due to the boundary condition $u = 0$ on the boundary, we only need to solve for the interior points where $i,j = 1,2,\ldots,N-1$. This gives us $(N-1) \times (N-1)$ unknown values.

### Finite Difference Approximation
For the Laplacian, we use the standard central difference approximations:

$$\frac{\partial^2 u}{\partial x_1^2} \approx \frac{u(x_1+h,x_2) - 2u(x_1,x_2) + u(x_1-h,x_2)}{h^2}$$

$$\frac{\partial^2 u}{\partial x_2^2} \approx \frac{u(x_1,x_2+h) - 2u(x_1,x_2) + u(x_1,x_2-h)}{h^2}$$

Substituting into the Poisson equation $-\Delta u = f$:

$$-\left(\frac{u(x_1+h,x_2) - 2u(x_1,x_2) + u(x_1-h,x_2)}{h^2} + \frac{u(x_1,x_2+h) - 2u(x_1,x_2) + u(x_1,x_2-h)}{h^2}\right) = f(x_1,x_2)$$

Simplifying:

$\frac{4u(x_1,x_2) - u(x_1+h,x_2) - u(x_1-h,x_2) - u(x_1,x_2+h) - u(x_1,x_2-h)}{h^2} = -f(x_1,x_2)$

This gives us the discrete form of the Poisson equation at each interior grid point.

### Linear System Formulation
We can organize this into a linear system $Ay = b$ where:
- $y$ is a vector containing the unknown values of $u$ at all interior grid points
- $A$ is the coefficient matrix representing the finite difference approximation of the Laplacian
- $b$ is the right-hand side vector derived from $f$

#### Ordering of Grid Points
For the ordering of grid points, we'll use a row-major ordering (left to right, then top to bottom):
- Grid point $(i,j)$ (where $i,j = 1,2,\ldots,N-1$) is mapped to index $k = (j-1)(N-1) + (i-1)$
- This gives us a one-dimensional vector of $(N-1)^2$ unknowns

#### Structure of Matrix A
For each interior point $(i,j)$ with corresponding index $k$:
- The diagonal element $A_{k,k} = 4$ (coefficient of the central point)
- If the point has a left neighbor, $A_{k,k-1} = -1$
- If the point has a right neighbor, $A_{k,k+1} = -1$
- If the point has a bottom neighbor, $A_{k,k-(N-1)} = -1$
- If the point has a top neighbor, $A_{k,k+(N-1)} = -1$

The resulting matrix $A$ is sparse with at most 5 non-zero entries per row, symmetric, and positive definite.

#### Right-hand Side Vector
The right-hand side vector $b$ is computed as:
$b_k = -f(x_1,x_2) \cdot h^2$

where $(x_1,x_2)$ are the coordinates of the grid point corresponding to index $k$.

## 3.2 Serial Implementation of CG Algorithm

The Conjugate Gradient (CG) algorithm for solving $Ax = b$ is implemented as follows:

1. Initialize:
   - $x_0 = 0$ (initial guess)
   - $r_0 = b - Ax_0 = b$ (since $x_0 = 0$)
   - $p_0 = r_0$ (initial search direction)

2. Iterate until convergence:
   - Compute $\alpha_k = \frac{r_k^T r_k}{p_k^T A p_k}$
   - Update solution: $x_{k+1} = x_k + \alpha_k p_k$
   - Update residual: $r_{k+1} = r_k - \alpha_k A p_k$
   - Compute $\beta_k = \frac{r_{k+1}^T r_{k+1}}{r_k^T r_k}$
   - Update search direction: $p_{k+1} = r_{k+1} + \beta_k p_k$

3. Check convergence:
   - Stop when $\|r_k\| < \text{tol} \cdot \|r_0\|$

### Efficient Implementation Details
To optimize performance:-
1. **Matrix-free Implementation**: Instead of explicitly storing the full matrix $A$, the code uses the stencil directly in the matrix-vector multiplication function `poisson_matrix_vector_mult()`. This saves memory and computation time.
2. **Reusing Dot Products**: The dot product $r_k^T r_k$ is calculated once and reused for both $\alpha_k$ and $\beta_k$.
3. **Memory Allocation**: The code allocates memory for vectors only once before the iteration starts, preventing repeated memory allocation/deallocation in the loop.
4. **Relative Residual Stopping Criterion**: The algorithm stops when the relative residual $\|r_k\|/\|r_0\|$ falls below the specified tolerance, which is more robust than using an absolute tolerance.

### Performance Results
The code solves the Poisson problem for the function $f(x,y) = 2\pi^2\sin(\pi x)\sin(\pi y)$ with increasing grid sizes. The exact solution is $u(x,y) = \sin(\pi x)\sin(\pi y)$.

For grid sizes N = 8, 16, 32, 64, 128, 256, the code tabulates:
- Number of CG iterations
- Time to solution
- Maximum error compared to the exact solution


## 3.3 Convergence of CG

### Dense Linear System
In this part, we study the convergence of CG for a dense linear system $Ay = b$ where:
- $A \in \mathbb{R}^{N \times N}$ with entries $(A)_{ij} = \frac{N - |i-j|}{N}$
- $b \in \mathbb{R}^N$ with all entries equal to 1

### Stopping Criterion
We implement the stopping criterion:
$\|r_k\| \leq \max(\text{reltol} \cdot \|r_0\|, \text{abstol})$

with:
- absolute tolerance: abstol = 0
- relative tolerance: reltol = $\sqrt{\epsilon_{\text{double}}}$

where $\epsilon_{\text{double}}$ is the machine epsilon for double precision.

### Convergence Analysis
The theoretical convergence bound for CG is:

$\|e_k\| \leq 2 \left(\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}\right)^k \|e_0\|$

where $\kappa = \kappa(A) = \frac{\lambda_{\max}}{\lambda_{\min}}$ is the spectral condition number.

The code computes and saves the residual norm at each iteration to study the convergence rate. It then compares the observed convergence rate with the theoretical bound.

For each matrix size N = 10², 10³, 10⁴, and 10⁵, the code:
1. Generates the dense matrix $A$ and right-hand side $b$
2. Solves the system using CG
3. Records the residual at each iteration
4. Compares the observed convergence rate with the theoretical prediction

---
## File Structure
1. poisson_fd.c - Sets up the Poisson problem as a linear system (Question 3.1)
2. cg_solver.c - Serial implementation of CG for the Poisson problem (Question 3.2)
3. cg_convergence.c - Study of CG convergence for a dense linear system (Question 3.3)
4. Makefile - Compiles and runs all code parts
### Compile all programs
make all
### Run each part individually:
make run_q1  # For Question 3.1
make run_q2  # For Question 3.2
make run_q3  # For Question 3.3

###  Run all parts and generate plots (requires gnuplot):
make run_all

## Part 3.1: The Poisson Problem

The code in `poisson_fd.c` demonstrates how to convert the Poisson problem into a linear system Ay = b. Key features:

- Uses finite difference approximation for the Laplacian
- Sets up a structured grid on the unit square
- Creates the sparse matrix A and right-hand side vector b
- For N=4, displays the full matrix and vector for verification

The main concept here is that the 5-point stencil for the Laplacian creates a sparse matrix with a specific pattern: diagonal entries are 4/h², and off-diagonal entries for neighboring grid points are -1/h². 

## Part 3.2: Serial Implementation of CG

The code in `cg_solver.c` implements the Conjugate Gradient method to solve the Poisson problem. Key features:
- Matrix-free implementation (no explicit storage of A)
- Solves for increasing grid sizes: N = 8, 16, 32, 64, 128, 256
- Measures performance (iterations and time)
- Calculates error compared to the exact solution
- Saves the solution for visualization

The implementation is optimized by:
- Avoiding explicit storage of the matrix
- Minimizing memory allocation/deallocation
- Reusing dot products where possible

## Part 3.3: Convergence of CG

The code in `cg_convergence.c` studies the convergence properties of CG for a dense linear system. Key features:
- Creates the dense matrix A with entries (A)ᵢⱼ = (N - |i-j|)/N
- Uses stopping criterion with relative tolerance = √eps(double)
- Records residuals at each iteration
- Compares observed convergence with theoretical prediction
- Tests for matrix sizes N = 10², 10³, 10⁴.

**Note**: For N = 10⁵, Not able to compute on Callan system.

## Visualization

After running `make run_all`, two plot files will be generated:
1. solution_plot.png - 3D plot of the solution to the Poisson equation
2. residuals_plot.png - Convergence history of residuals for different matrix sizes

## Understanding the Results

For Question 3.2, observe how:
![[solution_plot.png]]

- The number of iterations grows with grid size
- The time to solution scales with problem size
- The error decreases with finer grids

For Question 3.3, observe how:

![[residuals_plot 1.png]]
- The convergence rate depends on the condition number
- The residual decreases approximately linearly on a log scale
- The observed convergence matches the theoretical prediction

## Important Notes
1. For larger problem sizes (especially N = 10⁵ in Q3.3), we may need substantial memory and computation time. (Not able to execute on Callan, process getting killed automatically).
2. The plotting scripts require gnuplot.
