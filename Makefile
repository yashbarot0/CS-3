# Makefile for CG solver assignment

CC = gcc
CFLAGS = -Wall -O3 -march=native -ffast-math -fopenmp
LDFLAGS = -lm

# Main targets
all: poisson_fd cg_solver cg_convergence

# Poisson finite difference setup (Question 3.1)
poisson_fd: poisson_fd.c
	$(CC) $(CFLAGS) -o poisson_fd poisson_fd.c $(LDFLAGS)

# CG solver for Poisson problem (Question 3.2)
cg_solver: cg_solver.c
	$(CC) $(CFLAGS) -o cg_solver cg_solver.c $(LDFLAGS)

# CG convergence study (Question 3.3)
cg_convergence: cg_convergence.c
	$(CC) $(CFLAGS) -o cg_convergence cg_convergence.c $(LDFLAGS)

# Generate plotting script
plot_script:
	@echo 'set term png size 800,600' > plot_solution.gp
	@echo 'set output "solution_plot.png"' >> plot_solution.gp
	@echo 'set title "Solution of Poisson Equation"' >> plot_solution.gp
	@echo 'set xlabel "x"' >> plot_solution.gp
	@echo 'set ylabel "y"' >> plot_solution.gp
	@echo 'set zlabel "u(x,y)"' >> plot_solution.gp
	@echo 'set pm3d' >> plot_solution.gp
	@echo 'set view 45,45' >> plot_solution.gp
	@echo 'splot "solution_N64.txt" using 1:2:3 with pm3d title "u(x,y)"' >> plot_solution.gp
	@echo 'Created plotting script plot_solution.gp'
	
	@echo 'set term png size 800,600' > plot_residuals.gp
	@echo 'set output "residuals_plot.png"' >> plot_residuals.gp
	@echo 'set title "CG Convergence"' >> plot_residuals.gp
	@echo 'set xlabel "Iteration"' >> plot_residuals.gp
	@echo 'set ylabel "Residual (log scale)"' >> plot_residuals.gp
	@echo 'set logscale y' >> plot_residuals.gp
	@echo 'set grid' >> plot_residuals.gp
	@echo 'plot "residuals_n100.txt" using 1:2 with lines title "N=100", \' >> plot_residuals.gp
	@echo '     "residuals_n1000.txt" using 1:2 with lines title "N=1000", \' >> plot_residuals.gp
	@echo '     "residuals_n10000.txt" using 1:2 with lines title "N=10000"' >> plot_residuals.gp
	@echo 'Created plotting script plot_residuals.gp'

# Clean up
clean:
	rm -f poisson_fd cg_solver cg_convergence
	rm -f *.o *.txt *.png *.gp

# Run all tests (warning: cg_convergence for N=10^5 can take a long time)
run_all: all plot_script
	@echo "Running Poisson finite difference setup..."
	./poisson_fd
	@echo "\nRunning CG solver for Poisson problem..."
	./cg_solver
	@echo "\nRunning CG convergence study..."
	./cg_convergence
	@echo "\nGenerating plots (requires gnuplot)..."
	gnuplot plot_solution.gp
	gnuplot plot_residuals.gp
	@echo "Done. Check solution_plot.png and residuals_plot.png for results."

# Run individual parts
run_q1: poisson_fd
	./poisson_fd

run_q2: cg_solver
	./cg_solver

run_q3: cg_convergence
	./cg_convergence

.PHONY: all clean run_all run_q1 run_q2 run_q3 plot_script