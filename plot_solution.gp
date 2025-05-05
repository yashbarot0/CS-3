set term png size 800,600
set output "solution_plot.png"
set title "Solution of Poisson Equation"
set xlabel "x"
set ylabel "y"
set zlabel "u(x,y)"
set pm3d
set view 45,45
splot "solution_N64.txt" using 1:2:3 with pm3d title "u(x,y)"
