set term png size 800,600
set output "residuals_plot.png"
set title "CG Convergence"
set xlabel "Iteration"
set ylabel "Residual (log scale)"
set logscale y
set grid
plot "residuals_n100.txt" using 1:2 with lines title "N=100", \
     "residuals_n1000.txt" using 1:2 with lines title "N=1000", \
     "residuals_n10000.txt" using 1:2 with lines title "N=10000"
