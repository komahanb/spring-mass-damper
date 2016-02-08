echo "set xlabel 'Time'; set ylabel 'Solution X'; set term png; set output 'solution.png'; plot 'solution.dat' u 1:3 w l title 'BDF', 'solution.dat' u 1:8 w l title 'Exact'" | gnuplot
