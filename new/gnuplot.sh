echo "set xlabel 'Time'; set ylabel 'q'; set term png; set output 'stage$1.png'; plot 'erk$1.dat' u 1:2 w lp, 'irk$1.dat' u 1:2 w lp, 'dirk$1.dat' u 1:2 w lp" | gnuplot
