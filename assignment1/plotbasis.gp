set title "Tridiagonal Laguerre basis functions (l=0, alpha=2)"
set xlabel "r (a.u)"
set ylabel ""

set terminal pngcairo
set output "basisfns.png"

plot "eigenfns.txt" using 1:2 with lines title "k=1", \
     "" using 1:3 with lines title "k=2", \
     "" using 1:4 with lines title "k=3", \
     "" using 1:5 with lines title "k=4" 
