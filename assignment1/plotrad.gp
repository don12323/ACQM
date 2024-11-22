set title "Radial Wavefunctions N=50"
set xlabel "r (a.u)"
set ylabel ""

set terminal pngcairo
set output "radfns.png"
set grid
plot "radfns_l0.txt" using 1:2 with lines title "1s", \
     "radfns_l0.txt" using 1:3 with lines title "2s", \
     "radfns_l0.txt" using 1:4 with lines title "3s", \
     "radfns_l1.txt" using 1:2 with lines title "2p", \
     "radfns_l1.txt" using 1:3 with lines title "3p", \
     "radfns_l1.txt" using 1:4 with lines title "4p", \
