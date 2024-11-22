set title "H_{2} Vibrational Wavefunctions N=100"
set xlabel "R (a.u)"
set ylabel ""
set yrange[-0.61:-0.5]
set xrange[0:6]
set terminal pngcairo
set output "vibfns.png"
set grid
plot for [i=2:13] "vibfns_l0.txt" using 1:i with lines notitle, \
    "PEC.1ssg" using 1:2 w l title "1s{/Symbol s}_{g}"
