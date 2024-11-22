set title "H_{2} Vibrational Wavefunctions N=100"
set xlabel "R (a.u)"
set ylabel "Energy(Ha)"
set yrange[-0.61:-0.5]
set xrange[0:6]
set terminal pngcairo
set output "vibfns.png"
set grid
plot for [i=2:13] "vibfns_l0.txt" using 1:i with lines notitle, \
    "PEC.1ssg" using 1:2 w l title "1s{/Symbol s}_{g}"

set output "veq0wfs.png"
set title "v=0 vibrational wavefunction"
unset yrange
set xrange[0:4.5] 
plot "mu_2748.46079.txt" using 1:(abs($2)) w l title "T_{2}", \
     "mu_918.07635.txt" using 1:(abs($2)) w l title "H_{2}", \
     "mu_1223.89925.txt" using 1:(abs($2)) w l title "HD", \
     "mu_1376.39236.txt" using 1:(abs($2)) w l title "HT", \
     "mu_1835.24151.txt" using 1:(abs($2)) w l title "D_{2}", \
     "mu_2200.87999.txt" using 1:(abs($2)) w l title "DT"
 
