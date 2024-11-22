set title "H_{2} Dissociative Wavefunctions"
set xlabel "R (a.u)"
set ylabel "Energy (Ha)"
set yrange[-0.61:0]
set xrange[0:10]
set terminal pngcairo
set output "diswfns.png"
set grid
set key outside right
plot for [i=2:9] "contwfs.txt" using 1:i with lines notitle, \
     "PEC.2psu" using 1:2 w l title "2p{/Symbol s}_{u}", \
     "PEC.1ssg" using 1:2 w l title "1s{/Symbol s}_{g}"
set output "franckCon.png"
set title "Franck-Condon Apporoximation"
set ylabel "KER distribution (arb.units)"
set xlabel "Kinetic Energy release(eV)"
set xrange[0:20]
unset yrange
plot "KED.txt" using 1:2 with lines title "v_{i} = 0", \
     "KED.txt" using 1:3 with lines title "v_{i} = 3", \
     "KED.txt" using 1:4 with lines title "v_{i} = 6", \
     "KED.txt" using 1:5 with lines title "v_{i} = 9"

                       
