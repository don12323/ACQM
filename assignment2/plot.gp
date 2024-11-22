set title "Potential energy curves H_{2}. l_{max}=5, N=10, \alpha =2.0 "
set xrange[0:10]
set yrange[-1:2]
set ylabel "Energy(Ha)"
set xlabel "R(a_{0})"
plot for [i=2:3] "Energies.txt" using 1:i w p notitle,\
     "PEC.1ssg" using 1:2 w l title "1s{/Symbol s}_{g}",\
     "PEC.2psu" using 1:2 w l title "2p{/Symbol s}_{u}"
pause -1
