set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
output_dir = '/mnt/c/Users/Imesh/Desktop/honours/ACQM/assignment5/images/'
data = '/mnt/c/Users/Imesh/Desktop/honours/ACQM/assignment5/data/'
set output output_dir . 'V1s10ev.png'


set title "Fully Off Shell Vmatrix Elements E=10eV"
set xlabel "k'"
set ylabel "k"
set zlabel "V(k', k)"
set grid
set view 65, 38

set view 60, 60
set hidden3d
set xyplane at -1 

splot [0:5][0:5][-0.8:1] data . 'Vdir_1s1s.txt' using 1:2:3 with lines lc rgb "red" title "V-direct 1s-1s", \
      data . 'Vex_1s1s.txt' using 1:2:3 with lines ls 2 lc rgb "blue" title "V-exchange 1s-1s", \

set output output_dir . 'V2s10ev.png'
set view 60, 60
splot [0:5][0:5][-0.8:1] data . 'Vdir_1s2s.txt' using 1:2:3 with lines lc rgb "red" title "V-direct 1s-2s", \
      data . 'Vex_1s2s.txt' using 1:2:3 with lines ls 2 lc rgb "blue" title "V-exchange 1s-2s", \

set view 60, 60
set output output_dir . 'V3s10ev.png'
splot [0:5][0:5][-0.8:1] data . 'Vdir_1s3s.txt' using 1:2:3 with lines lc rgb "red" title "V-direct 1s-3s", \
      data . 'Vex_1s3s.txt' using 1:2:3 with lines lc rgb "blue" title "V-exchange 1s-3s", \

#PLOT 1eV
set title "Fully Off Shell Vmatrix Elements E=1eV"
set output output_dir . 'V1s1ev.png'
splot [0:5][0:5][-0.8:1] data . 'Vdir_1s1s1ev.txt' using 1:2:3 with lines lc rgb "red" title "V-direct 1s-1s", \
      data . 'Vex_1s1s1ev.txt' using 1:2:3 with lines ls 2 lc rgb "blue" title "V-exchange 1s-1s", \

set output output_dir . 'V2s1ev.png'
splot [0:5][0:5][-0.8:1] data . 'Vdir_1s2s1ev.txt' using 1:2:3 with lines lc rgb "red" title "V-direct 1s-2s", \
      data . 'Vex_1s2s1ev.txt' using 1:2:3 with lines ls 2 lc rgb "blue" title "V-exchange 1s-2s", \

set output output_dir . 'V3s1ev.png'
splot [0:5][0:5][-0.8:1] data . 'Vdir_1s3s1ev.txt' using 1:2:3 with lines lc rgb "red" title "V-direct 1s-3s", \
      data . 'Vex_1s3s1ev.txt' using 1:2:3 with lines lc rgb "blue" title "V-exchange 1s-3s", \


#------------------------------------------------------------------
set output output_dir . 'OS1s10ev.png'
set title "On shell values"
set xlabel "incident momentum k"
set ylabel "On shell V matrix elements"
direct_1s1s(k) = (2.0/pi) * (-1.0/4.0 * k**2 / (k**2 + 1.0) - 1.0/4.0 * log(1.0 + k**2)) 
direct_1s2s(k) = (2.0/pi)*(16.0 * k * (4.0 * k**2 + 3.0) * sqrt(8.0 * k**2 - 6.0)) / (81.0 * (4.0 * k**2 + 1.0)**2)
direct_1s3s(k) = (2.0/pi) *(9.0 * sqrt(3.0)*k* (135.0 * k**4 + 87.0 * k**2 - 4.0) * sqrt(9.0 * k**2 - 8.0)) / (128.0 * (9.0 * k**2 + 1.0)**3)
exchange_1s1s(k) = -k**2 * (2.0/pi) *(k**2 - 3.0) / (k**2 + 1.0)**3
exchange_1s2s(k) = -((2.0/pi) * 16.0 * k * (16.0 * k**4 - 72.0 * k**2 + 13.0) * sqrt(8.0 * k**2 - 6.0)) / (9.0 * (4.0 * k**2 + 1.0)**4)
exchange_1s3s(k) =-(2.0/pi) * (9.0 * sqrt(3.0)* k * (1701.0 * k**6 - 8208.0 * k**4 + 2325.0 * k**2 - 70.0) * sqrt(9.0 * k**2 - 8.0)) / (8.0 * (9.0 * k**2 + 1.0)**5)


plot data . 'OnSh_dir1s1s.txt' using 1:2 with points pt 7 ps 0.8 lc rgb 'black' title "Direct 1s-1s", \
     direct_1s1s(x) with lines title "1s-1s direct analytical", \
     data . 'OnSh_ex1s1s.txt' using 1:2 with points pt 7 ps 0.8 title "Exchange 1s-1s", \
     exchange_1s1s(x) with lines title "1s-1s exchange analytical"

set output output_dir . 'OS2s10ev.png'
plot data . 'OnSh_dir1s2s.txt' using 1:2 with points pt 7 ps 0.8 lc rgb 'black' title "Direct 1s-2s", \
     direct_1s2s(x) with lines title "1s-2s direct analytical", \
     data . 'OnSh_ex1s2s.txt' using 1:2 with points pt 7 ps 0.8 title "Exchange 1s-2s", \
     exchange_1s2s(x) with lines title "1s-2s exchange analytical"

set output output_dir . 'OS3s10ev.png'
plot data . 'OnSh_dir1s3s.txt' using 1:2 with points pt 7 ps 0.8 lc rgb 'black' title "Direct 1s-3s", \
     direct_1s3s(x) with lines title "1s-3s direct analytical", \
     data . 'OnSh_ex1s3s.txt' using 1:2 with points pt 7 ps 0.8 title "Exchange 1s-3s", \
     exchange_1s3s(x) with lines title "1s-3s exchange analytical"

set output output_dir . 'wf.png'
plot 'wf.txt' using 1:2 w l notitle


set output output_dir . 'KVthet0.png'
unset title
set ylabel "Triplet 1s-1s matrix elements"
set xlabel "off-shell linear momentum k"
set yrange[-1.5:0.5]

set style line 1 lc rgb "black" pt 6 ps 2    
set style line 2 lc rgb "black" pt 5 ps 1    
plot data . 'k1V_thet0.txt' using 1:2 with lines lc rgb "red" title "V(k)", \
     data . 'k1V_thet0.txt' using 1:2 with points ls 1 title "kgrid 1 theta = 0", \
     data . 'k1K_thet0.txt' using 1:2 with lines lc rgb "blue" title "K(k)", \
     data . 'k1K_thet0.txt' using 1:2 with points ls 1 notitle, \
     data . 'k2V_thet0.txt' using 1:2 with lines lc rgb "red" notitle, \
     data . 'k2V_thet0.txt' using 1:2 with points ls 2 title "kgrid 2 theta = 0", \
     data . 'k2K_thet0.txt' using 1:2 with lines lc rgb "blue" notitle, \
     data . 'k2K_thet0.txt' using 1:2 with points ls 2 notitle


set output output_dir . 'KVthet12.png'
plot data . 'k1V_thet1.txt' using 1:2 with lines lc rgb "red" title "V(k)", \
     data . 'k1V_thet1.txt' using 1:2 with points ls 1 title "kgrid 1 theta = 1", \
     data . 'k1K_thet1.txt' using 1:2 with lines lc rgb "blue" title "K(k)", \
     data . 'k1K_thet1.txt' using 1:2 with points ls 1 notitle, \
     data . 'k2V_thet2.txt' using 1:2 with lines lc rgb "red" notitle, \
     data . 'k2V_thet2.txt' using 1:2 with points ls 2 title "kgrid 2 theta = 2", \
     data . 'k2K_thet2.txt' using 1:2 with lines lc rgb "blue" notitle, \
     data . 'k2K_thet2.txt' using 1:2 with points ls 2 notitle


