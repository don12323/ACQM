set terminal pngcairo enhanced font 'Arial,10' size 600,900
set output 'Elvls_alp_1_L0.png'

# Set title and labels
set title "Energy Level Diagram"
unset xlabel
set ylabel "Eigen-Energies"
set nokey

# Set the range of the plot
set xrange [3:55]
set style line 1 pointtype 7 pointsize 0.5
# Set up multiplot for positive and negative values
set multiplot layout 2,1

# Top plot with positive values
set title "Energy Diagram"
set logscale y
plot 'E4.txt' using 1:($2 > 0 ? $2 : 1/0) with points ls 1, \
     'E5.txt' using 1:($2 > 0 ? $2 : 1/0) with points ls 1, \
     'E10.txt' using 1:($2 > 0 ? $2 : 1/0) with points ls 1, \
     'E20.txt' using 1:($2 > 0 ? $2 : 1/0) with points ls 1, \
     'EN30l0.txt' using 1:($2 > 0 ? $2 : 1/0) with points ls 1, \
     'EN40l0.txt' using 1:($2 > 0 ? $2 : 1/0) with points ls 1, \
     'EN50l0.txt' using 1:($2 > 0 ? $2 : 1/0) with points ls 1, \
     'EN4l0a1.50.txt' using 1:($2 > 0 ? $2 : 1/0) with points ls 1 pointtype 6 linecolor rgb "blue",\
     'EN5l0a1.50.txt' using 1:($2 > 0 ? $2 : 1/0) with points ls 1 pointtype 6 linecolor rgb "blue",\
     'EN10l0a1.50.txt' using 1:($2 > 0 ? $2 : 1/0) with points ls 1 pointtype 6 linecolor rgb "blue",\
     'EN20l0a1.50.txt' using 1:($2 > 0 ? $2 : 1/0) with points ls 1 pointtype 6 linecolor rgb "blue",\
     'EN30l0a1.50.txt' using 1:($2 > 0 ? $2 : 1/0) with points ls 1 pointtype 6 linecolor rgb "blue",\
     'EN40l0a1.50.txt' using 1:($2 > 0 ? $2 : 1/0) with points ls 1 pointtype 6 linecolor rgb "blue",\
     'EN50l0a1.50.txt' using 1:($2 > 0 ? $2 : 1/0) with points ls 1 pointtype 6 linecolor rgb "blue"

# Bottom plot with negative values
unset title
set yrange [1:0.0003] reverse
set logscale y
set label "ASDASD"
set xlabel "Number of Basis Functions N"
set format y '-10^{%T}'

plot 'E4.txt' using 1:($2 < 0 ? abs($2) : 1/0) with points ls 1, \
     'E5.txt' using 1:($2 < 0 ? abs($2) : 1/0) with points ls 1, \
     'E10.txt' using 1:($2 < 0 ? abs($2) : 1/0) with points ls 1, \
     'E20.txt' using 1:($2 < 0 ? abs($2) : 1/0) with points ls 1, \
     'EN30l0.txt' using 1:($2 < 0 ? abs($2) : 1/0) with points ls 1, \
     'EN40l0.txt' using 1:($2 < 0 ? abs($2) : 1/0) with points ls 1, \
     'EN50l0.txt' using 1:($2 < 0 ? abs($2) : 1/0) with points ls 1, \
     'EN4l0a1.50.txt' using 1:($2 < 0 ? abs($2) : 1/0) with points ls 1 pointtype 6 linecolor rgb "blue",\
     'EN5l0a1.50.txt' using 1:($2 < 0 ? abs($2) : 1/0) with points ls 1 pointtype 6 linecolor rgb "blue",\
     'EN10l0a1.50.txt' using 1:($2 < 0 ? abs($2) : 1/0) with points ls 1 pointtype 6 linecolor rgb "blue",\
     'EN20l0a1.50.txt' using 1:($2 < 0 ? abs($2) : 1/0) with points ls 1 pointtype 6 linecolor rgb "blue",\
     'EN30l0a1.50.txt' using 1:($2 < 0 ? abs($2) : 1/0) with points ls 1 pointtype 6 linecolor rgb "blue",\
     'EN40l0a1.50.txt' using 1:($2 < 0 ? abs($2) : 1/0) with points ls 1 pointtype 6 linecolor rgb "blue",\
     'EN50l0a1.50.txt' using 1:($2 < 0 ? abs($2) : 1/0) with points ls 1 pointtype 6 linecolor rgb "blue" 

set key outside
#plot NaN with points ls 1 linecolor rgb "red" title "Alpha = 1", \
    #NaN with points ls 1 pointtype 6 linecolor rgb "blue" title "Alpha = 1.5"
unset multiplot

