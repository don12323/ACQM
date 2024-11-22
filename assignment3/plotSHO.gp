set terminal pngcairo enhanced size 800,600
set output 'SHO.png'

# Set titles and labels
set title " Quantum Harmonic Oscillator {/Symbol d}x=1"
set xlabel "x(a_{0})"
set ylabel "Energy(Ha)"
set xrange[-5:5]
set yrange[0:4]

# Set grid
set grid

# Define the analytical solutions, scaled by 0.75 and shifted up by their energies
f0(x) = 0.75*pi**(-1.0/4.0) * exp(-x**2 / 2.0) +0.5 
f1(x) =-0.75*pi**(-1.0/4.0) * sqrt(2) * x * exp(-x**2 / 2.0) + 1.5
f2(x) = -0.75*pi**(-1.0/4.0) * (1.0/sqrt(2.0)) * (2.0 * x**2 - 1.0) * exp(-x**2 / 2.0) + 2.5
f3(x) = 0.75*pi**(-1.0/4.0) * (1.0/sqrt(3.0)) * (2.0 * x**3 - 3.0 * x) * exp(-x**2 / 2.0) + 3.5

plot 'output.txt' using 1:2 w l title 'V(x)',\
     'wf_0.txt' using 1:2 w l title 'n=0',\
     'wf_1.txt' using 1:2 w l title 'n=1',\
     'wf_2.txt' using 1:2 w l title 'n=2',\
     'wf_3.txt' using 1:2 w l title 'n=3', \
     f0(x) with lines dt 2 lc rgb 'black' title 'Analytical', \
     f1(x) with lines dt 2 lc rgb 'black' notitle, \
     f2(x) with lines dt 2 lc rgb 'black' notitle, \
     f3(x) with lines dt 2 lc rgb 'black' notitle	

