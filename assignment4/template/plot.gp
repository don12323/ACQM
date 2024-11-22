set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'DCS.png'
set title "Differential cross section for Electron scattering E=50eV"
set xlabel "Scattering Angle(deg)"
set ylabel "Differential cross section (a_{0}^2/Sr)"
set grid

plot 'dcs_E50.txt' using 1:2 w l title "E=50.0eV", \

set output 'ICS_pos.png'
set title "Integrated Cross section for Positron scattering"
set xlabel "Energy(eV)"
set ylabel " Integrated Cross Section (a_{0}^2)"
set logscale y
set yrange[1.0e-4: 1.0e+4]
datafile='ics_results.txt'
stats datafile u 0 noout
plot for [i=2:STATS_columns] datafile using 1:i with lines title columnhead


set output 'total_ICS_pos.png'

# Set title and labels for total ICS
set title "Total Integrated Cross Section for Positron scattering"
set xlabel "Incident Particle Energy (eV)"
set ylabel "Total Integrated Cross Section (a_{0}^2)"
unset yrange
unset logscale y
# Generate dynamic sum of all columns for total ICS
sum_columns = ""
do for [i=2:STATS_columns] {
    sum_columns = sum_columns . (i > 2 ? "+" : "") . "column(" . i . ")"
}

# Plot total ICS by summing all y-value columns
plot datafile using 1:(@sum_columns) with lines title 'Total ICS'
