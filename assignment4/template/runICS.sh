#!/bin/bash

# Ensure the output file is empty before starting
output_file="ics_results.txt"
echo -n > $output_file
#define parameters
E_start=0.5
E_end=50
dE=0.5
L=3


header="E"
for l in $(seq 0 $L); do
	header="$header l=$l"
done
echo $header >> $output_file

for E in $(seq $E_start $dE $E_end); do
    echo "$E" > data.in
    echo "20000, 0.001" >> data.in
    echo "1, 0, $L" >> data.in  
    
    ./main
    
    ICS_values=$(awk '{print $2}' ics.txt | paste -sd " ")
    echo "$E $ICS_values" >> $output_file
done

