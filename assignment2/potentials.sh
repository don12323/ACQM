#!/bin/bash
make clean
make main
rm Energies.txt


for R in 0.1 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.5 7.0 7.5 8.0
do
	sed s/RRRR/${R}/ input2.txt > input.txt
	./main
		

done
