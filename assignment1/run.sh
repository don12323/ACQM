#!/bin/bash

make clean

make

if [ $? -eq 0 ]; then
	./main
	gnuplot plotbasis.gp
	#gnuplot plotrad.gp
	#mv EN* Eplots/
else
	echo "Error" 
fi
