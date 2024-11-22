
FC := gfortran

#Or switch to ifort if you have it
#FC := ifort 

#UNCOMMENT THE BELOW IF YOU IMPLEMENT OPENMP
#CFLAGS := -fopenmp   #(gfortran)
#CFLAGS := -qopenmp   #(ifort)

#UNCOMMENT THE BELOW FOR DEBUGGING:
#CFLAGS := $(CFLAGS) -g -Og -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow   #(gfortran)
#CFLAGS := $(CFLAGS) -g -O0 -debug extended -check all -traceback -ftrapuv   	       #(ifort)

main : constants.o dgesv.o main.o kgrid_igor.o iqpackd.o plql.o
	$(FC) $(CFLAGS) main.o constants.o dgesv.o kgrid_igor.o iqpackd.o plql.o -o main

main.o : main.f90
	$(FC) $(CFLAGS) -c main.f90

constants.o : constants.f90
	$(FC) $(CFLAGS) -c constants.f90

kgrid_igor.o : kgrid_igor.f
	$(FC) $(CFLAGS) -c kgrid_igor.f

iqpackd.o : iqpackd.f
	$(FC) $(CFLAGS) -c iqpackd.f

plql.o : iqpackd.f
	$(FC) $(CFLAGS) -c plql.f

dgesv.o : dgesv.f
	$(FC) $(CFLAGS) -c dgesv.f

clean : 
	rm *.o main
