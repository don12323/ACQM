FC := gfortran

CFLAGS := -O3 -ffree-line-length-512 -fdefault-real-8 -fdefault-double-8

#UNCOMMENT THE BELOW LINE FOR DEBUGGING:
CFLAGS := $(CFLAGS) -g -Og -fcheck=all -fbacktrace
all: main main2
main : 
	$(FC) $(CFLAGS) main.f90 -o main
main2 : 
	$(FC) $(CFLAGS) rsg.f intp.f main2.f90 -o main2

clean : 
	rm -f main2 main
