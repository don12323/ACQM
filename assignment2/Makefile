FC := gfortran

CFLAGS := -O3 -ffree-line-length-512 -fdefault-real-8 -fdefault-double-8 

#UNCOMMENT THE BELOW LINE FOR DEBUGGING:
CFLAGS := $(CFLAGS) -g -Og -fcheck=all -fbacktrace

#Write your program in main.f90
all: main main2

main : main.f90 rsg.f wigner.f 
	$(FC) $(CFLAGS) main.f90 rsg.f wigner.f -o main

main2 : main2.f90 rsg.f intp.f 
	$(FC) $(CFLAGS) main2.f90 rsg.f intp.f -o main2

clean : 
	rm -f main main2 
