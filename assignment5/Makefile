FC := gfortran
CFLAGS := -O3 -ffree-line-length-512 -fdefault-real-8 -fdefault-double-8
# Uncomment the below line for debugging:
# CFLAGS := $(CFLAGS) -g -Og -fcheck=all -fbacktrace

all: main main2

main: constants.o main1.o
	$(FC) $(CFLAGS) constants.o main1.o rsg.f -o main

main1.o: main.f90
	$(FC) $(CFLAGS) -c main.f90 -o main1.o

main2: constants.o dgesv.o main2.o kgrid_igor.o iqpackd.o plql.o
	$(FC) $(CFLAGS) main2.o constants.o dgesv.o kgrid_igor.o iqpackd.o plql.o -o main2

main2.o: main2.f90
	$(FC) $(CFLAGS) -c main2.f90 -o main2.o

constants.o: constants.f90
	$(FC) $(CFLAGS) -c constants.f90

dgesv.o: dgesv.f
	$(FC) $(CFLAGS) -c dgesv.f

kgrid_igor.o: kgrid_igor.f
	$(FC) $(CFLAGS) -c kgrid_igor.f

iqpackd.o: iqpackd.f
	$(FC) $(CFLAGS) -c iqpackd.f

plql.o: iqpackd.f
	$(FC) $(CFLAGS) -c plql.f

clean:
	rm -f main main2 *.o

