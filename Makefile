all:
	gfortran -O3 -ffree-line-length-none -c dfl.f90
	gfortran -O3 -ffree-line-length-none ode.f90 dfl.o -o ode
	gfortran -O3 -ffree-line-length-none example.f90 dfl.o -o example
clean:
	rm -f example
	rm -f ode
	rm -f *.o
	rm -f *.mod
