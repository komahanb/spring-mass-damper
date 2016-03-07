default:
	gfortran -c functions.f90
	gfortran -c abstract_runge_kutta.f90
	gfortran -c explicit_runge_kutta.f90
	gfortran -c implicit_runge_kutta.f90
	gfortran -c test_runge_kutta.f90
	gfortran -o test *.o -llapack

clean: 
	rm *.o *.mod *~
