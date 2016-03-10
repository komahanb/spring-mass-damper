default:
	f95 -c functions.f90
	f95 -c abstract_runge_kutta.f90
	f95 -c explicit_runge_kutta.f90
	f95 -c implicit_runge_kutta.f90
	f95 -c compare_formulations.f90
	f95 -o test *.o -llapack

clean: 
	rm test *.o *.mod *~ 
