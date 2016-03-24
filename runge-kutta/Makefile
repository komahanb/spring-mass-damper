default:
	gfortran -c -g -Wall functions.f90
	gfortran -c -g -Wall abstract_runge_kutta.f90
	gfortran -c -g -Wall explicit_runge_kutta.f90
	gfortran -c -g -Wall implicit_runge_kutta.f90
#	gfortran -c -g compare_formulations.f90
#	gfortran -c -g order_conv_study.f90	
	gfortran -c -g  -Wall test_runge_kutta.f90
#	gfortran -o compare compare_formulations.o functions.o abstract_runge_kutta.o explicit_runge_kutta.o implicit_runge_kutta.o -llapack
#	gfortran -o order order_conv_study.o functions.o abstract_runge_kutta.o explicit_runge_kutta.o implicit_runge_kutta.o -llapack
#	gfortran -o test test_runge_kutta.o functions.o abstract_runge_kutta.o explicit_runge_kutta.o  -llapack
	gfortran -o test test_runge_kutta.o functions.o abstract_runge_kutta.o explicit_runge_kutta.o implicit_runge_kutta.o -llapack
clean: 
	rm test order compare *.o *.mod *~ 
