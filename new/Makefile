default:
	f95 -c functions.f90
	f95 -c abstract_runge_kutta.f90
	f95 -c explicit_runge_kutta.f90
	f95 -c implicit_runge_kutta.f90
	f95 -c compare_formulations.f90
	f95 -c order_conv_study.f90	
	f95 -c test_runge_kutta.f90
	f95 -o compare compare_formulations.o functions.o abstract_runge_kutta.o explicit_runge_kutta.o implicit_runge_kutta.o -llapack
	f95 -o order order_conv_study.o functions.o abstract_runge_kutta.o explicit_runge_kutta.o implicit_runge_kutta.o -llapack
	f95 -o test test_runge_kutta.o functions.o abstract_runge_kutta.o explicit_runge_kutta.o implicit_runge_kutta.o -llapack
clean: 
	rm test order compare *.o *.mod *~ 
