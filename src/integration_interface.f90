!=====================================================================!
! An interface module for all the integration schemes. The idea is to
! define a common interface that shall allow the endpoints to easily
! swap integration schmes at their discretion.
!=====================================================================!

module integration_interface
  
  use precision
  use abstract_nonlinear_solve
  
  implicit none
  
  ! Variables instances with getters and setters
  type, abstract :: integrator_descriptor

     ! Any implementation of the steady state solver: Newton solver, Newton-Krylov 
     class(abstract_nonlinear_solve), pointer :: nonlinear_solver => null()
     
     contains
       
       procedure :: set_nonlinear_solver

  end type integrator_descriptor

  type, extends(integrator_descriptor) :: abstract_integrator

  end type abstract_integrator

contains

  !-------------------------------------------------------------------!
  ! The endpoints call this function to set the preferred nonlinear
  ! solution method
  !-------------------------------------------------------------------!
  
  pure function set_nonlinear_solver(this, solver)
    
    class(integrator_descriptor) :: this
    class(abstract_nonlinear_solve), target :: solver
    
    this % nonlinear_solver => target_solver
    
  end function set_nonlinear_solver

end module integration_interface
