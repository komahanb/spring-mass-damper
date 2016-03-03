!=====================================================================!
! Module that implements explicit, implicit and semi-implicit
! Runge-Kutta integration schemes
!=====================================================================!

module runge_kutta_class
  
  implicit none

  private

  public :: DIRK, IRK, ERK

  !-------------------------------------------------------------------!
  ! Abstract Runge-Kutta type
  !-------------------------------------------------------------------!

  type, abstract :: RK
     
     integer :: num_stages = 1  ! default number of stages
     real(8) :: h = 0.1d0       ! default step size ( will reconsider
                                ! when implementing adaptive step
                                ! size)
     integer :: order           ! order of accuracy
     real(8), dimension(:,:), allocatable :: A ! forms the coeff matrix
     real(8), dimension(:)  , allocatable :: B ! multiplies the state derivatives
     real(8), dimension(:)  , allocatable :: C ! multiplies the time

   contains

     ! put all the common procedures
     procedure(integrate_interface), deferred :: integrate
     procedure(buthcher_interface), deferred  :: setup_butcher_tableau
     procedure :: initialize, finalize
     procedure :: update_states

  end type RK
  
  !-------------------------------------------------------------------!
  ! Interfaces for deferred specialized procedures 
  !-------------------------------------------------------------------!

  interface

     !----------------------------------------------------------------!
     ! Interface for integration logic
     !----------------------------------------------------------------!

     subroutine integrate_interface(this, q, qdot, N)
       import RK
       class(RK) :: this
       real(8), intent(inout), dimension(:) :: q, qdot
       integer, intent(in) :: N 
     end subroutine integrate_interface

     !----------------------------------------------------------------!
     ! Interface for setting the Butcher tableau for each type of RK
     ! scheme
     !----------------------------------------------------------------!

     subroutine buthcher_interface(this)
       import RK
       class(RK) :: this
     end subroutine buthcher_interface
     
  end interface

  !-------------------------------------------------------------------!  
  ! Explicit Runge-Kutta
  !-------------------------------------------------------------------!  

  type, extends(RK) :: ERK

   contains

     procedure :: integrate => integrateERK
     procedure :: setup_butcher_tableau => ButcherERK

  end type ERK

  !-------------------------------------------------------------------!  
  ! Diagonally implicit Runge-Kutta
  !-------------------------------------------------------------------!  

  type, extends(RK) :: DIRK

   contains

     procedure :: integrate => integrateDIRK
     procedure :: setup_butcher_tableau => ButcherDIRK

     procedure :: get_first_stage_deriv
     procedure :: get_approx_q
     procedure :: newton_solve

  end type DIRK

  !-------------------------------------------------------------------!  
  ! Implicit Runge-Kutta  
  !-------------------------------------------------------------------!  

  type, extends(DIRK) :: IRK

   contains

     procedure :: integrate => integrateIRK
     procedure :: setup_butcher_tableau => ButcherIRK

  end type IRK
  
contains

  !-------------------------------------------------------------------!
  ! Butcher's tableau for ERK 
  !-------------------------------------------------------------------!
  
  subroutine ButcherERK(this)

    class(ERK) :: this
    real(8), parameter :: half = 1.0d0/2.0d0
    real(8), parameter :: onethird = 1.0d0/3.0d0
    real(8), parameter :: onesixth = 1.0d0/6.0d0
    real(8), parameter :: oneeight = 1.0d0/8.0d0
    real(8), parameter :: one  = 1.0d0
    real(8), parameter :: alpha = half

    ! put the entries into the tableau
    if (this % num_stages .eq. 1) then 
       
       ! Explicit Euler
       
       this % B(1) = one
       this % C(1) = 0.0d0

       this % order = 1

    else if (this % num_stages .eq. 2) then 

       ! Midpoint method alpha = 1.0d0/2.0d0
       ! Raltson method alpha = 2.0d0/3.0d0
       ! Heun's Method alpha = 1.0d0
       
       this % A(2,1) = alpha

       this % B(1) = (one-half/alpha)
       this % B(2) = half/alpha

       this % C(1) = 0.0d0
       this % C(2) = alpha

       this % order = 2 

    else if (this % num_stages .eq. 3) then 

       ! Kutta's third order method
       
       this % A(2,1) = half
       this % A(3,1) = -one
       this % A(3,2)  = 2.0d0
       
       this % B(1) = onesixth
       this % B(2) = 2.0d0*onethird
       this % B(3) = onesixth

       this % C(1) = 0.0d0
       this % C(2) = half
       this % C(3) = one

       this % order = 3

    else if (this % num_stages .eq. 4) then 
       
       ! Kutta's formula more accurate (not the classical Runge formula)
       
       this % A(2,1) = onethird
       this % A(3,1) = -onethird
       this % A(4,1) = one
       this % A(3,2) = one
       this % A(4,2) = -one       
       this % A(4,3) = one

       this % B(1) = oneeight
       this % B(2) = 3.0d0*oneeight
       this % B(3) = this % B(2)
       this % B(4) = oneeight

       this % C(1) = 0.0d0
       this % C(2) = onethird
       this % C(3) = 2.0d0 * this % C(2)
       this % C(4) = one

       this % order = 4
       
!!$       ! The classical RK tableau
!!$       this % A(2,1) =  half
!!$       this % A(3,2) =  half
!!$       this % A(4,3) =  one
!!$
!!$       this % B(1) = onesixth
!!$       this % B(2) = onethird
!!$       this % B(3) = onethird
!!$       this % B(4) = onesixth
!!$
!!$       this % C(1) = 0.0d0
!!$       this % C(2) = half
!!$       this % C(3) = half
!!$       this % C(4) = one

    else
       
       print *, this % num_stages
       stop "ERK Butcher tableau is not implemented for the requested order"

    end if

  end subroutine ButcherERK
 
  !-------------------------------------------------------------------!
  ! Butcher's tableau for DIRK 
  !-------------------------------------------------------------------!

  subroutine ButcherDIRK(this)

    class(DIRK) :: this
    real(8), parameter :: PI = 22.0d0/7.0d0
    real(8), parameter :: tmp  = 1.0d0/(2.0d0*sqrt(3.0d0))
    real(8), parameter :: half = 1.0d0/2.0d0
    real(8), parameter :: one  = 1.0d0
    real(8), parameter :: alpha = 2.0d0*cos(PI/18.0d0)/sqrt(3.0d0)

    ! put the entries into the tableau (ROGER ALEXANDER 1977)
    if (this % num_stages .eq. 1) then 

       ! Implicit mid-point rule (A-stable)

       this % A(1,1) = half
       this % B(1)   = one
       this % C(1)   = half
       
       this % order = 2

!!$       ! Implicit Euler (Backward Euler) but first order accurate
!!$       this % A(1,1) = one
!!$       this % B(1)   = one
!!$       this % C(1)   = one
!!$       this % order = 1
       

    else if (this % num_stages .eq. 2) then

       ! Crouzeix formula (A-stable)

       this % A(1,1) = half + tmp
       this % A(2,1) = -one/sqrt(3.0d0)
       this % A(2,2) = this % A(1,1)

       this % B(1)   = half
       this % B(2)   = half

       this % C(1)   = half + tmp
       this % C(2)   = half - tmp

       this % order = 3

    else if (this % num_stages .eq. 3) then

       ! Crouzeix formula (A-stable)

       this % A(1,1) = (one+alpha)*half
       this % A(2,1) = -half*alpha
       this % A(3,1) =  one + alpha
       this % A(2,2) = this % A(1,1)
       this % A(3,2) = -(one + 2.0d0*alpha)
       this % A(3,3) = this % A(1,1)
       
       this % B(1)   = one/(6.0d0*alpha**2)
       this % B(2)   = one - one/(3.0d0*alpha**2)
       this % B(3)   = this % B(1)

       this % C(1) = (one + alpha)*half
       this % C(2) = half
       this % C(3) = (one - alpha)*half
       
       this % order = 4

    else if (this % num_stages .eq. 4) then

       stop "Four stage DIRK formula does not exist"
       
    else
       
       print *, this % num_stages
       stop "DIRK Butcher tableau is not implemented for the requested order"

    end if

  end subroutine ButcherDIRK

  !-------------------------------------------------------------------!
  ! Butcher's tableau for IRK 
  !-------------------------------------------------------------------!

  subroutine ButcherIRK(this)

    class(IRK) :: this
    
    ! put the entries into the tableau
    if (this % num_stages .eq. 1) then 

       ! Implicit mid-point rule (A-stable)

       this % A(1,1) = 0.5d0
       this % B(1)   = 1.0d0
       this % C(1)   = 0.5d0

       this % order = 2

!!$       ! Implicit Euler (Backward Euler) but first order accurate
!!$       this % A(1,1) = one
!!$       this % B(1)   = one
!!$       this % C(1)   = one
!!$       this % order = 1


    else if (this % num_stages .eq. 2) then 

       ! Radau II A scheme (2 step)

       this % A(1,1) = 5.0d0/12.0d0
       this % A(2,1) = 3.0d0/4.0d0
       this % A(1,2) = -1.0d0/12.0d0
       this % A(2,2) = 1.0d0/4.0d0

       this % B(1) = 3.0d0/4.0d0
       this % B(2) = 1.0d0/6.0d0

       this % C(1) = 1.0d0/3.0d0
       this % C(2) = 1.0d0

       this % order = 3

    else if (this % num_stages .eq. 3) then 

       ! Radau II A scheme (3 step)

       this % A(1,1) = 11.0d0/45.0d0 - 7.0d0*sqrt(6.0d0)/360.0d0
       this % A(2,1) = 37.0d0/225.0d0 + 169.0d0*sqrt(6.0d0)/1800.0d0
       this % A(3,1) = 4.0d0/9.0d0 - sqrt(6.0d0)/36.0d0

       this % A(1,2) = 37.0d0/225.0d0 - 169.0d0*sqrt(6.0d0)/1800.0d0
       this % A(2,2) = 11.0d0/45.0d0 + 7.0d0*sqrt(6.0d0)/360.0d0
       this % A(3,2) = 4.0d0/9.0d0 + sqrt(6.0d0)/36.0d0

       this % A(1,3) = -2.0d0/225.0d0 + sqrt(6.0d0)/75.0d0
       this % A(2,3) = -2.0d0/225.0d0 - sqrt(6.0d0)/75.0d0 
       this % A(3,3) = 1.0d0/9.0d0

       this % B(1) = 4.0d0/9.0d0 - sqrt(6.0d0)/36.0d0
       this % B(2) = 4.0d0/9.0d0 + sqrt(6.0d0)/36.0d0
       this % B(3) = 1.0d0/9.0d0

       this % C(1) = 2.0d0/5.0d0 - sqrt(6.0d0)/10.0d0
       this % C(2) = 2.0d0/5.0d0 + sqrt(6.0d0)/10.0d0
       this % C(3) = 1.0d0

       this % order = 4

    else

       print *, this % num_stages
       stop "IRK Butcher tableau is not implemented for the requested order"

    end if

  end subroutine ButcherIRK

  !-------------------------------------------------------------------!
  ! Time integration logic for ERK
  !-------------------------------------------------------------------!
  ! Input: 
  ! . state arrays q and qdot with initial conditions set at q(1)
  ! . number of steps N
  ! . step size h
  !-------------------------------------------------------------------!
  ! Output:
  ! . q, qdot arrays are modified by the routine
  !-------------------------------------------------------------------!
  
  subroutine IntegrateERK(this, q, qdot, N)

    class(ERK) :: this
    real(8), intent(inout), dimension(:) :: q, qdot
    integer, intent(in) :: N 

    real(8) :: time(10)
    real(8) :: ydot(this % num_stages) ! stage derivatives
    integer :: k

  end subroutine IntegrateERK

  !-------------------------------------------------------------------!
  ! Time integration logic for IRK
  !-------------------------------------------------------------------!
  ! Input: 
  ! . state arrays q and qdot with initial conditions set at q(1)
  ! . number of steps N
  ! . step size h
  !-------------------------------------------------------------------!
  ! Output:
  ! . q, qdot arrays are modified by the routine
  !-------------------------------------------------------------------!

  subroutine IntegrateIRK(this, q, qdot, N)

    class(IRK) :: this
    real(8), intent(inout), dimension(:) :: q, qdot
    integer, intent(in) :: N 

    real(8) :: time(10)
    real(8) :: ydot(this % num_stages) ! stage derivatives
    integer :: k

  end subroutine IntegrateIRK

  !-------------------------------------------------------------------!
  ! Time integration logic for DIRK
  !-------------------------------------------------------------------!
  ! Input: 
  ! . state arrays q and qdot with initial conditions set at q(1)
  ! . number of steps N
  ! . step size h
  !-------------------------------------------------------------------!
  ! Output:
  ! . q, qdot arrays are modified by the routine
  !-------------------------------------------------------------------!

  subroutine IntegrateDIRK(this, q, qdot, N)
    
    class(DIRK) :: this

    real(8), intent(inout), dimension(:) :: q, qdot
    integer, intent(in) :: N 
    real(8) :: ydot(this % num_stages) ! stage derivatives
    real(8) :: time(10)
    integer :: k

    ! Integration Logic
    
!!$    ! March in time
!!$    march: do k = 2, N + 1
!!$
!!$       q(k) = q(k-1) 
!!$
!!$       time(k) = time(k-1) + this % h
!!$       
!!$       ! find the stage derivatives ydot
!!$!       qdot  = q + h*c(1)*(ydot+small)
!!$
!!$ !      call this % newton_solve(k, time(k), q(k), ydot(k))
!!$       
!!$       call this % update_states(k, q, ydot)
!!$
!!$    end do march
   
  end subroutine IntegrateDIRK
  
  !-------------------------------------------------------------------!
  ! Newton solve to solve the linear system to get the stage
  ! derivatives at each time step
  !-------------------------------------------------------------------!

  subroutine newton_solve(this, k, time, q, qdot, ydot)
    
    class(dirk) :: this
    integer, intent(in) :: k ! current time step    
    real(8), intent(in) :: time
    real(8), intent(inout) :: q, qdot, ydot
    real(8), dimension(this % num_stages) :: R ! residuals
    real(8), dimension(this % num_stages, this % num_stages) :: dR !jacobians
    integer :: n, jj
    integer :: max_newton = 20
    real(8) :: tmp(this % num_stages)
    
    newton: do n = 1, max_newton
       
       ! Make as many residual and jacobian calls as the number of stages      
!!$       forall(jj = 1:this % num_stages) &
!!$            & R(jj) = residual(time + this % h * C(jj), &
!!$            & q + this % h * this % B(jj) * ydot(jj), &
!!$            & qdot) !qdot or ydot
       
!!$       forall(jj = 1:this % num_stages) R(jj) = residual(time, q, qdot+small)
!!$       call residual(ydot+small, q + h*c(1)*(ydot+small), t+h*b(1), tmp)

       ! FD approximation of the residual
       do jj = 1, this % num_stages
!          R(jj) = residual(time, q, qdot)
       end do

    end do newton

  end subroutine newton_solve
  
  !-------------------------------------------------------------------!
  ! Residual of the govenrning equations
  !-------------------------------------------------------------------!

  real(8) pure function residual(time, q, qdot)
    
    real(8), intent(in)  :: time
    real(8), intent(in)  :: q, qdot ! actual states

    residual = qdot + cos(q) - sin(time)
    
  end function residual

  !-------------------------------------------------------------------!
  ! Get the stage derivatives by solving the nonlinear system using
  ! Newton's method
  !-------------------------------------------------------------------!

  function get_first_stage_deriv(this) result(ydot)
    
    class(dirk) :: this
    real(8) :: q, qdot(this % num_stages)
    integer :: r
    real(8) :: ydot(this % num_stages)

    ! solve the nonlinear system to get the stage derivatives
    ! call newton1(1, time(I), q(i), K(1,i), b, c)

  end function get_first_stage_deriv
  
  !-------------------------------------------------------------------!
  ! Approximate q based on the runge-kutta scheme
  !-------------------------------------------------------------------!
  
  real(8) function get_approx_q(this)
    
    class(dirk) :: this
    real(8) :: q, qdot(this % num_stages)
    integer :: r, i, j
    real(8) :: ydot(this % num_stages)
   
    ydot = this % get_first_stage_deriv()
    
    do i = 1, this % num_stages
       do  j = 1, this % num_stages
          q = q + this % h * this % A(j,i)*ydot(j)
       end do
    end do
    
  end function get_approx_q

  !-------------------------------------------------------------------!
  ! Initialize the dirk datatype and construct the tableau
  !-------------------------------------------------------------------!

  subroutine initialize(this, num_stages, h)

    class(rk) :: this
    integer, OPTIONAL, intent(in) :: num_stages
    real(8), OPTIONAL, intent(in) :: h

    ! set the order of integration
    if (present(num_stages)) this % num_stages = num_stages
    
    ! set the user supplied initial step size
    if (present(h)) this % h = h 
    
    ! allocate space for the tableau
    allocate(this % A(this % num_stages, this % num_stages))
    this % A = 0.0d0

    allocate(this % B(this % num_stages))    
    this % B = 0.0d0

    allocate(this % C(this % num_stages))
    this % C = 0.0d0

    call this % setup_butcher_tableau()

  end subroutine initialize

  !-------------------------------------------------------------------!
  ! Deallocate the tableau entries
  !-------------------------------------------------------------------!

  subroutine finalize(this)

    class(rk) :: this

    if(allocated(this % A)) deallocate(this % A)
    if(allocated(this % B)) deallocate(this % B)
    if(allocated(this % C)) deallocate(this % C)

  end subroutine finalize

  !-------------------------------------------------------------------!
  ! Update the states based on RK Formulae
  !-------------------------------------------------------------------!
  
  subroutine update_states(this, k, q, ydot, qdot, yddot)

    class(RK) :: this
    integer, intent(in) :: k ! current time step
    real(8), intent(inout), dimension(:) :: q ! actual states
    real(8), intent(in), dimension(:) :: ydot ! first stage derivative
    real(8), OPTIONAL, intent(inout), dimension(:) :: qdot ! actual state
    real(8), OPTIONAL, intent(in), dimension(:) :: yddot ! second stage derivative

    ! update q (for first order ODE)
    q(k) = q(k-1) + this % h*sum(this % B(:)*ydot(:))

    ! update qdot (for second order ODE)
    if (present(qdot) .and. present(yddot))then
       qdot(k) = qdot(k-1) + this % h*sum(this % B(:)*yddot(:))
    end if

  end subroutine update_states

end module runge_kutta_class

program main

  use runge_kutta_class

  implicit none

  integer, parameter :: N=  100
  real(8), parameter :: h = 1.0d-2

  type(DIRK) :: DIRK1
  type(IRK)  :: IRK1
  type(ERK)  :: ERK1

  real(8) :: q(N+1) = 0.0d0, qdot(N+1) = 0.0d0

  print *, "> Beginning execution"

  ! set initial condition
  q(1) = 1.0

  print *, " > Explicit Runge Kutta"
  call ERK1  % initialize()
  call ERK1  % integrate(q, qdot, N)
  call ERK1  % finalize()
  
  print *, " > Implicit Runge Kutta"
  call IRK1  % initialize()
  call IRK1  % integrate(q, qdot, N)
  call IRK1  % finalize()

  print *, " > Diagonally-Implicit Runge Kutta"
  call DIRK1 % initialize()
  call DIRK1 % integrate(q, qdot, N)
  call DIRK1 % finalize()
  
  print *, "> Execution complete"

end program main










