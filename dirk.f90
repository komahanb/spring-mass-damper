!=====================================================================!
! Module that implements Explicit, Implicit and Semi-implicit
! Runge-Kutta integration schemes
!=====================================================================!
! o This module is suitable for the Explicit First order form
!                 qdot = f(q(t),t)

! o User needs to provide the implementation for the function 'f' in
!   the module
!=====================================================================!
! Author: Komahan Boopathy, komahan@gatech.edu
!=====================================================================!
  
module runge_kutta_class
  
  implicit none

  private

  public :: DIRK, IRK, ERK, F

  !-------------------------------------------------------------------!
  ! Abstract Runge-Kutta type
  !-------------------------------------------------------------------!

  type, abstract :: RK
     
     integer :: num_stages = 1  ! default number of stages
     integer :: order           ! order of accuracy
     real(8) :: h = 0.1d0       ! default step size ( will reconsider
                                ! when implementing adaptive step
                                ! size)

     ! The Butcher Tableau 
     real(8), dimension(:,:), allocatable :: A ! forms the coeff matrix
     real(8), dimension(:)  , allocatable :: B ! multiplies the state derivatives
     real(8), dimension(:)  , allocatable :: C ! multiplies the time

     ! The stage derivatives
     real(8), dimension(:)  , allocatable :: K ! the stage derivatives

   contains

     ! Deferred common procedures
     procedure(integrate_interface), deferred :: integrate
     procedure(stage_derivative_interface), deferred :: get_stage_derivatives
     procedure(buthcher_interface), deferred  :: setup_butcher_tableau

     ! Implemented common procedures
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
     ! Interface for finding the stage derivatives at each time step
     !----------------------------------------------------------------!
     
     subroutine stage_derivative_interface(this, k, q, qdot)
       import RK
       class(RK) :: this
       integer, intent(in) :: k 
       real(8), intent(in), dimension(:) :: q, qdot
     end subroutine stage_derivative_interface

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
     procedure :: get_stage_derivatives =>get_stage_derivativesERK
     procedure :: setup_butcher_tableau => ButcherERK

  end type ERK

  !-------------------------------------------------------------------!  
  ! Diagonally implicit Runge-Kutta
  !-------------------------------------------------------------------!  

  type, extends(RK) :: DIRK

   contains

     procedure :: integrate => integrateDIRK
     procedure :: setup_butcher_tableau => ButcherDIRK
     procedure :: get_stage_derivatives =>get_stage_derivativesDIRK

     procedure :: newton_solve

  end type DIRK

  !-------------------------------------------------------------------!  
  ! Implicit Runge-Kutta  
  !-------------------------------------------------------------------!  

  type, extends(DIRK) :: IRK

   contains

     procedure :: integrate => integrateIRK
     procedure :: setup_butcher_tableau => ButcherIRK
     procedure :: get_stage_derivatives =>get_stage_derivativesIRK

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
    real(8), parameter :: one = 1.0d0
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
       stop "IRK Butcher tableau is not implemented for the requested order/stage"

    end if

  end subroutine ButcherIRK

  !-------------------------------------------------------------------!
  ! Time integration logic for ERK
  !-------------------------------------------------------------------!
  ! Input: 
  ! o state arrays q and qdot with initial conditions set at q(1)
  ! o number of steps N
  ! o step size h
  !-------------------------------------------------------------------!
  ! Output:
  ! o q, qdot arrays are modified by the routine
  !-------------------------------------------------------------------!
  
  subroutine IntegrateERK(this, q, qdot, N)

    class(ERK) :: this
    real(8), intent(inout), dimension(:) :: q, qdot
    integer, intent(in) :: N 
    integer :: k

    ! March in time
    march: do k = 2, N + 1
       
       ! find the stage derivatives at the current step
       call this % get_stage_derivatives(k, q, qdot)
       
       ! advance the state to the current step
       call this % update_states(k, q, qdot)

    end do march

  end subroutine IntegrateERK

  !-------------------------------------------------------------------!
  ! Get the stage derivative array for the current step and states ERK
  !-------------------------------------------------------------------!
  
  subroutine get_stage_derivativesERK(this, k, q, qdot)

    class(ERK) :: this
    integer, intent(in) :: k 
    real(8), intent(in), dimension(:) :: q, qdot
    real(8) :: tau, Y
    integer :: j
    
    do j = 1, this % num_stages

       ! stage time
       tau = dble(k-2)*this % h + this % C(j)*this % h

       ! state q
       Y = q(k-1) + this % h*sum(this % A(j,:)*this % K(:))
       
       !qdot = f(t,q)
       this % K(j) =  F(tau, Y)

    end do
    
  end subroutine get_stage_derivativesERK

  !-------------------------------------------------------------------!
  ! Get the stage derivative array for the current step and states DIRK
  !-------------------------------------------------------------------!
  
  subroutine get_stage_derivativesDIRK(this, k, q, qdot)

    class(DIRK) :: this
    integer, intent(in) :: k 
    real(8), intent(in), dimension(:) :: q, qdot
    integer :: j

  end subroutine get_stage_derivativesDIRK


  !-------------------------------------------------------------------!
  ! Get the stage derivative array for the current step and states IRK
  !-------------------------------------------------------------------!
  
  subroutine get_stage_derivativesIRK(this, k, q, qdot)

    class(IRK) :: this
    integer, intent(in) :: k 
    real(8), intent(in), dimension(:) :: q, qdot
    integer :: j

  end subroutine get_stage_derivativesIRK

  !-------------------------------------------------------------------!
  ! Time integration logic for IRK
  !-------------------------------------------------------------------!
  ! Input: 
  ! o state arrays q and qdot with initial conditions set at q(1)
  ! o number of steps N
  ! o step size h
  !-------------------------------------------------------------------!
  ! Output:
  ! o q, qdot arrays are modified by the routine
  !-------------------------------------------------------------------!

  subroutine IntegrateIRK(this, q, qdot, N)

    class(IRK) :: this
    real(8), intent(inout), dimension(:) :: q, qdot
    integer, intent(in) :: N 
    integer :: k

  end subroutine IntegrateIRK

  !-------------------------------------------------------------------!
  ! Time integration logic for DIRK
  !-------------------------------------------------------------------!
  ! Input: 
  ! o state arrays q and qdot with initial conditions set at q(1)
  ! o number of steps N
  ! o step size h
  !-------------------------------------------------------------------!
  ! Output:
  ! o q, qdot arrays are modified by the routine
  !-------------------------------------------------------------------!

  subroutine IntegrateDIRK(this, q, qdot, N)
    
    class(DIRK) :: this

    real(8), intent(inout), dimension(:) :: q, qdot
    integer, intent(in) :: N 
    integer :: k
   
  end subroutine IntegrateDIRK
  
  !-------------------------------------------------------------------!
  ! Newton solve to solve the linear system to get the stage
  ! derivatives at each time step
  !-------------------------------------------------------------------!

  subroutine newton_solve(this, k, time, q, qdot)
    
    class(dirk) :: this
    integer, intent(in) :: k ! current time step    
    real(8), intent(in) :: time
    real(8), intent(inout) :: q, qdot
    integer :: max_newton = 20
    integer :: n, jj

    newton: do n = 1, max_newton
       
    end do newton

  end subroutine newton_solve
  
  !-------------------------------------------------------------------!
  ! F of the govenrning equations
  !-------------------------------------------------------------------!

  real(8) pure function F(time, q)
    
    real(8), intent(in)  :: time
    real(8), intent(in)  :: q

    F = - sin(time)
    
    ! F = qdot + cos(q) - sin(time)
    ! F = qdot - cos(time)
    ! F = cos(q) - sin(time)
    ! F = exp(time)

  end function F
  
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

    ! allocate space for the stage derivatives
    allocate(this % K(this % num_stages))
    this % K = 0.0d0
    
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
    if(allocated(this % K)) deallocate(this % K)

  end subroutine finalize

  !-------------------------------------------------------------------!
  ! Update the states based on RK Formulae
  !-------------------------------------------------------------------!
  
  subroutine update_states(this, k, q, qdot)

    class(RK) :: this
    integer, intent(in) :: k ! current time step
    real(8), intent(inout), dimension(:) :: q ! actual states
    real(8), intent(inout), dimension(:) :: qdot ! actual state

    ! update q (for first order ODE)
    q(k) = q(k-1) + this % h*sum(this % B(:)*this % K(:))
    
    ! update the qdot value
    qdot(k) = sum(this % B(:)*this % K(:))
    
  end subroutine update_states

end module runge_kutta_class

program main

  use runge_kutta_class

  implicit none

  integer, parameter :: N = 100
  real(8), parameter :: h = 1.0d-1

  type(DIRK) :: DIRKOBJ
  type(IRK)  :: IRKOBJ
  type(ERK)  :: ERKOBJ
  
  real(8) :: q(4, N+1) = 0.0d0, qdot(4, N+1) = 0.0d0, error(4, N+1) = 0.0d0
  integer :: i, kk
  
  !-------------------------------------------------------------------!
  ! Explicit Runge Kutta
  !-------------------------------------------------------------------!

  q = 0.0d0; q(:,1) = 1.0

  do kk = 1, 4

     ! Test for each RK stage
     call ERKOBJ  % initialize(kk)
     call ERKOBJ  % integrate(q(kk,:), qdot(kk,:), N)
     call ERKOBJ  % finalize()
     
     ! Find the error
     do i = 1, N +1
        error(kk,i) = abs(q(kk,i) - cos(dble(i-1)*ERKOBJ % h))
     end do

  end do

  write (*, '(4E15.8)') (norm2(error(i,:)), i = 1,4)
  
  !-------------------------------------------------------------------!
  ! Implicit Runge Kutta
  !-------------------------------------------------------------------!
  
  q= 0.0d0; q(:,1) = 1.0
  
  do kk = 1, 3

     ! Test for each RK stage
     call IRKOBJ  % initialize(kk)
     call IRKOBJ  % integrate(q(kk,:), qdot(kk,:), N)
     call IRKOBJ  % finalize()
     
     ! Find the error
     do i = 1, N +1
        error(kk,i) = abs(q(kk,i) - cos(dble(i-1)*IRKOBJ % h))
     end do

  end do

  write (*, '(4E15.8)') (norm2(error(i,:)), i = 1, 3)
  
  !-------------------------------------------------------------------!
  ! Diagonally Implicit Runge Kutta
  !-------------------------------------------------------------------!
  
  q= 0.0d0; q(:,1) = 1.0

  do kk = 1, 3

     ! Test for each RK stage
     call DIRKOBJ  % initialize(kk)
     call DIRKOBJ  % integrate(q(kk,:), qdot(kk,:), N)
     call DIRKOBJ  % finalize()

     ! Find the error
     do i = 1, N +1
        error(kk,i) = abs(q(kk,i) - cos(dble(i-1)*DIRKOBJ % h))
     end do

  end do

  write (*, '(4E15.8)') (norm2(error(i,:)), i = 1, 3)

end program main

!!$ 
!!$  write (*, '(4E15.8)', advance='yes') (dble(i-1)*ERKOBJ % h, &
!!$       & (q(i)), abs(q(i)-cos(dble(i-1)*ERKOBJ % h)), &
!!$       & F(dble(i-1)*ERKOBJ % h, &
!!$       & q(i), qdot(i)), i = 1, N+1)
!!$  
  !print *, " > Implicit Runge Kutta"
!!$  call IRKOBJ  % initialize()
!!$  call IRKOBJ  % integrate(q, qdot, N)
!!$  call IRKOBJ  % finalize()
!!$
!!$  !print *, " > Diagonally-Implicit Runge Kutta"
!!$  call DIRKOBJ % initialize()
!!$  call DIRKOBJ % integrate(q, qdot, N)
!!$  call DIRKOBJ % finalize()
!!$  !-------------------------------------------------------------------!
!!$  ! Get the stage derivatives by solving the nonlinear system using
!!$  ! Newton's method
!!$  !-------------------------------------------------------------------!
!!$
!!$  function get_first_stage_deriv(this) result(ydot)
!!$    
!!$    class(dirk) :: this
!!$    real(8) :: q, qdot(this % num_stages)
!!$    integer :: r
!!$    real(8) :: ydot(this % num_stages)
!!$
!!$    ! solve the nonlinear system to get the stage derivatives
!!$    ! call newton1(1, time(I), q(i), K(1,i), b, c)
!!$
!!$  end function get_first_stage_deriv
