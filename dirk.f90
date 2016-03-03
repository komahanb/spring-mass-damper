!=====================================================================!
! Module that implements explicit, implicit and semi-implicit
! Runge-Kutta integration schemes
!=====================================================================!

module runge_kutta_class
  
  implicit none

  private

  public :: DIRK, IRK, ERK

  ! Abstract Runge-Kutta type
  type, abstract :: RK
     ! put all the common variables
   contains
     ! put all the common procedures
     ! procedure :: integrate
  end type RK

  ! Explicit Runge-Kutta
  type, extends(RK) :: ERK
  end type ERK

  ! Implicit Runge-Kutta  
  type, extends(RK) :: IRK
  end type IRK
  
  ! Diagonally implicit Runge-Kutta
  type, extends(RK) :: DIRK

     integer :: dirk_stages = 1 ! default number of stages
     real(8) :: h = 0.1d0       ! default step size ( will reconsider
                                ! when implementing adaptive step
                                ! size)
     real(8), dimension(:,:), allocatable :: A ! forms the coeff matrix
     real(8), dimension(:)  , allocatable :: B ! multiplies the state derivatives
     real(8), dimension(:)  , allocatable :: C ! multiplies the time

   contains

     procedure :: integrate
     procedure :: init, finalize
     procedure :: get_first_stage_deriv
     procedure :: get_approx_q
     procedure :: update_states
     procedure :: newton_solve

     !procedure :: get_second_stage_deriv
     !procedure :: approx_qdot

  end type DIRK

contains

  !--------------------------------------------------------------------!
  ! Procedure that will be called by the end user to march in time
  ! Input: 
  ! state arrays q and qdot with initial conditions set at q(1)
  ! number of steps N
  ! step size h
  !--------------------------------------------------------------------!

  subroutine integrate(this, q, qdot, N)
    
    class(dirk) :: this

    real(8), intent(inout), dimension(:) :: q, qdot
    real(8) :: time(10)
    integer, intent(in) :: N 
    integer :: k
    real(8) :: ydot(this % dirk_stages) ! stage derivatives
    
    ! March in time
    march: do k = 2, N + 1

       q(k) = q(k-1) 

       time(k) = time(k-1) + this % h
       
       ! find the stage derivatives ydot
       qdot  = q + h*c(1)*(ydot+small)

       call this % newton_solve(k, time(k), q(k), ydot(k))
       
       call this % update_states(k, q, ydot)

    end do march
    
  end subroutine integrate

  ! Newton solve to solve the linear system to get the stage derivatives at each time step
  subroutine newton_solve(this, k, time, q, qdot, ydot)
    
    class(dirk) :: this
    integer, intent(in) :: k ! current time step    
    real(8), intent(in) :: time
    real(8), intent(inout) :: q, qdot, ydot
    real(8), dimension(this % dirk_stages) :: R ! residuals
    real(8), dimension(this % dirk_stages, this % dirk_stages) :: dR !jacobians
    integer :: n, jj
    integer :: max_newton = 20
    real(8) :: tmp(this % dirk_stages)
    
    newton: do n = 1, max_newton
       
       ! Make as many residual and jacobian calls as the number of stages      
       forall(jj = 1:this % dirk_stages) &
            & R(jj) = residual(time + this % h * C(jj), &
            & q + this % h * this % B(jj) * ydot(jj), &
            & qdot) !qdot or ydot
       
!!$       forall(jj = 1:this % dirk_stages) R(jj) = residual(time, q, qdot+small)
!!$       call residual(ydot+small, q + h*c(1)*(ydot+small), t+h*b(1), tmp)

       ! FD approximation of the residual
       do jj = 1, this % dirk_stages
          R(jj) = residual(time, q, qdot)
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

  subroutine update_states(this, k, q, ydot, qdot, yddot)

    class(dirk) :: this
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

  !-------------------------------------------------------------------!
  ! Get the stage derivatives by solving the nonlinear system using
  ! Newton's method
  !-------------------------------------------------------------------!

  function get_first_stage_deriv(this) result(ydot)
    
    class(dirk) :: this
    real(8) :: q, qdot(this % dirk_stages)
    integer :: r
    real(8) :: ydot(this % dirk_stages)

    ! solve the nonlinear system to get the stage derivatives
    ! call newton1(1, time(I), q(i), K(1,i), b, c)

  end function get_first_stage_deriv

  ! Approximate q based on the runge-kutta scheme
  real(8) function get_approx_q(this)
    
    class(dirk) :: this
    real(8) :: q, qdot(this % dirk_stages)
    integer :: r, i, j
    real(8) :: ydot(this % dirk_stages)
   
    ydot = this % get_first_stage_deriv()
    
    do i = 1, this % dirk_stages
       do  j = 1, this % dirk_stages
          q = q + this % h * this % A(j,i)*ydot(j)
       end do
    end do
    
  end function get_approx_q

  ! Initialize the dirk datatype and construct the tableau
  subroutine init(this, dirk_stages, h)

    class(dirk) :: this
    integer, OPTIONAL, intent(in) :: dirk_stages
    real(8), OPTIONAL, intent(in) :: h

    ! set the order of integration
    if (present(dirk_stages)) this % dirk_stages = dirk_stages
    
    ! set the user supplied initial step size
    if (present(h)) this % h = h 
    
    ! allocate space for the tableau
    allocate(this % A(this % dirk_stages, this % dirk_stages))
    allocate(this % B(this % dirk_stages))    
    allocate(this % C(this % dirk_stages))

    ! put the entries into the tableau
    if (this % dirk_stages .eq. 1) then 
       this % A(1,1) = 0.5d0
       this % B(1)   = 1.0d0
       this % C(1)   = 0.5d0
    else
       stop"yet to add entries into the tableau"
    end if
    
  end subroutine init

  ! Deallocate the tableau entries
  subroutine finalize(this)

    class(dirk) :: this

    if(allocated(this % A)) deallocate(this % A)
    if(allocated(this % B)) deallocate(this % B)
    if(allocated(this % C)) deallocate(this % C)

  end subroutine finalize

end module runge_kutta_class

program main

  use runge_kutta_class

  implicit none

  type(DIRK) :: dirk1, dirk2
  real(8) :: h = 1.0d-2

  !real :: a(3), b(3)
  !a = (/1.,2.,3./)
  !b = (/4.,5.,6./)
  !print *, sum(a*b)

  print *, "Beginning execution"
  
  ! Initialize DIRK instance and create the Butcher tableau
  call dirk1 % init()
  
  ! call dirk % integrate()

  ! Deallocate the Butcher tableau
  call dirk1 % finalize()

  print *, "Executing complete"

end program main
