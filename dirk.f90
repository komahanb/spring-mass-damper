!=====================================================================!
! Module that implements Explicit, Implicit and Semi-implicit
! Runge-Kutta integration schemes
!=====================================================================!
! o This module is suitable for the Explicit First order form
!                 qdot = f(q(t),t)

! o User needs to provide the implementation for the function 'f' in
!   the module
! 
! o Yet to implement the multivariate case
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

     ! The stage time and its corresponding derivatives
     real(8), dimension(:)  , allocatable :: T ! the corresponding stage time
     real(8), dimension(:)  , allocatable :: Y ! the corresponding state
     real(8), dimension(:)  , allocatable :: K ! the stage derivatives K = F(T,Y)

     real(8), dimension(:)  , allocatable :: R ! stage residual
     real(8), dimension(:,:), allocatable :: J ! stage jacobian
          
   contains

     ! Deferred common procedures
     procedure(compute_stage_values_interface), deferred :: compute_stage_values
     procedure(buthcher_interface), deferred  :: setup_butcher_tableau

     ! Implemented common procedures
     procedure :: initialize, finalize
     procedure :: update_states
     procedure :: integrate
     procedure :: reset_stage_values

  end type RK
  
  !-------------------------------------------------------------------!
  ! Interfaces for deferred specialized procedures 
  !-------------------------------------------------------------------!

  interface

     !----------------------------------------------------------------!
     ! Interface for finding the stage derivatives at each time step
     !----------------------------------------------------------------!
     
     subroutine compute_stage_values_interface(this, k, q)
       import RK
       class(RK) :: this
       integer, intent(in) :: k 
       real(8), intent(in), dimension(:) :: q
     end subroutine compute_stage_values_interface

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

     procedure :: compute_stage_values =>compute_stage_valuesERK
     procedure :: setup_butcher_tableau => ButcherERK

  end type ERK

  !-------------------------------------------------------------------!  
  ! Diagonally implicit Runge-Kutta
  !-------------------------------------------------------------------!  

  type, extends(RK) :: DIRK

   contains

     procedure :: setup_butcher_tableau => ButcherDIRK
     procedure :: compute_stage_values =>compute_stage_valuesDIRK

     procedure :: newton_solve

     procedure :: compute_stage_residual
     procedure :: compute_stage_jacobian

  end type DIRK

  !-------------------------------------------------------------------!  
  ! Implicit Runge-Kutta  
  !-------------------------------------------------------------------!  

  type, extends(DIRK) :: IRK

   contains

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
  ! Get the stage derivative array for the current step and states ERK
  !-------------------------------------------------------------------!
  
  subroutine compute_stage_valuesERK(this, k, q)

    class(ERK) :: this
    integer, intent(in) :: k 
    real(8), intent(in), dimension(:) :: q
    real(8) :: tmp
    integer :: j, i 
    
    ! Stage derivatives are explicitly found at each iteration

    do j = 1, this % num_stages

       ! stage time
       this % T(j) = dble(k-2)*this % h + this % C(j)*this % h

       ! stage Y
       !this % Y(j) = q(k-1) + this % h*sum(this % A(j,:)*this % K(:))
       tmp = 0.0d0
       do i = 1, j - 1
          tmp = tmp + this % A(j,i)*this % K(j)
       end do
       this % Y(j) = q(k-1) +  this % h * tmp

       ! stage derivative
       this % K(j) =  F(this % T(j), this % Y(j))

    end do
    
  end subroutine compute_stage_valuesERK

  !-------------------------------------------------------------------!
  ! Get the stage derivative array for the current step and states DIRK
  !-------------------------------------------------------------------!
  
  subroutine compute_stage_valuesDIRK(this, k, q)

    class(DIRK) :: this
    integer, intent(in) :: k 
    real(8), intent(in), dimension(:) :: q
    integer :: j

    ! Find the stage times
    do j = 1, this % num_stages
       this % T(j) = dble(k-2)*this % h + this % C(j)*this % h
    end do
    
    ! Guess the solution for stage states
    this % Y(:) = 1.0d0
    
    ! solve the non linear stage equations using Newton's method for
    ! the actual stage states 
    ! S function calls
    ! S jacobian calls
    ! Times the number of newton iterations (N)
    call this % newton_solve(q(k-1))

    ! at this point Y, K, T are finalized--just like ERK
    
  end subroutine compute_stage_valuesDIRK

  !-------------------------------------------------------------------!
  ! Solve nonlinear stage equations using Newton's method at each time
  ! step
  !-------------------------------------------------------------------!

  subroutine newton_solve(this, qk)
    
    class(dirk) :: this
    real(8), intent(in) :: qk
    real(8) :: dQ(this % num_stages)
    integer :: max_newton = 20
    integer :: n, jj
    integer :: ipiv(this%num_stages), info
    logical :: conv = .false.
    
    newton: do n = 1, max_newton
       
       ! Get the residual of the function
       call this % compute_stage_residual(qk)
       
       ! Get the jacobian matrix
       call this % compute_stage_jacobian()

       ! call lapack to solve the system
       dq = - this % R
       call DGESV(this % num_stages, 1, this % J, this % num_stages, &
            & IPIV, dq, this % num_stages, INFO)

       ! check stop (change this to norm2 when implementing
       ! multivariate case)
       if(norm2(this % R) .le. 1.0d-12 .or. norm2(dq) .le. 1.0d-12) then
          conv = .true.
          exit newton
       end if
       
       ! update q
       this % Y = this % Y + dQ
       
    end do newton

    if (.not.conv) print *, "Newton solve failed after", n

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
  ! DFDY of the function
  !-------------------------------------------------------------------!

  real(8) pure function DFDY(time, Y)
    
    real(8), intent(in)  :: time
    real(8), intent(in)  :: Y

    DFDY = 0.0d0
    
    ! F = qdot + cos(q) - sin(time)
    ! F = qdot - cos(time)
    ! F = cos(q) - sin(time)
    ! F = exp(time)
    
  end function DFDY

  !-------------------------------------------------------------------!
  ! Computes the stage residual and sets into the same instance
  !-------------------------------------------------------------------!
  
  subroutine compute_stage_residual(this, qk)

    class(DIRK) :: this
    real(8) :: qk
    integer :: i

    do i = 1, this % num_stages

       ! this could be the last newton iteration, so we store the
       ! function call value; may be there is better way to handle
       ! this logic
       
       this % K(i) = F(this % T(i), this % Y(i))

       this % R(i) = this % Y(i) - qk - this % h * this % K(i)
       
    end do
    
  end subroutine compute_stage_residual

  !-------------------------------------------------------------------!
  ! Computes the stage jacobian and sets into the same instance  
  !          J[s x s] = [I(s)-h A(i,i) DFDY(T(i),Y(i))]
  !-------------------------------------------------------------------!

  subroutine compute_stage_jacobian(this)

    class(DIRK) :: this
    integer :: i, j
    
    select type (this)

    type is (DIRK)

       ! Jacobian is a lower triangle matrix
       do i = 1, this % num_stages
          do j = 1, i
             this % J(j,i) = 1.0d0 - this % h * this % A(j,i) &
                  &* DFDY(this % T(i), this % Y(i))
          end do
       end do

    type is (IRK)

       ! Jacobian is a FULL  matrix
       do i = 1, this % num_stages
          do j = 1, this % num_stages
             this % J(j,i) = 1.0d0 - this % h * this % A(j,i) &
                  &* DFDY(this % T(i), this % Y(i))
          end do
       end do

    end select

    !this % J(i,i) = 1.0d0 - this % h * this % A(i,i)
    !* DFDY(this % T(i), this % Y(i))
    
  end subroutine compute_stage_jacobian

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

    ! allocate space for the stage state
    allocate(this % Y(this % num_stages))
    this % Y = 0.0d0

    ! allocate space for the stage time
    allocate(this % T(this % num_stages))
    this % T = 0.0d0

    ! allocate space for the stage time
    allocate(this % R(this % num_stages))
    this % R = 0.0d0

    ! allocate space for the stage time
    allocate(this % J(this % num_stages, this % num_stages))
    this % J = 0.0d0

    call this % setup_butcher_tableau()

  end subroutine initialize

  !-------------------------------------------------------------------!
  ! Deallocate the tableau entries
  !-------------------------------------------------------------------!

  subroutine finalize(this)

    class(rk) :: this

    ! clear butcher's tableau
    if(allocated(this % A)) deallocate(this % A)
    if(allocated(this % B)) deallocate(this % B)
    if(allocated(this % C)) deallocate(this % C)

    ! clear stage values
    if(allocated(this % K)) deallocate(this % K)
    if(allocated(this % T)) deallocate(this % T)
    if(allocated(this % Y)) deallocate(this % Y)

    ! clear the stage residual and jacobian
    if(allocated(this % R)) deallocate(this % R)
    if(allocated(this % J)) deallocate(this % J)
    
  end subroutine finalize

  !-------------------------------------------------------------------!
  ! Time integration logic
  !-------------------------------------------------------------------!
  ! Input: 
  ! o state arrays q and qdot with initial conditions set at q(1)
  ! o number of steps N
  ! o step size h
  !-------------------------------------------------------------------!
  ! Output:
  ! o q, qdot arrays are modified by the routine
  !-------------------------------------------------------------------!
  
  subroutine Integrate(this, q, qdot, N)

    class(RK) :: this
    real(8), intent(inout), dimension(:) :: q, qdot
    integer, intent(in) :: N 
    integer :: k

    ! March in time
    march: do k = 2, N + 1

       ! find the stage derivatives at the current step
       call this % compute_stage_values(k, q)

       ! advance the state to the current step
       call this % update_states(k, q, qdot)

       ! set the stage values to zero
       call this % reset_stage_values()

    end do march

  end subroutine Integrate

  !-------------------------------------------------------------------!
  ! Update the states based on RK Formulae
  !-------------------------------------------------------------------!
  
  subroutine update_states(this, k, q, qdot)

    class(RK) :: this
    integer, intent(in) :: k ! current time step
    real(8), intent(inout), dimension(:) :: q ! actual states
    real(8), intent(inout), dimension(:) :: qdot ! actual state

    ! update the qdot value
    qdot(k) = sum(this % B(:)*this % K(:))

    ! update q (for first order ODE)
    q(k) = q(k-1) + this % h * qdot(k)
    
  end subroutine update_states

  !-------------------------------------------------------------------!
  ! Reset the array to store new stage values at each time step
  !-------------------------------------------------------------------!
  
  subroutine reset_stage_values(this)

    class(RK) :: this

    this % K = 0.0d0
    this % Y = 0.0d0
    this % T = 0.0d0

    this % R = 0.0d0
    this % J = 0.0d0

  end subroutine reset_stage_values

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

  write (*, '(a,4E15.8)') "ERK :",(norm2(error(i,:)), i = 1,4)
  
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

  write (*, '(a,4E15.8)') "IRK :", (norm2(error(i,:)), i = 1, 3)

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

  write (*, '(a,4E15.8)') "DIRK:", (norm2(error(i,:)), i = 1, 3)

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
