!=====================================================================!
! Module that manages the precision of real variables
!=====================================================================!

module precision

  !-------------------------------------------------------------------!
  !  Define constants to manage precision [TUNABLE]
  !-------------------------------------------------------------------!

  integer, parameter :: sp = kind(0.0)    ! single precision
  integer, parameter :: dp = kind(0.0d0)  ! double precision

end module precision


!=====================================================================!
! Module that stores the state of the function
!=====================================================================!

module variables

  use precision

  implicit none

  integer :: ndim   = 0
  integer :: nsteps = 1

  !-------------------------------------------------------------------!
  ! The state variables
  !-------------------------------------------------------------------!

  real(dp), dimension(:,:)   , allocatable :: q
  real(dp), dimension(:,:)   , allocatable :: qdot
  real(dp), dimension(:,:)   , allocatable :: qddot

  real(dp), dimension(:,:)   , allocatable :: R 
  real(dp), dimension(:,:,:) , allocatable :: dR
  real(dp), dimension(:)     , allocatable :: dq

  !-------------------------------------------------------------------!
  ! The system parameters 
  !-------------------------------------------------------------------!

  real(dp) :: M = 1.0_dp
  real(dp) :: C = 0.02_dp
  real(dp) :: K = 5.0_dp

  !-------------------------------------------------------------------!
  !  Integration constants
  !-------------------------------------------------------------------!  

  integer, parameter :: max_bdf_order = 3 
  integer, parameter :: max_fo_terms = max_bdf_order + 1
  integer, parameter :: max_so_terms = 2*max_bdf_order + 1
  
  real(dp)  :: alpha(max_bdf_order, max_fo_terms) = 0.0_dp
  real(dp)  :: beta(max_bdf_order, max_so_terms) = 0.0_dp
  
  real(dp), parameter :: tinit  = 0.0_dp
  real(dp), parameter :: tfinal = 10.0_dp

  real(dp) :: dt, dt2

contains
  
  ! initialization tasks for the variables in simulation
  subroutine initialize(ndimin, nstepsin)

    integer, intent(in) :: ndimin, nstepsin

    write(*, *) "Initializing variables"

    ! set the number of dimensions or variables
    ndim = ndimin
    if (ndim .eq. 0) stop "Wrong dimension"

    ! set the number of time steps
    nsteps = nstepsin
    if (nsteps .le. 0) stop "Wrong number of time steps"

    if (allocated(q)) deallocate(q)
    allocate(q(nsteps, ndim))
    q = 0.0_dp

    if (allocated(qdot)) deallocate(qdot)
    allocate(qdot(nsteps, ndim))
    qdot = 0.0_dp

    if (allocated(qddot)) deallocate(qddot)
    allocate(qddot(nsteps, ndim))
    qddot = 0.0_dp

    if (allocated(R)) deallocate(R)
    allocate(R(nsteps, ndim))
    R = 0.0_dp

    if (allocated(dq)) deallocate(dq)
    allocate(dq(ndim))
    dq = 0.0_dp

    if (allocated(dR)) deallocate(dR)
    allocate(dR(nsteps, ndim, ndim))
    dR = 0.0_dp
    
    !-------------------------------------------------------------------!
    ! Setup BDF coefficients
    !-------------------------------------------------------------------!
    
    ! set the BDF coefficeints for first derivative
    alpha(1, 1:2) = (/ 1.0, -1.0 /)
    alpha(2, 1:3) = (/ 1.5_dp, -2.0_dp, 0.5_dp /)
    alpha(3, 1:4) = (/ 11.0_dp/6.0_dp, -3.0_dp, 1.5_dp, -1.0_dp/3.0_dp /)

    ! set the BDF coefficient for second derivative
    beta(1, 1:3) = (/ 1.0_dp, -2.0_dp, 1.0_dp /)
    beta(2, 1:5) = (/ 2.25_dp, -6.0_dp, 5.5_dp, -2.0_dp, 0.25_dp /)

    dt  = (tfinal-tinit)/dble(nsteps)
    dt2 =  dt*dt

    ! initial condition
    q(1, :) = 1.0_dp
    
  end subroutine initialize

  ! finalization tasks
  subroutine finalize()

    write(*, *) "Finalize variables"

    if (allocated(q))     deallocate(q)
    if (allocated(qdot))  deallocate(qdot)
    if (allocated(qddot)) deallocate(qddot)
    if (allocated(R))     deallocate(R)
    if (allocated(dR))    deallocate(dR)
    if (allocated(dq))    deallocate(dq)

  end subroutine finalize
  
  ! implementation for residual
  subroutine residual(step_num)
    integer  :: step_num
    R(step_num, :) = M*qddot(step_num,:) + C*qdot(step_num,:) + K*q(step_num,:)
  end subroutine residual
  
  ! implementation for the jacobian
  subroutine jacobian(step_num)
    integer :: step_num
    dR(step_num, :, :) = K + C/dt + M/dt2
  end subroutine jacobian
  
  ! update the states at each newton iteration
  subroutine update_states(step_num)
    
    integer :: step_num

    !-------------------------------------------------------------!
    ! Update the state
    !-------------------------------------------------------------!
 
    q(step_num,:) = q(step_num,:) + dq(:)
    
    !-------------------------------------------------------------!
    ! FD approximation to I derivative
    !-------------------------------------------------------------!

    if (step_num .eq. 2) then
       ! first order approximation to I derivative
       qdot(step_num,:) &
            &= alpha(1, 1)*q(step_num,:)/dt &
            &+ alpha(1, 2)*q(step_num-1,:)/dt
    else if (step_num .eq. 3) then
       ! second order approximation to I derivative
       qdot(step_num,:) &
            &= alpha(2, 1)*q(step_num,:)/dt &
            &+ alpha(2, 2)*q(step_num-1,:)/dt &
            &+ alpha(2, 3)*q(step_num-2,:)/dt
    else
       ! third order approximation to I derivative
       qdot(step_num,:) &
            &= alpha(3, 1)*q(step_num,:)/dt &
            &+ alpha(3, 2)*q(step_num-1,:)/dt &
            &+ alpha(3, 3)*q(step_num-2,:)/dt &
            &+ alpha(3, 4)*q(step_num-3,:)/dt
    end if

    !-------------------------------------------------------------!
    ! FD approximation to II derivative
    !-------------------------------------------------------------!
    
    if (step_num .le. 3) then
       ! first order approx to II derivative
       qddot(step_num,:) &
            &= alpha(1, 1)*qdot(step_num,:)/dt &
            &+ alpha(1, 2)*qdot(step_num-1,:)/dt
    else
       ! second order approx to II derivative
       qddot(step_num,:) &
            &= alpha(1, 1)*qdot(step_num,:)/dt &
            &+ alpha(1, 2)*qdot(step_num-1,:)/dt &
            &+ alpha(1, 3)*qdot(step_num-2,:)/dt
    end if

    print*, dq,  q(step_num, :), qdot(step_num, :),qddot(step_num, :)

  end subroutine update_states

end module variables

!=====================================================================!
! A module that wraps all the data used in Newton solve
!
! Author :  Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module newton_solve_bean_class

  use precision

  implicit none

  private
  public :: newton_solve_bean

  type newton_solve_bean

!     private

     !----------------------------------------------------------------!
     ! Basic setup variables
     !----------------------------------------------------------------!

     integer :: max_newton_iters = 50
     integer :: nvars = 1
     integer :: iter_num = 0

     !----------------------------------------------------------------!
     ! Stopping criteria
     !----------------------------------------------------------------!

     ! absolute tolerances
     real(dp) :: atol_rnrm = 1.0d-12
     real(dp) :: atol_unrm = 1.0d-12

     ! relative tolerances
     real(dp) :: rtol_rnrm = 1.0d-6
     real(dp) :: rtol_unrm = 1.0d-6

     ! actual norm values
     real(dp), dimension(:), allocatable :: rnrm, unrm

     !----------------------------------------------------------------!
     ! Print and solution write control
     !----------------------------------------------------------------!

     integer :: filenum              = 6 
     logical :: write_newton_details = .true.
     logical :: exit_on_failure      = .false. 
     
     logical :: converged = .false. 

   contains

     procedure :: get_max_newton_iters, set_max_newton_iters
     procedure :: get_num_vars, set_num_vars

     procedure :: get_atol_unrm, set_atol_unrm
     procedure :: get_atol_rnrm, set_atol_rnrm

     procedure :: get_rtol_unrm, set_rtol_unrm
     procedure :: get_rtol_rnrm, set_rtol_rnrm

     procedure :: set_file_number
     procedure :: set_write_newton_details
     procedure :: set_exit_on_failure

!     procedure :: set_init_x, set_init_xdot

  end type newton_solve_bean

contains

  !-------------------------------------------------------------------!
  ! Get max_newton_iters
  !-------------------------------------------------------------------!

  integer function get_max_newton_iters(this)

    class(newton_solve_bean) :: this

    get_max_newton_iters = this % max_newton_iters

  end function get_max_newton_iters

  !-------------------------------------------------------------------!
  ! Set max_newton_iters
  !-------------------------------------------------------------------!

  subroutine set_max_newton_iters(this, max_newton_iters)

    class(newton_solve_bean) :: this
    integer :: max_newton_iters

    this % max_newton_iters =  max_newton_iters

  end subroutine set_max_newton_iters

  !-------------------------------------------------------------------!
  ! Get number of variables
  !-------------------------------------------------------------------!

  integer function get_num_vars(this)

    class(newton_solve_bean) :: this

    get_num_vars = this % nvars

  end function get_num_vars

  !-------------------------------------------------------------------!
  ! Set number of variables
  !-------------------------------------------------------------------!

  subroutine set_num_vars(this, nvars)

    class(newton_solve_bean) :: this
    integer :: nvars

    this % nvars = nvars

  end subroutine set_num_vars

  !-------------------------------------------------------------------!
  ! Get absolute tolerance of update norm
  !-------------------------------------------------------------------!

  real(dp) function get_atol_unrm(this)

    class(newton_solve_bean) :: this

    get_atol_unrm  = this % atol_unrm

  end function get_atol_unrm

  !-------------------------------------------------------------------!
  ! Set absolute tolerance of update norm
  !-------------------------------------------------------------------!

  subroutine set_atol_unrm(this, atol_unrm)

    class(newton_solve_bean) :: this
    real(dp) :: atol_unrm

    this % atol_unrm = atol_unrm

  end subroutine set_atol_unrm

  !-------------------------------------------------------------------!
  ! Get absolute tolerance of residual norm
  !-------------------------------------------------------------------!

  real(dp) function get_atol_rnrm(this)

    class(newton_solve_bean) :: this

    get_atol_rnrm  = this % atol_rnrm

  end function get_atol_rnrm

  !-------------------------------------------------------------------!
  ! Set absolute tolerance of residual norm
  !-------------------------------------------------------------------!

  subroutine set_atol_rnrm(this, atol_rnrm)

    class(newton_solve_bean) :: this
    real(dp) :: atol_rnrm

    this % atol_rnrm = atol_rnrm

  end subroutine set_atol_rnrm

  !-------------------------------------------------------------------!
  ! Get relative tolerance of update norm
  !-------------------------------------------------------------------!

  real(dp) function get_rtol_unrm(this)

    class(newton_solve_bean) :: this

    get_rtol_unrm  = this % rtol_unrm

  end function get_rtol_unrm

  !-------------------------------------------------------------------!
  ! Set relative tolerance of update norm
  !-------------------------------------------------------------------!

  subroutine set_rtol_unrm(this, rtol_unrm)

    class(newton_solve_bean) :: this
    real(dp) :: rtol_unrm

    this % rtol_unrm = rtol_unrm

  end subroutine set_rtol_unrm

  !-------------------------------------------------------------------!
  ! Get relative tolerance of residual norm
  !-------------------------------------------------------------------!

  real(dp) function get_rtol_rnrm(this)

    class(newton_solve_bean) :: this

    get_rtol_rnrm  = this % rtol_rnrm

  end function get_rtol_rnrm

  !-------------------------------------------------------------------!
  ! Set relative tolerance of residual norm
  !-------------------------------------------------------------------!

  subroutine set_rtol_rnrm(this, rtol_rnrm)

    class(newton_solve_bean) :: this
    real(dp) :: rtol_rnrm

    this % rtol_rnrm = rtol_rnrm

  end subroutine set_rtol_rnrm

  !-------------------------------------------------------------------!
  ! Set the output file number
  !-------------------------------------------------------------------!

  subroutine set_file_number(this, filenum)

    class(newton_solve_bean) :: this
    integer :: filenum

    this % filenum = filenum

  end subroutine set_file_number

  !-------------------------------------------------------------------!
  ! Set the print control for writing the details of newton solve
  !-------------------------------------------------------------------!

  subroutine set_write_newton_details(this, write_newton_details)

    class(newton_solve_bean) :: this
    logical :: write_newton_details

    this % write_newton_details = write_newton_details

  end subroutine set_write_newton_details

  !-------------------------------------------------------------------!
  ! Set whether or not to exit when not converged
  !-------------------------------------------------------------------!

  subroutine set_exit_on_failure(this, exit_on_failure)

    class(newton_solve_bean) :: this
    logical :: exit_on_failure

    this % exit_on_failure = exit_on_failure

  end subroutine set_exit_on_failure

  !-------------------------------------------------------------------!
  ! Set the initial value of x (starting point)
  !-------------------------------------------------------------------!

!!$  subroutine set_init_x(this, init_q)
!!$
!!$    class(newton_solve_bean) :: this
!!$    real(dp) :: init_q
!!$
!!$    this % init_q = init_q
!!$
!!$  end subroutine set_init_x
!!$
!!$  !-------------------------------------------------------------------!
!!$  ! Set the initial value of x (starting point)
!!$  !-------------------------------------------------------------------!
!!$
!!$  subroutine set_init_xdot(this, init_qdot)
!!$
!!$    class(newton_solve_bean) :: this
!!$    real(dp) :: init_qdot
!!$
!!$    this % init_qdot = init_qdot
!!$
!!$  end subroutine set_init_xdot

end module newton_solve_bean_class

!=====================================================================!
! A module that uses Newton's root finding method to approximate the
! solution to equations of the form R(q,qdot,qddot) = 0
!
! The module is generic, which means it can be used to find solutions
! (a) R = 0 (using R and dR values)
! (b) R'= 0 (using R' and d2R values)
!
! Usage : In a main program,
!
! type(newton_solve) :: newton ! create an instance of the solve
!
! call newton%init()
! call newton%solve()
!
! Author :  Komahan Boopathy (komahan@gatech.edu)
! =====================================================================!

module newton_solve_class

  ! import dependencies
  use precision
  use newton_solve_bean_class
  use variables

  ! No implicit varaible definitions
  implicit none

  ! All routines and varaibles are private by default
  private

  ! Expose datatypes
  public :: newton_solve

  ! A type that contains the logic for Newton's method
  type, extends(newton_solve_bean) :: newton_solve

contains

  ! private procedures
  procedure, private :: init  => init
  procedure, private :: finish => finish
  procedure, private :: work => work
  procedure, private :: linear_solve
  procedure, private :: check_stop
  procedure, private :: extrapolate

  ! public procedures
  procedure, public :: solve => solve

end type newton_solve

!!$! user interface for implementing the residual
!!$interface 
!!$   function get_residual(q, qdot, qddot) result(func_val)
!!$     use precision
!!$     implicit none
!!$     real(dp), optional :: q
!!$     real(dp), optional :: qdot
!!$     real(dp), optional :: qddot
!!$     real(dp) :: func_val
!!$   end function get_residual
!!$end interface

contains

  ! solve the linear system to find the newton update
  subroutine linear_solve(this, step_num)
    class(newton_solve) :: this
    integer :: step_num
    real(dp):: rtmp, drtmp
    
    if (ndim .eq. 1) then
       rtmp  = R(step_num, 1)
       drtmp = dR(step_num, 1, 1)
       dq    = -rtmp/drtmp
       !       print*, rtmp, drtmp, dq
    else
       stop"linear solve not implemented for ndim>1"
    end if

  end subroutine linear_solve

  ! solve the linear system to find the newton update
  subroutine check_stop(this, step_num)
    class(newton_solve) :: this
    integer :: step_num

    this % rnrm (this % iter_num) = norm2(R(step_num,:))
    this % unrm (this % iter_num) = norm2(dq)

    !    print*, this % iter_num, this % rnrm (this % iter_num), this % unrm (this % iter_num)

    if ((this % rnrm (this % iter_num) .le. this % get_atol_rnrm() ) .or. &
         &(this % unrm(this % iter_num) .le. this % get_atol_unrm() )) &
         & this % converged = .true.

  end subroutine check_stop

  ! extrapolate to next time step
  subroutine extrapolate(this, step_num)
    class(newton_solve) :: this
    integer :: step_num
    
    if (step_num .gt. 1) then
       q(step_num,:) = q(step_num-1,:) + qdot(step_num-1,:)*dt  &
            &+ qddot(step_num-1,:)*dt2/2.0_dp
    end if
    
  end subroutine extrapolate


  !-------------------------------------------------------------------!
  ! Return the residual of the function
  !-------------------------------------------------------------------!
!!$
!!$function get_residual(q, qdot, qddot) result(func_val)
!!$
!!$ use precision
!!$ 
!!$ real(dp), optional :: q
!!$ real(dp), optional :: qdot
!!$ real(dp), optional :: qddot
!!$
!!$ real(dp) :: func_val
!!$
!!$ real(dp) :: M = 1.0_dp
!!$ real(dp) :: C = 0.02_dp
!!$ real(dp) :: K = 5.0_dp
!!$
!!$ func_val = M*qddot + C*qdot + K*q
!!$
!!$end function get_residual

!-------------------------------------------------------------------!
! Routine for initialization tasks
!-------------------------------------------------------------------!

subroutine init(this)

 class(newton_solve) :: this

 print *, "Initializing Newton Solve"
 
 if (allocated(this%rnrm)) deallocate(this%rnrm)
 allocate(this%rnrm(this%get_max_newton_iters()))
 this%rnrm = 0.0_dp

 if (allocated(this%unrm)) deallocate(this%unrm)
 allocate(this%unrm(this%get_max_newton_iters()))
 this%unrm = 0.0_dp

 ! initialize module variables
 call initialize(this%get_num_vars(), 2)
! call this % extrapolate(step_num)

end subroutine init

!-------------------------------------------------------------------!
! Routine that wraps the logic of newton solve                                                          
!-------------------------------------------------------------------!

subroutine work(this)

  class(newton_solve) :: this

  real(dp) :: R, dR
  integer  :: n

  print *, "Executing Newton solve"
  
  call this % extrapolate(2)

  newton: do n = 1, this % get_max_newton_iters()

     ! call this % set_iter_num(n)
     this % iter_num = this% iter_num + 1

     ! call this % extrapolate()

     call residual(2)

     call jacobian(2)

     call this % linear_solve(2)

     call update_states(2)

     call this % check_stop(2)
     
     if (this % converged) exit newton
     
  end do newton

end subroutine work


!-------------------------------------------------------------------!
! Routine that performs finishing tasks
!-------------------------------------------------------------------!

subroutine finish(this)

  class(newton_solve) :: this

  print *, "Finish Newton solve"
  
  if (allocated(this%rnrm)) deallocate(this%rnrm)
  if (allocated(this%unrm)) deallocate(this%unrm)
 
  call finalize

end subroutine finish

!-------------------------------------------------------------------!
! Routine that performs newton solve                                                                   
!-------------------------------------------------------------------!

subroutine solve(this)

  class(newton_solve) :: this

  ! perform initialization tasks
  call this % init()

  ! perform newton solve
  call this % work()

  ! perform finalization tasks
  call this % finish()

end subroutine solve

end module newton_solve_class

!=====================================================================!
! Program to test newton_solve_class module
!=====================================================================!

program test_newton_solve_class

use newton_solve_class
!use variables

implicit none  

type(newton_solve) :: newton

integer :: nvars  = 1
integer :: nsteps = 1000

!call newton % init()

! Set optional parameters
call newton % set_num_vars(nvars)
call newton % set_exit_on_failure(.true.)

! solve the problem
call newton % solve()

!call newton % finish()

!call finalize


contains

end program test_newton_solve_class
