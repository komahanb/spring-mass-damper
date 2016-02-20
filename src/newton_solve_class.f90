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

module system_class

  use precision
  
  implicit none

  private

  public :: abstract_residual, abstract_jacobian
  public :: system, state, mesh, initialize_simulation, set_initial_condition
  
  type system_descriptor
     
     private 

     integer  :: nvars = 1
     integer  :: num_time_steps  = 1
     integer  :: current_time_step = 1

     real(dp) :: start_time = 0.0_dp
     real(dp) :: end_time = 1.0_dp
     real(dp) :: current_time = 0.0_dp
     real(dp) :: dt = 1.0e-3_dp

     logical  :: unsteady = .false.
     
   contains
     
      procedure :: initialize_system_variables
      
      procedure :: get_num_vars, set_num_vars
      procedure :: get_num_time_steps, set_num_time_steps
      procedure :: get_current_time_step, set_current_time_step
      procedure :: get_current_time, set_current_time

      procedure :: get_step_size, set_step_size
      procedure :: get_start_time, set_start_time
      procedure :: get_end_time, set_end_time

      procedure :: set_unsteady, is_unsteady

  end type system_descriptor

  !-------------------------------------------------------------------!
  ! The state variables
  !-------------------------------------------------------------------!

  type :: state_variables

     real(dp), dimension(:,:), allocatable :: q
     real(dp), dimension(:,:), allocatable :: qdot
     real(dp), dimension(:,:), allocatable :: qddot

     real(dp), dimension(:),   allocatable :: dq

     real(dp), dimension(:,:), allocatable :: R
     real(dp), dimension(:,:,:), allocatable :: dR

   contains

     procedure :: set_initial_condition
     procedure :: initialize_state_variables
     procedure :: finalize_state_variables

  end type state_variables
  
  !-------------------------------------------------------------------!
  ! The mesh variables
  !-------------------------------------------------------------------!
  
  type :: mesh_variables

   contains
     
     procedure :: initialize_mesh_variables
     procedure :: finalize_mesh_variables

  end type mesh_variables
  
  !-------------------------------------------------------------------!
  ! Abstract type for residual vector and its interface
  !-------------------------------------------------------------------!
  
  type, abstract :: abstract_residual

     
   contains

     procedure(get_residual_interface), deferred :: get_residual

  end type abstract_residual

  !-------------------------------------------------------------------!
  ! Interface for implementing the residual
  !-------------------------------------------------------------------!

  interface

     subroutine get_residual_interface(this)

       import abstract_residual

       class(abstract_residual) :: this

     end subroutine get_residual_interface

  end interface

  !-------------------------------------------------------------------!
  ! Abstract type for jacobian matrix and its interface
  !-------------------------------------------------------------------!

  type, abstract :: abstract_jacobian

   contains

     procedure(get_jacobian_interface), deferred :: get_jacobian
     
  end type abstract_jacobian

  !-------------------------------------------------------------------!
  ! Interface for implementing the jacobian
  !-------------------------------------------------------------------!

  interface

     subroutine get_jacobian_interface(this)

       import abstract_jacobian

       class(abstract_jacobian) :: this

     end subroutine get_jacobian_interface

  end interface

  !-------------------------------------------------------------------!
  ! The variables that are shared all over and contains main data
  !-------------------------------------------------------------------!

  type(system_descriptor) :: system
  type(state_variables)   :: state
  type(mesh_variables)    :: mesh
  
contains

  !-------------------------------------------------------------------!
  ! Initialization tasks for ALL  variables in the simulation
  !-------------------------------------------------------------------!
  
  subroutine initialize_simulation()

    call system % initialize_system_variables()
    call state  % initialize_state_variables()
    call mesh   % initialize_mesh_variables()

  end subroutine initialize_simulation

  !-------------------------------------------------------------------!
  ! Set the number of variables
  !-------------------------------------------------------------------!  

  subroutine set_num_vars(this, nvars)

    class(system_descriptor) :: this
    integer :: nvars

    this % nvars = nvars
    
  end subroutine set_num_vars

  !-------------------------------------------------------------------!
  ! Set the number of time steps
  !-------------------------------------------------------------------!  

  subroutine set_num_time_steps(this, nsteps)

    class(system_descriptor) :: this
    integer :: nsteps

    this % num_time_steps = nsteps
    
  end subroutine set_num_time_steps  

  !-------------------------------------------------------------------!
  ! Set the set current time step number
  !-------------------------------------------------------------------!  
  
  subroutine set_current_time_step(this, step_num)

    class(system_descriptor) :: this
    integer :: step_num

    this % current_time_step = step_num

  end subroutine set_current_time_step

  !-------------------------------------------------------------------!
  ! Set the set current time
  !-------------------------------------------------------------------!  
  
  subroutine set_current_time(this, time)

    class(system_descriptor) :: this
    real(dp) :: time

    this % current_time = time
    
  end subroutine set_current_time
  
  !-------------------------------------------------------------------!
  ! Set the step size
  !-------------------------------------------------------------------!  
  
  subroutine set_step_size(this, dt)
    
    class(system_descriptor) :: this
    real(dp) :: dt

    this % dt = dt
    
  end subroutine set_step_size

  !-------------------------------------------------------------------!
  ! Set the start time
  !-------------------------------------------------------------------!  
  
  subroutine set_start_time(this, tstart)
    
    class(system_descriptor) :: this
    real(dp) :: tstart

    this % start_time = tstart
    
  end subroutine set_start_time
  
  !-------------------------------------------------------------------!
  ! Set the end time
  !-------------------------------------------------------------------!  
  
  subroutine set_end_time(this, tend)

    class(system_descriptor) :: this
    real(dp) :: tend

    this % end_time = tend

  end subroutine set_end_time
  
  !-------------------------------------------------------------------!
  ! Is this an unsteady simulation
  !-------------------------------------------------------------------!  
  
  subroutine set_unsteady(this, unsteady)

    class(system_descriptor) :: this
    logical :: unsteady

    this % unsteady = unsteady

  end subroutine set_unsteady

  !-------------------------------------------------------------------!
  ! Get the number of variables
  !-------------------------------------------------------------------!  
  
  integer pure function get_num_vars(this)
    
    class(system_descriptor), intent(in) :: this

    get_num_vars = this % nvars
    
  end function get_num_vars
  
  !-------------------------------------------------------------------!
  ! Get the number of time steps
  !-------------------------------------------------------------------!  
  
  integer function get_num_time_steps(this)
    
    class(system_descriptor) :: this
    
    get_num_time_steps = this % num_time_steps
    
  end function get_num_time_steps

  !-------------------------------------------------------------------!
  ! Get the current time step
  !-------------------------------------------------------------------!  
  
  integer function get_current_time_step(this)
    
    class(system_descriptor) :: this
    
    get_current_time_step  = this % current_time_step

  end function get_current_time_step

  !-------------------------------------------------------------------!
  ! Get the current time
  !-------------------------------------------------------------------!  
  
  real(dp) function get_current_time(this)
    
    class(system_descriptor) :: this
    
    get_current_time = this % current_time

  end function get_current_time
  
  !-------------------------------------------------------------------!
  ! Get the step size
  !-------------------------------------------------------------------!  
  
  real(dp) function get_step_size(this)

    class(system_descriptor) :: this

    get_step_size = this % dt
    
  end function get_step_size

  !-------------------------------------------------------------------!
  ! Get the start time
  !-------------------------------------------------------------------!  
  
  real(dp) function get_start_time(this)

    class(system_descriptor) :: this

    get_start_time = this % start_time

  end function get_start_time
  
  !-------------------------------------------------------------------!
  ! Get the end time
  !-------------------------------------------------------------------!  
  
  real(dp) function get_end_time(this)

    class(system_descriptor) :: this
    
    get_end_time = this % end_time
    
  end function get_end_time
  
  !-------------------------------------------------------------------!
  ! Is this an unsteady simulation
  !-------------------------------------------------------------------!  
  
  logical function is_unsteady(this)

    class(system_descriptor) :: this

    is_unsteady = this % unsteady

  end function is_unsteady
  
  !-------------------------------------------------------------------!
  ! Initialization tasks for system descriptor variables
  !-------------------------------------------------------------------!
  
  subroutine initialize_system_variables(this)
    
    class(system_descriptor) :: this
    integer :: num_steps

    ! Use namelist system to read these params
    call system % set_num_vars(1)
    call system % set_start_time(0.0_dp)
    call system % set_end_time(1.0_dp)
    call system % set_step_size(1.0e-3_dp)
    call system % set_unsteady(.true.)
   
    num_steps = int((system % get_end_time()-system%get_start_time()) &
         &/system % get_step_size())+1
    
    call system % set_num_time_steps(num_steps)
    
  end subroutine initialize_system_variables
  
  !-------------------------------------------------------------------!
  ! Initialization tasks for mesh variables
  !-------------------------------------------------------------------!
  
  subroutine initialize_mesh_variables(this)

    class(mesh_variables) :: this

    write(*, *) "Initialize mesh variables"
    
  end subroutine initialize_mesh_variables
  
  !-------------------------------------------------------------------!
  ! Set initial condition
  !-------------------------------------------------------------------!
  
  subroutine set_initial_condition(this, qinit, qdotinit)

    class(state_variables) :: this
    real(dp), dimension(:) :: qinit, qdotinit

    write(*,*) "Setting initial condition"

    this % q(1,:) = qinit
    this % qdot(1,:) = qdotinit
    
  end subroutine set_initial_condition

  !-------------------------------------------------------------------!
  ! Initialization tasks for state variables in the simulation
  !-------------------------------------------------------------------!
  
  subroutine initialize_state_variables(this)

    class(state_variables) :: this

    write(*, *) "Initializing state variables"
    
    !-----------------------------------------------------------------!
    ! allocate and initialize q
    !-----------------------------------------------------------------!
    if (allocated(this % q)) deallocate(this % q)
    allocate(this % q(system % num_time_steps, system % nvars))
    this % q = 0.0_dp

    !-----------------------------------------------------------------!
    ! allocate and initialize qdot
    !-----------------------------------------------------------------!
    if (allocated(this % qdot)) deallocate(this % qdot)
    allocate(this % qdot(system % num_time_steps, system % nvars))
    this % qdot = 0.0_dp

    !-----------------------------------------------------------------!
    ! allocate and initialize qddot
    !-----------------------------------------------------------------!
    if (allocated(this % qddot)) deallocate(this % qddot)
    allocate(this % qddot(system % num_time_steps, system % nvars))
    this % qddot = 0.0_dp

    !-----------------------------------------------------------------!
    ! allocate and initialize dq (update)
    !-----------------------------------------------------------------!
    if (allocated(this % dq)) deallocate(this % dq)
    allocate(this % dq(system % nvars))
    this % dq = 0.0_dp

    if (allocated(this % R)) deallocate(this % R)
    allocate(this % R(system % num_time_steps, system % nvars))
    this % R = 0.0_dp

    if (allocated(this % dR)) deallocate(this % dR)
    allocate(this % dR(system % num_time_steps, system % nvars, system % nvars))
    this % dR = 0.0_dp

  end subroutine initialize_state_variables
  
  !-------------------------------------------------------------------!
  ! Finalize the state variables and freeup memory
  !-------------------------------------------------------------------!
  
  subroutine finalize_state_variables(this)

    class(state_variables) :: this

    write(*, *) "Finalize variables"

    if (allocated(this % q))     deallocate(this % q)
    if (allocated(this % qdot))  deallocate(this % qdot)
    if (allocated(this % qddot)) deallocate(this % qddot)

  end subroutine finalize_state_variables

  !-------------------------------------------------------------------!
  ! Finalize the mesh variables and freeup memory
  !-------------------------------------------------------------------!
  
  subroutine finalize_mesh_variables(this)

    class(mesh_variables) :: this

    write(*,*) "Finalize mesh variables"

  end subroutine finalize_mesh_variables

end module system_class

!=====================================================================!
! Module that wraps the backward difference integration logic
!=====================================================================!

module backward_difference

  use precision
  use system_class

  implicit none

  private                                                               

  public get_bdf_coeffs, update_states, extrapolate

contains

  !-------------------------------------------------------------------!
  ! Returns the bdf coeffs (unscaled with respect to the step size h)
  !-------------------------------------------------------------------!
  ! Input:
  !-------------------------------------------------------------------!
  ! d : d-th derivative
  ! m : order of accuracy
  !-------------------------------------------------------------------!
  ! Output:
  !-------------------------------------------------------------------!
  ! Array of size (d+m) containing the coefficients
  !-------------------------------------------------------------------!

  function get_bdf_coeffs(d, m) result(c)

    integer, intent(in)    :: d
    integer, intent(in)    :: m
    integer                :: n, k

    real(dp)               :: c(d+m), ctmp(d+m)
    real(dp)               :: x(d+m)
    real(dp), parameter    :: h = 1.0_dp

    ! number of data points needed for desired derivative degree and
    ! order of accuracy
    n = m + d
    call differ_backward ( h, d, m, ctmp, x )

    ! flip the array
    forall(k = 1:n) c(k) = ctmp(n-k+1)

  end function get_bdf_coeffs

  !-------------------------------------------------------------------!
  ! Update the states at each newton iteration
  !-------------------------------------------------------------------!

  subroutine update_states()

    implicit none

    real(dp) :: alpha(20), beta(20), dt, dt2

    integer  :: current_time_step
    integer  :: order, k

    current_time_step = system % get_current_time_step()
    dt = system % get_step_size()
    dt2 = dt * dt

    !-----------------------------------------------------------------!
    ! Update the state
    !-----------------------------------------------------------------!

    state % q(current_time_step, :) = state % q(current_time_step, :) &
         &+ state % dq(:)

    !-----------------------------------------------------------------!
    ! FD approximation to I derivative
    !-----------------------------------------------------------------!

    ! Order of accuracy for first derivative
    order = current_time_step - 1
    if (order .gt. 3) order = 3

    ! Get the BDF coefficients
    alpha(1:order+1) = get_bdf_coeffs(1, order)

    ! Approximate qdot
    forall(k = 1:order+1) state % qdot(k,:) = state % qdot(k,:) &
         &+ alpha(k)*state % q(k,:)/dt

    !-----------------------------------------------------------------!
    ! FD approximation to II derivative
    !-----------------------------------------------------------------!

    ! Order of accuracy for second derivative
    if (current_time_step .le. 3) then
       order = 1
    else
       order = 2
    end if

    ! Get the BDF coefficients
    beta(1:order+2) = get_bdf_coeffs(2, order)

    ! Approximate qddot
    forall(k = 1:order+1) state % qddot(k,:) = state % qddot(k,:) &
         &+ beta(k)*state % q(k,:)/dt2

  end subroutine update_states

  !-------------------------------------------------------------------!
  ! Extrapolate to next time step
  !-------------------------------------------------------------------!

  subroutine extrapolate()

    integer  :: step_num
    real(dp) :: dt, dt2

    dt = system % get_step_size()
    dt2 = dt * dt

    step_num = system % get_current_time_step()

    if (step_num .gt. 1) then
       state % q(step_num,:) = state % q(step_num-1,:) &
            &+ state % qdot(step_num-1,:)*dt &
            &+ state % qddot(step_num-1,:)*dt2/2.0_dp
    end if
    
  end subroutine extrapolate

!!$
!!$  subroutine integrate(this)
!!$
!!$    type(state_variables) :: this
!!$
!!$    integer :: k
!!$
!!$    ! Set optional parameters
!!$    ! call newton % set_num_vars(nvars)
!!$    ! call newton % set_exit_on_failure(.true.)
!!$
!!$    ! set initial condition
!!$    this % q(1, :) = 1.0_dp
!!$
!!$    ! march in time
!!$    time: do k = 2, this % num_time_steps + 1
!!$
!!$       !call newton % solve()
!!$
!!$    end do time
!!$
!!$  end subroutine integrate

end module backward_difference

!=====================================================================!
! The user should provide implementation for the residual subroutine 
!=====================================================================!

module residual_class

  use precision
  use system_class

  implicit none

  type, extends(abstract_residual) :: residual

     real(dp) :: M = 1.0_dp
     real(dp) :: C = 0.1_dp
     real(dp) :: K = 5.0_dp

   contains

     procedure :: get_residual => user_residual

  end type residual

contains
  
  subroutine user_residual(this)
    
    class(residual) :: this
    integer :: tstep

    tstep = system % get_current_time_step()

    state % R(tstep,:) = this % M * state%qddot(tstep,:) &
         &+ this % C * state % qdot(tstep,:) &
         &+ this % K * state % q(tstep,:)


  end subroutine user_residual

end module residual_class

!=====================================================================!
! The user should provide implementation for the jacobian subroutine 
!=====================================================================!

module jacobian_class

  use precision
  use system_class

  implicit none
  
  type, extends(abstract_jacobian) :: jacobian_matrix

     real(dp) :: M = 1.0_dp
     real(dp) :: C = 0.1_dp
     real(dp) :: K = 5.0_dp

   contains

     procedure :: get_jacobian => user_jacobian

  end type jacobian_matrix

contains

  subroutine user_jacobian(this)

    !use backward_difference

    class(jacobian_matrix) :: this

    real(dp) :: dt, dt2
    integer :: tstep

    tstep = system % get_current_time_step()
   
    dt = system % get_step_size()
    dt2 = dt * dt

    state % dR(tstep,:,:) = this % K + this % C/dt + this % M/dt2

  end subroutine user_jacobian

end module jacobian_class

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
! call newton % init()
! call newton % solve()
!
! Author :  Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module newton_solve_class

  ! import dependencies
  use precision
  use system_class

  use residual_class
  use jacobian_class

  use backward_difference
  use newton_solve_bean_class

  ! no implicit varaible definitions
  implicit none

  ! all routines and varaibles are private by default
  private

  ! expose datatypes
  public :: newton_solve

  ! a type that contains the logic for Newton's method
  type, extends(newton_solve_bean) :: newton_solve

   contains

     ! public procedures
     procedure, public :: solve => solve

     ! private procedures
     procedure, private :: init => init
     procedure, private :: finish => finish
     procedure, private :: work => work

     procedure, private :: linear_solve
     procedure, private :: check_stop

  end type newton_solve

contains

  !-------------------------------------------------------------------!
  ! Solve the linear system to find the newton update
  !-------------------------------------------------------------------!

  subroutine linear_solve(this)

    class(newton_solve) :: this

    type(residual) :: res
    type(jacobian_matrix) :: jac
   
    real(dp):: rtmp, drtmp
    integer :: step_num

    step_num = system % get_current_time_step()

    if (system % get_num_vars() .eq. 1) then

       rtmp  = state % R(step_num, 1)
       drtmp = state % dR(step_num, 1, 1)
       
       state % dq = -rtmp/drtmp
       
    else

       !   state % dq    = - res % R(step_num,:)
       !  / jac % dR (step_num, :, :)
       stop"linear solve not implemented for ndim>1"
       !  lapack ()

    end if

  end subroutine linear_solve

  !-------------------------------------------------------------------!
  ! Solve the linear system to find the newton update
  !-------------------------------------------------------------------!

  subroutine check_stop(this)

    class(newton_solve) :: this

    type(residual) :: res
    integer :: step_num
    
    step_num =  system % get_current_time_step()

    this % rnrm (this % iter_num) = norm2(state % R(step_num,:))
    this % unrm (this % iter_num) = norm2(state % dq)
    
    if ((this % rnrm (this % iter_num).le.this % get_atol_rnrm()) .or.&
         & (this % unrm(this % iter_num).le.this % get_atol_unrm()))  &
         & this % converged = .true.
    
  end subroutine check_stop

  !-------------------------------------------------------------------!
  ! Routine for initialization tasks
  !-------------------------------------------------------------------!
  
  subroutine init(this)

    class(newton_solve) :: this

    print *, "Initializing Newton Solve"

    if (allocated(this % rnrm)) deallocate(this % rnrm)
    allocate(this % rnrm(this % get_max_newton_iters()))
    this%rnrm = 0.0_dp

    if (allocated(this % unrm)) deallocate(this % unrm)
    allocate(this % unrm(this % get_max_newton_iters()))
    this%unrm = 0.0_dp

  end subroutine init
  
  !-------------------------------------------------------------------!
  ! Routine that wraps the logic of newton solve                                                          
  !-------------------------------------------------------------------!
  
  subroutine work(this)

    class(newton_solve) :: this

    type(jacobian_matrix) :: jac
    type(residual) :: res

    integer :: n, step_num

    print *, "Executing Newton solve"

    call extrapolate()

    newton: do n = 1, this % get_max_newton_iters()

       this % iter_num = this % iter_num + 1

      call res % get_residual()
      
      !      print *, state % R(2,:)

      call jac % get_jacobian()
      
      !     print *, state % dR(2,:,:)
      
      call this % linear_solve()
      
      call update_states()
      
      print*, state % q(2,:), state % qdot(2,:), state % qddot(2,:)
      
      call this % check_stop()
      
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
! Test program for the above module
!=====================================================================!

program test

  use precision
  use system_class
  use newton_solve_class
  
  real(dp) :: qinit(1), qdotinit(1)
  type(newton_solve) :: newton

  qinit(1) = 1.0_dp
  qdotinit(1) = 0.0_dp

  call initialize_simulation()
  call state % set_initial_condition(qinit, qdotinit)
  call system % set_current_time_step(2)
  call newton % solve()

contains

end program test





