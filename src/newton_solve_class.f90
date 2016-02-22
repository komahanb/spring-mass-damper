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
  public :: system, state, mesh, initialize_simulation, set_initial_state

  type :: system_descriptor

     private 

     integer  :: nvars = 1
     integer  :: num_time_steps  = 1
     integer  :: current_time_step = 1

     real(dp) :: start_time = 0.0_dp
     real(dp) :: end_time = 1.0_dp
     real(dp) :: current_time = 0.0_dp
     real(dp) :: dt = 1.0e-1_dp

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

     procedure :: set_initial_state
     procedure :: set_init_q, set_init_qdot

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
  ! Set the initial value of x (starting point)
  !-------------------------------------------------------------------!

  subroutine set_init_q(this, init_q)

    class(state_variables) :: this
    real(dp), dimension(:) :: init_q

    this % q(1,:) = init_q

  end subroutine set_init_q

  !-------------------------------------------------------------------!
  ! Set the initial value of x (starting point)
  !-------------------------------------------------------------------!

  subroutine set_init_qdot(this, init_qdot)

    class(state_variables) :: this
    real(dp), dimension(:) :: init_qdot

    this % qdot(1,:) = init_qdot

  end subroutine set_init_qdot

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
    call system % set_end_time(10.0_dp)
    call system % set_step_size(1.0e-3_dp)
    call system % set_unsteady(.true.)

    num_steps = int((system % get_end_time()-system%get_start_time()) &
         &/system % get_step_size()) + 1

    call system % set_num_time_steps(num_steps)

  end subroutine initialize_system_variables

  !-------------------------------------------------------------------!
  ! Initialization tasks for mesh variables
  !-------------------------------------------------------------------!

  subroutine initialize_mesh_variables(this)

    class(mesh_variables) :: this

  end subroutine initialize_mesh_variables

  !-------------------------------------------------------------------!
  ! Set initial condition
  !-------------------------------------------------------------------!

  subroutine set_initial_state(this, qinit, qdotinit)

    class(state_variables) :: this
    real(dp), dimension(:) :: qinit, qdotinit

    call this % set_init_q(qinit)
    call this % set_init_qdot(qdotinit)

  end subroutine set_initial_state
  
  !-------------------------------------------------------------------!
  ! Initialization tasks for state variables in the simulation
  !-------------------------------------------------------------------!

  subroutine initialize_state_variables(this)

    class(state_variables) :: this

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

    !-----------------------------------------------------------------!
    ! allocate and initialize R
    !-----------------------------------------------------------------!
    if (allocated(this % R)) deallocate(this % R)
    allocate(this % R(system % num_time_steps, system % nvars))
    this % R = 0.0_dp

    !-----------------------------------------------------------------!
    ! allocate and initialize dR
    !-----------------------------------------------------------------!
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

  public aa, bb
  public get_bdf_coeffs, update_states, extrapolate, approximate_derivatives

  ! Coefficients for linearising the residual

  real(dp) :: alpha(6+1), beta(2*6+1)
  real(dp) :: aa = 0.0_dp, bb = 0.0_dp

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
    call differ_backward ( h, d, m, ctmp, x )

    n = m + d

    ! flip the array
    forall(k=1:n) c(k) = ctmp(n-k+1)

    if(d .eq. 1) aa = c(1)
    if(d .eq. 2) bb = c(1)

  end function get_bdf_coeffs

  !-------------------------------------------------------------------!
  ! Updates the state vector and its derivatives at the k-th time step
  !-------------------------------------------------------------------!

  subroutine update_states()

    integer :: k
    real(dp) ::dt, dt2

    k = system % get_current_time_step()

    dt = system % get_step_size()
    dt2 = dt * dt

    state % q(k,:) = state % q(k,:) + state % dq(:) 
    state % qdot(k, :) = state % qdot(k,:) + aa*state % dq(:)/dt
    state % qddot(k, :) = state % qddot(k,:) + bb*state % dq(:)/dt2  

  end subroutine update_states

  !-------------------------------------------------------------------!
  ! Approximate first and second derivative for the next iteration
  !-------------------------------------------------------------------!

  subroutine approximate_derivatives()

    implicit none

    real(dp) :: dt, dt2
    integer  :: forder, sorder, k, i

    k = system % get_current_time_step()

    dt = system % get_step_size()
    dt2 = dt * dt

    !-----------------------------------------------------------------!
    ! Approximate state at the current step (copy the previous state)
    !-----------------------------------------------------------------!

    ! state % q(k, :) = state % q(k-1, :)

    !-----------------------------------------------------------------!
    ! FD approximation to I derivative
    !-----------------------------------------------------------------!

    ! Determine the order of accuracy for first derivative
    forder = k - 1
    if (forder .gt. 6) forder = 6

    ! Get the BDF coefficients for first derivative
    alpha(1:forder+1) = get_bdf_coeffs(1, forder)

    ! Apply the BDF formula
    do i = 0, forder ! m+1 points
       state % qdot(k,:) = state % qdot(k,:) + &
            & alpha(i+1) * state % q(k-i,:)/dt 
    end do

    !-----------------------------------------------------------------!
    ! FD approximation to II derivative
    !-----------------------------------------------------------------!

    ! Order of accuracy for second derivative
    sorder = (k-1)/2
    if (sorder .gt. 6) sorder = 6

    if (sorder .gt. 0) then

       ! Get the BDF coefficients for second derivative   
       beta(1:2*sorder+1) = get_bdf_coeffs(2, sorder)

       ! Apply the BDF formula for second derivative
       do i = 0, 2*sorder ! 2m+1 points
          state % qddot(k,:) = state % qddot(k,:) + &
               & beta(i+1) * state % q(k-i,:)/dt2
       end do

    else

       !  we dont have enought points yet
       bb = 1.0_dp
       state % qddot(k,:) = (state%qdot(k-1,:) - state%qdot(k,:))/dt

    end if

  end subroutine approximate_derivatives

  !-------------------------------------------------------------------!
  ! Extrapolate to next time step
  !-------------------------------------------------------------------!

  subroutine extrapolate()

    integer  :: k
    real(dp) :: dt, dt2

    k = system % get_current_time_step()

    dt = system % get_step_size()
    dt2 = dt*dt

    ! Extrapolate using gradient and Hessian information
    if (k .gt. 1) then
       state % q(k,:) = state % q(k-1,:) &
            &+ state % qdot(k-1,:)*dt &
            &+ state % qddot(k-1,:)*dt2/2.0_dp
    end if

  end subroutine extrapolate

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
     real(dp) :: C = 0.02_dp
     real(dp) :: K = 5.0_dp

   contains

     procedure :: get_residual => user_residual

  end type residual

contains
  
  subroutine user_residual(this)

    class(residual) :: this
    integer :: n

    n = system % get_current_time_step()

    state % R(n,:) = this % M * state % qddot(n,:) &
         &+ this % C * state % qdot(n,:) &
         &+ this % K * state % q(n,:)

  end subroutine user_residual

end module residual_class

!=====================================================================!
! The user should provide implementation for the jacobian subroutine 
!=====================================================================!

module jacobian_class
  
  use precision
  use system_class
  use backward_difference, only: aa, bb

  implicit none

  type, extends(abstract_jacobian) :: jacobian_matrix

     real(dp) :: M = 1.0_dp
     real(dp) :: C = 0.02_dp
     real(dp) :: K = 5.0_dp

   contains

     procedure :: get_jacobian => user_jacobian

  end type jacobian_matrix

contains
  
  subroutine user_jacobian(this)

    class(jacobian_matrix) :: this

    real(dp) :: dt, dt2
    integer :: n

    n = system % get_current_time_step()

    dt = system % get_step_size()
    dt2 = dt * dt

    state % dR(n,:,:) = this % K + this % C*aa/dt + this % M*bb/dt2

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
     integer :: iter_num = 1

     !----------------------------------------------------------------!
     ! Stopping criteria
     !----------------------------------------------------------------!

     ! absolute tolerances
     real(dp) :: atol_rnrm = 1.0d-12
     real(dp) :: atol_unrm = 1.0d-12

     ! relative tolerances
     real(dp) :: rtol_rnrm = 1.0d-12
     real(dp) :: rtol_unrm = 1.0d-12

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

     procedure :: get_file_number, set_file_number
     procedure :: get_write_newton_details, set_write_newton_details
     procedure :: get_exit_on_failure, set_exit_on_failure

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
  ! Get the output stream filenum
  !-------------------------------------------------------------------!
  
  integer function get_file_number(this)

    class(newton_solve_bean) :: this

    get_file_number = this % filenum

  end function get_file_number

  !-------------------------------------------------------------------!
  ! Should write the newton iteration details at each time step?
  !-------------------------------------------------------------------!
  
  logical function get_write_newton_details(this)

    class(newton_solve_bean) :: this

    get_write_newton_details = this % write_newton_details

  end function get_write_newton_details

  !-------------------------------------------------------------------!
  ! Should write the newton iteration details at each time step?
  !-------------------------------------------------------------------!
  
  logical function get_exit_on_failure(this)
    
    class(newton_solve_bean) :: this

    get_exit_on_failure = this % exit_on_failure

  end function get_exit_on_failure

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
     procedure, private :: write_solution

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
    integer :: k

    k = system % get_current_time_step()

    if (system % get_num_vars() .eq. 1) then

       rtmp  = state % R(k, 1)
       drtmp = state % dR(k, 1, 1)

       state % dq = -rtmp/drtmp

    else

       !   state % dq    = - res % R(k,:)
       !  / jac % dR (k, :, :)
       stop"linear solve not implemented for ndim>1"
       !  lapack ()

    end if

  end subroutine linear_solve
  
  !-------------------------------------------------------------------!
  ! Check if the desired tolerance has been reached
  !-------------------------------------------------------------------!
  
  subroutine check_stop(this)

    class(newton_solve) :: this

    type(residual) :: res
    integer :: k, n

    k =  system % get_current_time_step()
    n = this % iter_num

    this % rnrm(n) = norm2(state % R(k,:))
    this % unrm(n) = norm2(state % dq)

    if ((this % rnrm(n).le.this % get_atol_rnrm()) .or.&
         & (this % unrm(n).le.this % get_atol_unrm()))  &
         & this % converged = .true.

  end subroutine check_stop

  !-------------------------------------------------------------------!
  ! Write the solution to output file or screen
  !-------------------------------------------------------------------!
  
  subroutine write_solution(this)
    
    class(newton_solve) :: this
    integer :: k, n
    
    k = system % get_current_time_step()
    
    write(this% get_file_number(),*) dble(k)*system%get_step_size(),&
         &this % iter_num, state % q(k,:),  state % qdot(k,:),&
         &state % qddot(k,:),  state % R(k,:), state % dR(k,:,:),&
         &state % dq
    
  end subroutine write_solution

  !-------------------------------------------------------------------!
  ! Routine for initialization tasks
  !-------------------------------------------------------------------!

  subroutine init(this)

    class(newton_solve) :: this

    this % converged = .false.

    this % iter_num = 1

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

    integer :: n=0, k=0

    k = system % get_current_time_step()
    
    call extrapolate()
    
    call approximate_derivatives()

    newton: do n = 1, this % get_max_newton_iters()

       call res % get_residual()

       call jac % get_jacobian()

       call this % linear_solve()

       call this % write_solution()
       
       call this % check_stop()

       if (this % converged) then
          exit newton
       else
          call update_states()
       end if

       this % iter_num = this % iter_num + 1

    end do newton

  end subroutine work
  
  !-------------------------------------------------------------------!
  ! Routine that performs finishing tasks
  !-------------------------------------------------------------------!
  
  subroutine finish(this)

    class(newton_solve) :: this

    if (allocated(this%rnrm)) deallocate(this%rnrm)
    if (allocated(this%unrm)) deallocate(this%unrm)
    
    ! close the file
    if(this% get_file_number() .ne. 6) then
       close(this% get_file_number())
    end if
    
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
  use backward_difference

  real(dp) :: qinit(1), qdotinit(1)
  type(newton_solve) :: newton

  qinit(1) = 1.0_dp
  qdotinit(1) = 0.0_dp

  call initialize_simulation()
  call state % set_initial_state(qinit, qdotinit)

  !-------------------------------------------------------------------!
  ! March in time
  !-------------------------------------------------------------------!

  time: do k = 2, system % get_num_time_steps()

     call system % set_current_time_step(k)
     call newton % solve()

     if (.not.newton % converged) exit time

  end do time

contains

end program test
