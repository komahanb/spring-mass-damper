!=====================================================================!
! Module containing all the interfaces for the system
!=====================================================================!

module interfaces

  use precision

  implicit none

  private

  public :: abstract_unsteady_descriptor
  public :: abstract_state_variables
  public :: abstract_mesh_variables
  public :: abstract_residual
  public :: abstract_jacobian

  !-------------------------------------------------------------------!
  ! The abstract unsteady system descriptor
  !-------------------------------------------------------------------!

  type, abstract :: abstract_unsteady_descriptor

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

     procedure :: get_num_vars, set_num_vars
     procedure :: get_num_time_steps, set_num_time_steps
     procedure :: get_current_time_step, set_current_time_step
     procedure :: get_current_time, set_current_time

     procedure :: get_step_size, set_step_size
     procedure :: get_start_time, set_start_time
     procedure :: get_end_time, set_end_time

     procedure :: set_unsteady, is_unsteady

  end type abstract_unsteady_descriptor

  !-------------------------------------------------------------------!
  ! The abstract state variables
  !-------------------------------------------------------------------!
  
  type, abstract :: abstract_state_variables

     real(dp), dimension(:,:), allocatable :: q
     real(dp), dimension(:,:), allocatable :: qdot
     real(dp), dimension(:,:), allocatable :: qddot

     real(dp), dimension(:),   allocatable :: dq

     real(dp), dimension(:,:), allocatable :: R
     real(dp), dimension(:,:,:), allocatable :: dR

   contains

     procedure :: set_initial_state
     procedure :: set_init_q, set_init_qdot

  end type abstract_state_variables

  !-------------------------------------------------------------------!
  ! The abstract  mesh variables
  !-------------------------------------------------------------------!
  
  type, abstract :: abstract_mesh_variables
     !   contains
     !     procedure :: initialize_mesh_variables
     !     procedure :: finalize_mesh_variables
  end type abstract_mesh_variables
  
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

! All the default implementations go next
contains

  !-------------------------------------------------------------------!
  ! Set the initial value of x (starting point)
  !-------------------------------------------------------------------!

  subroutine set_init_q(this, init_q)

    class(abstract_state_variables) :: this
    real(dp), dimension(:) :: init_q

    this % q(1,:) = init_q

  end subroutine set_init_q

  !-------------------------------------------------------------------!
  ! Set the initial value of x (starting point)
  !-------------------------------------------------------------------!

  subroutine set_init_qdot(this, init_qdot)

    class(abstract_state_variables) :: this
    real(dp), dimension(:) :: init_qdot

    this % qdot(1,:) = init_qdot

  end subroutine set_init_qdot

  !-------------------------------------------------------------------!
  ! Set the number of variables
  !-------------------------------------------------------------------!  

  subroutine set_num_vars(this, nvars)

    class(abstract_unsteady_descriptor) :: this
    integer :: nvars

    this % nvars = nvars

  end subroutine set_num_vars

  !-------------------------------------------------------------------!
  ! Set the number of time steps
  !-------------------------------------------------------------------!  

  subroutine set_num_time_steps(this, nsteps)

    class(abstract_unsteady_descriptor) :: this
    integer :: nsteps

    this % num_time_steps = nsteps

  end subroutine set_num_time_steps

  !-------------------------------------------------------------------!
  ! Set the set current time step number
  !-------------------------------------------------------------------!  

  subroutine set_current_time_step(this, step_num)

    class(abstract_unsteady_descriptor) :: this
    integer :: step_num

    this % current_time_step = step_num

  end subroutine set_current_time_step

  !-------------------------------------------------------------------!
  ! Set the set current time
  !-------------------------------------------------------------------!  

  subroutine set_current_time(this, time)

    class(abstract_unsteady_descriptor) :: this
    real(dp) :: time

    this % current_time = time

  end subroutine set_current_time

  !-------------------------------------------------------------------!
  ! Set the step size
  !-------------------------------------------------------------------!  

  subroutine set_step_size(this, dt)

    class(abstract_unsteady_descriptor) :: this
    real(dp) :: dt

    this % dt = dt

  end subroutine set_step_size

  !-------------------------------------------------------------------!
  ! Set the start time
  !-------------------------------------------------------------------!  

  subroutine set_start_time(this, tstart)

    class(abstract_unsteady_descriptor) :: this
    real(dp) :: tstart

    this % start_time = tstart

  end subroutine set_start_time

  !-------------------------------------------------------------------!
  ! Set the end time
  !-------------------------------------------------------------------!  

  subroutine set_end_time(this, tend)

    class(abstract_unsteady_descriptor) :: this
    real(dp) :: tend

    this % end_time = tend

  end subroutine set_end_time

  !-------------------------------------------------------------------!
  ! Is this an unsteady simulation
  !-------------------------------------------------------------------!  

  subroutine set_unsteady(this, unsteady)

    class(abstract_unsteady_descriptor) :: this
    logical :: unsteady

    this % unsteady = unsteady

  end subroutine set_unsteady

  !-------------------------------------------------------------------!
  ! Get the number of variables
  !-------------------------------------------------------------------!  

  integer pure function get_num_vars(this)

    class(abstract_unsteady_descriptor), intent(in) :: this

    get_num_vars = this % nvars

  end function get_num_vars

  !-------------------------------------------------------------------!
  ! Get the number of time steps
  !-------------------------------------------------------------------!  

  integer function get_num_time_steps(this)

    class(abstract_unsteady_descriptor) :: this

    get_num_time_steps = this % num_time_steps

  end function get_num_time_steps

  !-------------------------------------------------------------------!
  ! Get the current time step
  !-------------------------------------------------------------------!  

  integer function get_current_time_step(this)

    class(abstract_unsteady_descriptor) :: this

    get_current_time_step  = this % current_time_step

  end function get_current_time_step

  !-------------------------------------------------------------------!
  ! Get the current time
  !-------------------------------------------------------------------!  

  real(dp) function get_current_time(this)

    class(abstract_unsteady_descriptor) :: this

    get_current_time = this % current_time

  end function get_current_time

  !-------------------------------------------------------------------!
  ! Get the step size
  !-------------------------------------------------------------------!  

  real(dp) function get_step_size(this)

    class(abstract_unsteady_descriptor) :: this

    get_step_size = this % dt

  end function get_step_size

  !-------------------------------------------------------------------!
  ! Get the start time
  !-------------------------------------------------------------------!  

  real(dp) function get_start_time(this)

    class(abstract_unsteady_descriptor) :: this

    get_start_time = this % start_time

  end function get_start_time

  !-------------------------------------------------------------------!
  ! Get the end time
  !-------------------------------------------------------------------!  

  real(dp) function get_end_time(this)

    class(abstract_unsteady_descriptor) :: this

    get_end_time = this % end_time

  end function get_end_time

  !-------------------------------------------------------------------!
  ! Is this an unsteady simulation
  !-------------------------------------------------------------------!  

  logical function is_unsteady(this)

    class(abstract_unsteady_descriptor) :: this

    is_unsteady = this % unsteady

  end function is_unsteady

  !-------------------------------------------------------------------!
  ! Set initial condition
  !-------------------------------------------------------------------!

  subroutine set_initial_state(this, qinit, qdotinit)

    class(abstract_state_variables) :: this
    real(dp), dimension(:) :: qinit, qdotinit

    call this % set_init_q(qinit)
    call this % set_init_qdot(qdotinit)

  end subroutine set_initial_state

end module interfaces
