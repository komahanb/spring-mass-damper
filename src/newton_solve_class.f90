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
  
  type system_descriptor

     integer  :: nvars
     integer  :: num_time_steps
     integer  :: current_time_step

     real(dp) :: start_time
     real(dp) :: end_time
     real(dp) :: current_time
     real(dp) :: dt

     logical  :: unsteady = .false.
     
   contains
     
     ! type bound getters and setters

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
  type, extends(system_descriptor) :: state_variables

     real(dp), dimension(:,:), allocatable :: q
     real(dp), dimension(:,:), allocatable :: qdot
     real(dp), dimension(:,:), allocatable :: qddot

   contains
     
     procedure :: initialize_state_variables
     procedure :: finalize_state_variables

!    procedure :: set_initial_condition
     
  end type state_variables
  
  !-------------------------------------------------------------------!
  ! The mesh variables
  !-------------------------------------------------------------------!
  
  type, extends(state_variables) :: mesh_variables

   contains

     procedure :: initialize_mesh_variables

  end type mesh_variables
  
  !-------------------------------------------------------------------!
  ! Abstract type for residual vector and its interface
  !-------------------------------------------------------------------!
  
  type, abstract :: abstract_residual

     real(dp), dimension(:,:), allocatable :: R 
     
   contains

     procedure(get_residual_interface), deferred :: get_residual

  end type abstract_residual

  ! Interface for implementing the residual

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

     real(dp), dimension(:,:,:), allocatable :: dR

   contains

     procedure(get_jacobian_interface), deferred :: get_jacobian
     
  end type abstract_jacobian

  ! Interface for implementing the jacobian

  abstract interface

     subroutine get_jacobian_interface(this)

       import abstract_jacobian

       class(abstract_jacobian) :: this
       
     end subroutine get_jacobian_interface

  end interface

  !-------------------------------------------------------------------!

  !-------------------------------------------------------------------!

  type(system_descriptor) :: system
  type(state_variables)   :: state
  type(mesh_variables)    :: mesh
  
contains

  !===================================================================!
  ! Set the number of variables
  !===================================================================!  

  subroutine set_num_vars(this, nvars)

    class(system_descriptor) :: this
    integer :: nvars

    this % nvars = nvars
    
  end subroutine set_num_vars

  !===================================================================!
  ! Set the number of time steps
  !===================================================================!  

  subroutine set_num_time_steps(this, nsteps)

    class(system_descriptor) :: this
    integer :: nsteps

    this % num_time_steps = nsteps
    
  end subroutine set_num_time_steps  

  !===================================================================!
  ! Set the set current time step number
  !===================================================================!  
  
  subroutine set_current_time_step(this, step_num)

    class(system_descriptor) :: this
    integer :: step_num

    this % current_time_step = step_num

  end subroutine set_current_time_step

  !===================================================================!
  ! Set the set current time
  !===================================================================!  
  
  subroutine set_current_time(this, time)

    class(system_descriptor) :: this
    real(dp) :: time

    this % current_time = time
    
  end subroutine set_current_time
  
  !===================================================================!
  ! Set the step size
  !===================================================================!  
  
  subroutine set_step_size(this, dt)
    
    class(system_descriptor) :: this
    real(dp) :: dt

    this % dt = dt
    
  end subroutine set_step_size

  !===================================================================!
  ! Set the start time
  !===================================================================!  
  
  subroutine set_start_time(this, tstart)
    
    class(system_descriptor) :: this
    real(dp) :: tstart

    this % start_time = tstart
    
  end subroutine set_start_time
  
  !===================================================================!
  ! Set the end time
  !===================================================================!  
  
  subroutine set_end_time(this, tend)

    class(system_descriptor) :: this
    real(dp) :: tend

    this % end_time = tend

  end subroutine set_end_time
  
  !===================================================================!
  ! Is this an unsteady simulation
  !===================================================================!  
  
  subroutine set_unsteady(this, unsteady)

    class(system_descriptor) :: this
    logical :: unsteady

    this % unsteady = unsteady

  end subroutine set_unsteady

  !===================================================================!
  ! Get the number of variables
  !===================================================================!  
  
  integer function get_num_vars(this)
    
    class(system_descriptor) :: this

    get_num_vars = this % nvars
    
  end function get_num_vars
  
  !===================================================================!
  ! Get the number of time steps
  !===================================================================!  
  
  integer function get_num_time_steps(this)
    
    class(system_descriptor) :: this
    
    get_num_time_steps = this % num_time_steps
    
  end function get_num_time_steps

  !===================================================================!
  ! Get the current time step
  !===================================================================!  
  
  integer function get_current_time_step(this)
    
    class(system_descriptor) :: this
    
    get_current_time_step  = this % current_time_step

  end function get_current_time_step

  !===================================================================!
  ! Get the current time
  !===================================================================!  
  
  real(dp) function get_current_time(this)
    
    class(system_descriptor) :: this
    
    get_current_time = this % current_time

  end function get_current_time
  
  !===================================================================!
  ! Get the step size
  !===================================================================!  
  
  real(dp) function get_step_size(this)

    class(system_descriptor) :: this

    get_step_size = this % dt
    
  end function get_step_size

  !===================================================================!
  ! Get the start time
  !===================================================================!  
  
  real(dp) function get_start_time(this)

    class(system_descriptor) :: this

    get_start_time = this % start_time

  end function get_start_time
  
  !===================================================================!
  ! Get the end time
  !===================================================================!  
  
  real(dp) function get_end_time(this)

    class(system_descriptor) :: this
    
    get_end_time = this % end_time
    
  end function get_end_time
  
  !===================================================================!
  ! Is this an unsteady simulation
  !===================================================================!  
  
  logical function is_unsteady(this)

    class(system_descriptor) :: this

    is_unsteady = this % unsteady

  end function is_unsteady

  !===================================================================!
  ! Initialization tasks for system descriptor variables
  !===================================================================!
  
  subroutine initialize_system_variables(this)

    class(system_descriptor) :: this

    call system % set_num_vars(10)
    call system % set_num_time_steps(100)
    call system % set_start_time(0.0_dp)
    call system % set_end_time(1.0_dp)
    call system % set_unsteady(.true.)

  end subroutine initialize_system_variables
  
  !===================================================================!
  ! Initialization tasks for mesh variables
  !===================================================================!
  
  subroutine initialize_mesh_variables(this)

    class(mesh_variables) :: this
    
  end subroutine initialize_mesh_variables
  
  !===================================================================!
  ! Initialization tasks for ALL  variables in the simulation
  !===================================================================!

  subroutine initialize_simulation()

    call system % initialize_system_variables()
    call state  % initialize_state_variables()
    call mesh   % initialize_mesh_variables()

  end subroutine initialize_simulation
  
  !===================================================================!
  ! Initialization tasks for state variables in the simulation
  !===================================================================!
  
  subroutine initialize_state_variables(this)

    class(state_variables) :: this

    write(*, *) "Initializing state variables"

    ! allocate and initialize q
    if (allocated(this % q)) deallocate(this % q)
    allocate(this % q(this % num_time_steps, this % nvars))

    ! allocate and initialize qdot
    if (allocated(this % qdot)) deallocate(this % qdot)
    allocate(this % qdot(this % num_time_steps, this % nvars))

    ! allocate and initialize qddot
    if (allocated(this % qddot)) deallocate(this % qddot)
    allocate(this % qddot(this % num_time_steps, this % nvars))

  end subroutine initialize_state_variables
  
  !===================================================================!
  ! Finalize the state variables and freeup memory
  !===================================================================!
  
  subroutine finalize_state_variables(this)

    class(state_variables) :: this

    write(*, *) "Finalize variables"

    if (allocated(this % q))     deallocate(this % q)
    if (allocated(this % qdot))  deallocate(this % qdot)
    if (allocated(this % qddot)) deallocate(this % qddot)

  end subroutine finalize_state_variables

end module system_class

!=====================================================================!
! The user should provide implementation for the residual subroutine 
!=====================================================================!

module residual_class

  use system_class

  type, extends(abstract_residual) :: residual

contains

  procedure :: get_residual => user_residual

end type residual

contains

subroutine user_residual(this)

 class(residual) :: this

 stop"please impl residual"

end subroutine user_residual

end module residual_class

!=====================================================================!
! The user should provide implementation for the jacobian subroutine 
!=====================================================================!

module jacobian_class

  use system_class
  
  type, extends(abstract_jacobian) :: jacobian_matrix

   contains

     procedure :: get_jacobian => user_jacobian

  end type jacobian_matrix

contains

subroutine user_jacobian(this)

 class(jacobian_matrix) :: this

 stop"please impl jacobian"

end subroutine user_jacobian

end module jacobian_class

!=====================================================================!
! Module that wraps the backward difference integration logic
!=====================================================================!

module backward_difference

  use precision
  use system_class

  implicit none

  private                                                               
  
  public get_bdf_coeffs, update_states

contains
 
  subroutine integrate(this)

    type(state_variables) :: this

    integer :: k
    
    ! Set optional parameters
    ! call newton % set_num_vars(nvars)
    ! call newton % set_exit_on_failure(.true.)
    
    ! set initial condition
    this % q(1, :) = 1.0_dp
    
    ! march in time
    time: do k = 2, this % num_time_steps + 1
       
       !call newton % solve()
       
    end do time

  end subroutine integrate

  subroutine initialize1(this, nvars, start_time, end_time, num_time_steps)

    type(system_descriptor)  :: this
    integer  :: nvars
    real(dp) :: start_time, end_time
    integer  :: num_time_steps

    this % nvars = nvars
    this % start_time = start_time
    this % end_time   = end_time
    this % num_time_steps = num_time_steps
    this % dt = (end_time - start_time)/dble(num_time_steps)

    if (this % start_time .ne. this %end_time) this % unsteady = .true.

  end subroutine initialize1

  subroutine initialize2(this, nvars, start_time, end_time, dt)

    type(system_descriptor)  :: this
    integer  :: nvars
    real(dp) :: start_time, end_time, dt

    this % nvars = nvars
    this % start_time = start_time
    this % end_time   = end_time
    this % dt = dt
    this % num_time_steps = int(dt*(end_time-start_time))

    if (this % start_time .ne. this %end_time) this % unsteady = .true.
    
  end subroutine initialize2

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
    do k = 1, n
       c(k) = ctmp(n - k + 1)
    end do
    
  end function get_bdf_coeffs

  !-------------------------------------------------------------------!
  ! Update the states at each newton iteration
  !-------------------------------------------------------------------!
  
  subroutine update_states(this, current_time_step)

    type(system_descriptor) :: this
    
    real(dp) :: alpha(20), beta(20), dt, dt2

    integer  :: current_time_step
    integer  :: order, k

    current_time_step =  this % current_time_step
    dt       = this % dt

    !-----------------------------------------------------------------!
    ! Update the state
    !-----------------------------------------------------------------!
    
    !q(current_time_step, :) = q(current_time_step, :) + dq(:)

    !-----------------------------------------------------------------!
    ! FD approximation to I derivative
    !-----------------------------------------------------------------!

    ! Order of accuracy for first derivative
    order = current_time_step - 1
    if (order .gt. 3) order = 3

    ! Get the BDF coefficients
    alpha(1:order+1) = get_bdf_coeffs(1, order)

    ! Find qdot
    do k = 1, order + 1
     !  qdot(k,:) = qdot(k,:) + alpha(k)*q(k,:)/dt
    end do

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

    ! Find qddot
    do k = 1, order + 1
      ! qddot(k,:) = qddot(k,:) + beta(k)*q(k,:)/dt2
    end do

    !print*, dq,  q(current_time_step, :), qdot(current_time_step, :),qddot(current_time_step, :)

  end subroutine update_states

end module backward_difference

program test

  use system_class
  
!  class(system_class) :: mydata

  contains
!!$
!!$    ! implementation for residual
!!$    subroutine get_residual(this)
!!$
!!$      type(system_descriptor) :: this
!!$
!!$      real(dp) :: dt, dt2
!!$      integer  :: current_time_step
!!$
!!$      current_time_step =  this % current_time_step
!!$      dt       = this % dt
!!$      !R(current_time_step, :) = M*qddot(current_time_step,:) + C*qdot(current_time_step,:) + K*q(current_time_step,:)
!!$
!!$    end subroutine get_residual

!!$
!!$    ! implementation for the jacobian
!!$    subroutine get_jacobian(this)
!!$
!!$      type(system_descriptor) :: this
!!$
!!$      real(dp) :: dt, dt2
!!$      integer  :: current_time_step
!!$
!!$      current_time_step =  this % current_time_step
!!$      dt       = this % dt
!!$
!!$      !dR(current_time_step, :, :) = K + C/dt + M/dt2
!!$
!!$    end subroutine 
end program test

!!$
!!$    ! allocate and initialize qdot
!!$    if (allocated(qdot)) deallocate(qdot)
!!$    allocate(qdot(nsteps, nvars))
!!$    qdot = 0.0_dp
!!$
!!$    ! allocate and initialize qddot
!!$    if (allocated(qddot)) deallocate(qddot)
!!$    allocate(qddot(nsteps, nvars))
!!$    qddot = 0.0_dp
!!$
!!$    ! allocate and initialize R (residual)
!!$    if (allocated(R)) deallocate(R)
!!$    allocate(R(nsteps, nvars))
!!$    R = 0.0_dp
!!$
!!$    ! allocate and initialize dq (update)
!!$    if (allocated(dq)) deallocate(dq)
!!$    allocate(dq(nvars))
!!$    dq = 0.0_dp
!!$
!!$    ! allocate and initialize dR (jacobian)
!!$    if (allocated(dR)) deallocate(dR)
!!$    allocate(dR(nsteps, nvars, nvars))
!!$    dR = 0.0_dp
