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
     
!!$     procedure :: get_num_vars, set_num_vars
!!$     procedure :: get_num_time_steps, set_num_time_steps
!!$     procedure :: get_current_time_step, set_current_time_step
!!$     procedure :: get_current_time, set_current_time
!!$     
!!$     procedure :: get_step_size, set_step_size
!!$     procedure :: get_start_time, set_start_time
!!$     procedure :: get_end_time, set_end_time
!!$     
!!$     procedure :: set_unsteady, is_unsteady

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
     
  end type state_variables

  !-------------------------------------------------------------------!
  ! The mesh variables
  !-------------------------------------------------------------------!

  type, extends(state_variables) :: mesh_variables

  end type mesh_variables
  
  !-------------------------------------------------------------------!
  ! Abstract type for residual vector and its interface
  !-------------------------------------------------------------------!
  
  type, abstract :: abstract_residual

     real(dp), dimension(:,:), allocatable :: R 
     
   contains

     procedure(get_residual_interface), nopass, deferred :: get_residual

  end type abstract_residual

  ! Interface for implementing the residual

  abstract interface

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

     procedure(get_jacobian_interface), nopass, deferred :: get_jacobian
     
  end type abstract_jacobian

  ! Interface for implementing the jacobian

  abstract interface

     subroutine get_jacobian_interface(this)

       import abstract_jacobian

       class(abstract_jacobian) :: this
       
     end subroutine get_jacobian_interface

  end interface
  
contains
  
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

!!$
!!$module residual_class
!!$  
!!$  use system_class
!!$
!!$  type, extends(abstract_residual) :: residual
!!$  
!!$   contains
!!$
!!$     procedure :: get_residual => get_residual
!!$     
!!$  end type residual
!!$
!!$contains
!!$  
!!$  subroutine get_residual()
!!$
!!$  end subroutine get_residual
!!$
!!$! implement the residual here
!!$  
!!$end module residual_class

!=====================================================================!
! Module that wraps the backward difference integration logic
!=====================================================================!

module backward_difference

  use precision
  use system_class !, only: system_descriptor, get_residual, get_jacobian

  implicit none

  private                                                               
  
  public get_bdf_coeffs, update_states

  public system_descriptor
  
contains
 
  subroutine integrate(this)

    type(system_descriptor) :: this
!    type(newton_solve) :: newton
!!$
!!$    integer :: nvars  = 1
!!$    integer :: k
!!$
!!$    ! Set optional parameters
!!$ !   call newton % set_num_vars(nvars)
!!$ !   call newton % set_exit_on_failure(.true.)
!!$
!!$    ! set initial condition
!!$    q(1, :) = 1.0_dp
!!$
!!$    ! march in time
!!$    time: do k = 2, this % num_time_steps + 1
!!$
!!$       call newton % solve()
!!$
!!$    end do time

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
