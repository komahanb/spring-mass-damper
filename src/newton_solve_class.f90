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
  
  !-------------------------------------------------------------------!
  ! The system parameters 
  !-------------------------------------------------------------------!
  
  real(dp) :: M = 1.0_dp
  real(dp) :: C = 0.02_dp
  real(dp) :: K = 5.0_dp
  
  !  real(dp), dimension(:)     , allocatable :: dq
  
  type system_descriptor

     integer  :: ndim
     integer  :: step_num = 1
     integer  :: num_time_steps

     real(dp) :: dt, dt2
     real(dp) :: tstart
     real(dp) :: tend

     logical  :: unsteady = .false.

  end type system_descriptor

  !-------------------------------------------------------------------!
  ! The state variables
  !-------------------------------------------------------------------!
  type state_variables

     type(system_descriptor) :: bean

     real(dp), dimension(:,:), allocatable :: q
     real(dp), dimension(:,:), allocatable :: qdot
     real(dp), dimension(:,:), allocatable :: qddot

  end type state_variables

  type mesh_variables

     type(system_descriptor) :: bean

     ! any mesh based system
  end type mesh_variables
  
  ! type that wraps all the variables
  type blob 
     type(system_descriptor) :: bean
     type(state_variables)   :: state
     type(mesh_variables)    :: mesh
  end type blob

  type residual_vector
     type(system_descriptor) :: bean
     real(dp), dimension(:,:), allocatable :: R 
  end type residual_vector

  type jacobian_matrix
     type(system_descriptor) :: bean
     real(dp), dimension(:,:,:), allocatable :: dR
  end type jacobian_matrix

contains

  ! initialization tasks for the variables in simulation
  subroutine initialize()

!!$    integer, intent(in) :: ndimin, nstepsin
!!$
!!$    write(*, *) "Initializing variables"
!!$
!!$    ! set the number of dimensions or variables
!!$    ndim = ndimin
!!$    if (ndim .eq. 0) stop "Wrong dimension"
!!$
!!$    ! set the number of time steps
!!$    nsteps = nstepsin
!!$    if (nsteps .le. 0) stop "Wrong number of time steps"
!!$
!!$    ! allocate and initialize q
!!$    if (allocated(q)) deallocate(q)
!!$    allocate(q(nsteps, ndim))
!!$    q = 0.0_dp
!!$
!!$    ! allocate and initialize qdot
!!$    if (allocated(qdot)) deallocate(qdot)
!!$    allocate(qdot(nsteps, ndim))
!!$    qdot = 0.0_dp
!!$
!!$    ! allocate and initialize qddot
!!$    if (allocated(qddot)) deallocate(qddot)
!!$    allocate(qddot(nsteps, ndim))
!!$    qddot = 0.0_dp
!!$
!!$    ! allocate and initialize R (residual)
!!$    if (allocated(R)) deallocate(R)
!!$    allocate(R(nsteps, ndim))
!!$    R = 0.0_dp
!!$
!!$    ! allocate and initialize dq (update)
!!$    if (allocated(dq)) deallocate(dq)
!!$    allocate(dq(ndim))
!!$    dq = 0.0_dp
!!$
!!$    ! allocate and initialize dR (jacobian)
!!$    if (allocated(dR)) deallocate(dR)
!!$    allocate(dR(nsteps, ndim, ndim))
!!$    dR = 0.0_dp
!!$
  end subroutine initialize
!!$
!!$  ! finalization tasks
  subroutine finalize()
!!$
!!$    write(*, *) "Finalize variables"
!!$
!!$    if (allocated(q))     deallocate(q)
!!$    if (allocated(qdot))  deallocate(qdot)
!!$    if (allocated(qddot)) deallocate(qddot)
!!$    if (allocated(R))     deallocate(R)
!!$    if (allocated(dR))    deallocate(dR)
!!$    if (allocated(dq))    deallocate(dq)

  end subroutine finalize

end module system_class

!=====================================================================!
! Module that wraps the backward difference integration logic
!=====================================================================!

module backward_difference

  use precision
  use system_class

  implicit none

  private                                                               
  
  public get_bdf_coeffs, update_states, residual, jacobian

  public system_descriptor
  
contains
 
  ! implementation for residual
  subroutine residual(this)
    
    type(system_descriptor) :: this

    real(dp) :: dt, dt2
    integer  :: step_num

    step_num =  this % step_num
    dt       = this % dt
    dt2      = this % dt2

    !R(step_num, :) = M*qddot(step_num,:) + C*qdot(step_num,:) + K*q(step_num,:)
    
  end subroutine residual
  
  ! implementation for the jacobian
  subroutine jacobian(this)
    
    type(system_descriptor) :: this

    real(dp) :: dt, dt2
    integer  :: step_num

    step_num =  this % step_num
    dt       = this % dt
    dt2      = this % dt2
    
    !dR(step_num, :, :) = K + C/dt + M/dt2

  end subroutine jacobian

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

  subroutine initialize1(this, ndim, tstart, tend, num_time_steps)

    type(system_descriptor)  :: this
    integer  :: ndim
    real(dp) :: tstart, tend
    integer  :: num_time_steps

    this % ndim = ndim
    this % tstart = tstart
    this % tend   = tend
    this % num_time_steps = num_time_steps
    this % dt = (tend - tstart)/dble(num_time_steps)
    this % dt2 = this % dt * this % dt

    if (this % tstart .ne. this %tend) this % unsteady = .true.

  end subroutine initialize1

  subroutine initialize2(this, ndim, tstart, tend, dt)

    type(system_descriptor)  :: this
    integer  :: ndim
    real(dp) :: tstart, tend, dt

    this % ndim = ndim
    this % tstart = tstart
    this % tend   = tend
    this % dt = dt
    this % dt2 = dt*dt
    this % num_time_steps = int(dt*(tend-tstart))

    if (this % tstart .ne. this %tend) this % unsteady = .true.
    
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
  
  subroutine update_states(this, step_num)

    type(system_descriptor) :: this
    
    real(dp) :: alpha(20), beta(20), dt, dt2

    integer  :: step_num
    integer  :: order, k

    step_num =  this % step_num
    dt       = this % dt
    dt2      = this % dt2
    !-----------------------------------------------------------------!
    ! Update the state
    !-----------------------------------------------------------------!
    
    !q(step_num, :) = q(step_num, :) + dq(:)

    !-----------------------------------------------------------------!
    ! FD approximation to I derivative
    !-----------------------------------------------------------------!

    ! Order of accuracy for first derivative
    order = step_num - 1
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
    if (step_num .le. 3) then
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

    !print*, dq,  q(step_num, :), qdot(step_num, :),qddot(step_num, :)

  end subroutine update_states

end module backward_difference

program test
  print*, "Hello world!!!"
end program test
